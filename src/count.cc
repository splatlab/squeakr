/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "clipp.h"
#include "ProgOpts.h"
#include "SqueakrFS.h"
#include "squeakrconfig.h"
#include "gqf_cpp.h"
#include "chunk.h"
#include "kmer.h"
#include "reader.h"
#include "util.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define MAX_NUM_THREADS 64

typedef struct {
	CQF<KeyObject> *local_cqf;
	CQF<KeyObject> *main_cqf;
	uint32_t count {0};
	uint32_t ksize {28};
	bool exact;
	spdlog::logger* console{nullptr};
} flush_object;

/*create a multi-prod multi-cons queue for storing the chunk of fastq file.*/
boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files(64);
boost::atomic<int> num_files {0};

/* dump the contents of a local QF into the main QF */
static bool dump_local_qf_to_main(flush_object *obj)
{
	CQF<KeyObject>::Iterator it = obj->local_cqf->begin();
	do {
		KeyObject hash = it.get_cur_hash();
		int ret = obj->main_cqf->insert(hash, QF_WAIT_FOR_LOCK | QF_KEY_IS_HASH);
		if (ret == QF_NO_SPACE) {
			obj->console->error("The CQF is full. Please rerun the with a larger size.");
			return false;
		}
		++it;
	} while (!it.done());
	obj->local_cqf->reset();

	return true;
}

/* convert a chunk of the fastq file into kmers */
bool reads_to_kmers(chunk &c, flush_object *obj)
{
	auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		std::string read(fs, fe-fs);

start_read:
		if (read.length() < obj->ksize) // start with the next read if length is smaller than K
			goto next_read;
		{
			__int128_t first = 0;
			__int128_t first_rev = 0;
			__int128_t item = 0;
			for(uint32_t i = 0; i < obj->ksize; i++) { //First kmer
				uint8_t curr = Kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					if (i + 1 < read.length())
						read = read.substr(i + 1, read.length());
					else
						continue;
					goto start_read;
				}
				first = first | curr;
				first = first << 2;
			}
			first = first >> 2;
			first_rev = Kmer::reverse_complement(first, obj->ksize);

			if (Kmer::compare_kmers(first, first_rev))
				item = first;
			else
				item = first_rev;

			/*
			 * first try and insert in the main QF.
			 * If lock can't be acquired in the first attempt then
			 * insert the item in the local QF.
			 */
			KeyObject k(item, 0, 1);
			int ret = obj->main_cqf->insert(k, QF_TRY_ONCE_LOCK);
			if (ret == QF_NO_SPACE) {
				obj->console->error("The CQF is full. Please rerun the with a larger size.");
				exit(1);
			} else if (ret == -2) {
				obj->local_cqf->insert(k, QF_NO_LOCK);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					if (!dump_local_qf_to_main(obj))
						return false;
					obj->count = 0;
				}
			}

			uint64_t next = (first << 2) & BITMASK(2 * obj->ksize);
			uint64_t next_rev = first_rev >> 2;

			for(uint32_t i = obj->ksize; i < read.length(); i++) { //next kmers
				uint8_t curr = Kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					if (i + 1 < read.length())
						read = read.substr(i + 1, read.length());
					else
						continue;
					goto start_read;
				}
				next |= curr;
				uint64_t tmp = Kmer::reverse_complement_base(curr);
				tmp <<= (obj->ksize * 2 - 2);
				next_rev = next_rev | tmp;
				if (Kmer::compare_kmers(next, next_rev))
					item = next;
				else
					item = next_rev;

				/*
				 * first try and insert in the main QF.
				 * If lock can't be accuired in the first attempt then
				 * insert the item in the local QF.
				 */
				KeyObject k(item, 0, 1);
				ret = obj->main_cqf->insert(k, QF_TRY_ONCE_LOCK);
				if (ret == QF_NO_SPACE) {
					obj->console->error("The CQF is full. Please rerun the with a larger size.");
					exit(1);
				} else if (ret == -2) {
					obj->local_cqf->insert(k, QF_NO_LOCK);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						if (!dump_local_qf_to_main(obj))
							return false;
						obj->count = 0;
					}
				}

				next = (next << 2) & BITMASK(2*obj->ksize);
				next_rev = next_rev >> 2;
			}
		}

next_read:
		fs = ++fe;		// increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
		fs++; // increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
		fs++; // increment the pointer
	}
	free(c.get_reads());

	return true;
}

/* read a part of the fastq file, parse it, convert the reads to kmers, and
 * insert them in the CQF
 */
static bool fastq_to_uint64kmers_prod(flush_object* obj)
{
	file_pointer* fp;

	while (num_files) {
		while (ip_files.pop(fp)) {
			if (reader::fastq_read_parts(fp->mode, fp)) {
				ip_files.push(fp);
				chunk c(fp->part, fp->size);
				if (!reads_to_kmers(c, obj)) {
					obj->console->error("Insertion in the CQF failed.");
					abort();
				}
			} else {
				/* close the file */
				if (fp->mode == 0)
					fclose(fp->freader->in);
				else if (fp->mode == 1)
					gzclose(fp->freader->in_gzip);
				else if (fp->mode == 2)
					if (fp->freader->in) {
						BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
						fclose(fp->freader->in);
					}
				delete[] fp->part_buffer;
				delete fp;
				num_files--;
			}
		}
	}
	if (obj->count) {
		if (!dump_local_qf_to_main(obj)) {
			obj->console->error("Insertion in the CQF failed.");
			abort();
		}
		obj->count = 0;
	}

	return true;
}

/* main method */
int count_main(CountOpts &opts)
{
	int mode = 0;

	spdlog::logger* console = opts.console.get();

	if (opts.exact && opts.ksize > 32) {
		console->error("Does not support k-mer size > 32 for squeakr-exact.");
		return 1;
	}

	if (opts.qbits == 0)
		opts.qbits = 28;

	if (opts.numthreads == 0)
		opts.numthreads = std::thread::hardware_concurrency();

	enum qf_hashmode hash = QF_HASH_DEFAULT;
	int num_hash_bits = opts.qbits+8;	// we use 8 bits for remainders in the main QF
	if (opts.exact) {
		num_hash_bits = 2*opts.ksize; // Each base 2 bits.
		hash = QF_HASH_INVERTIBLE;
	}

	std::string ser_ext(".squeakr");
	std::string log_ext(".log");
	std::string cluster_ext(".cluster");
	std::string freq_ext(".freq");
	struct timeval start1, start2, end1, end2;
	struct timezone tzp;
	uint32_t OVERHEAD_SIZE = 65535;

	std::string filepath(opts.filenames.front());
	auto const pos = filepath.find_last_of('.');
	std::string input_ext = filepath.substr(pos + 1);
	if (input_ext == std::string("fastq") || input_ext == std::string("fq"))
		mode = 0;
	else if (input_ext == std::string("gz"))
		mode = 1;
	else if (input_ext == std::string("bz2"))
		mode = 2;
	else {
		console->error("Does not support this input file type.");
		return 1;
	}

	for( auto& fn : opts.filenames ) {
		auto* fr = new reader;
		if (reader::getFileReader(mode, fn.c_str(), fr)) {
			file_pointer* fp = new file_pointer;
			fp->mode = mode;
			fp->freader.reset(fr);
			fp->part_buffer = new char[OVERHEAD_SIZE];
			ip_files.push(fp);
			num_files++;
		} else {
			delete fr;
		}
	}

	std::string ds_file = opts.output_file;

	//Initialize the main  QF
	CQF<KeyObject> cqf(opts.qbits, num_hash_bits, hash, SEED);
	if (opts.numthreads == 1)
		cqf.set_auto_resize();
	CQF<KeyObject> *local_cqfs = (CQF<KeyObject>*)calloc(MAX_NUM_THREADS,
																											 sizeof(CQF<KeyObject>));

	boost::thread_group prod_threads;

	for (int i = 0; i < opts.numthreads; i++) {
		local_cqfs[i] = CQF<KeyObject>(QBITS_LOCAL_QF, num_hash_bits, hash, SEED);
		flush_object* obj = (flush_object*)malloc(sizeof(flush_object));
		obj->local_cqf = &local_cqfs[i];
		obj->main_cqf = &cqf;
		obj->ksize = opts.ksize;
		obj->exact = opts.exact;
		obj->count = 0;
		obj->console = console;
		prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod,
																							obj));
	}

	console->info("Reading from the fastq file and inserting in the CQF.");
	gettimeofday(&start1, &tzp);
	prod_threads.join_all();

	// Resize the CQF if:
	//      there is cutoff value greater than 1
	//      the final CQF doesn't need counts
	//      the default size if used
	if (opts.cutoff > 1 || opts.contains_counts == 0 || !opts.setqbits) {
		console->info("Filtering k-mers based on the cutoff.");
		uint64_t num_kmers{0}, estimated_size{0}, log_estimated_size{0};
		CQF<KeyObject>::Iterator it = cqf.begin();
		while (!it.done()) {
			KeyObject hash = it.get_cur_hash();
			if (hash.count >= (uint32_t)opts.cutoff)
				num_kmers++;
			//console->info("hash fraction: {}", k.key / (float)cqf.range());
			if (cqf.get_unique_index(hash, QF_KEY_IS_HASH) / (float)(1ULL <<
																														opts.qbits) >
					0.05) {
				estimated_size = num_kmers * (cqf.range() / hash.key);
				//console->info("estimated size: {}", estimated_size);
				if (opts.contains_counts == 1)
					estimated_size *= 3;    // to account for counts.
				log_estimated_size = ceil(log2(estimated_size));
				uint64_t total_slots = 1ULL << log_estimated_size;
				//console->info("estimated size after ceiling: {}", 1ULL << log_estimated_size);
				if ((total_slots-estimated_size)/(float)estimated_size < 0.1)
					log_estimated_size += 1;
				//console->info("estimated size after ceiling: {}", 1ULL << log_estimated_size);
				break;
			}
			++it;
		}
		console->info("Estimated size of the final CQF: {}", log_estimated_size);
		if (cqf.numslots() > (1ULL << log_estimated_size)) {
			CQF<KeyObject> filtered_cqf(log_estimated_size, num_hash_bits, hash, SEED);
			filtered_cqf.set_auto_resize();
			it = cqf.begin();
			uint64_t max_cnt = 0;
			while (!it.done()) {
				KeyObject hash = it.get_cur_hash();
				if (hash.count >= (uint32_t)opts.cutoff) {
					int ret;
					if (opts.contains_counts == 1)
						ret = filtered_cqf.insert(hash, QF_NO_LOCK | QF_KEY_IS_HASH);
					else {
						hash.count = 1;
						ret = filtered_cqf.insert(hash, QF_NO_LOCK | QF_KEY_IS_HASH);
					}
					if (ret == QF_NO_SPACE) {
						console->error("The CQF is full. Estimated size of the final CQF is wrong.");
						exit(1);
					}
				}
				if (max_cnt < hash.count)
					max_cnt = hash.count;
				++it;
			}
			cqf = filtered_cqf;
		}
	}
	// serialize the CQF to disk.
	cqf.serialize(ds_file);
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1, console);

	console->info("Calculating frequency distribution:");
	gettimeofday(&start2, &tzp);
	uint64_t max_cnt = 0;
	CQF<KeyObject>::Iterator it = cqf.begin();
	while (!it.done()) {
		KeyObject hash = it.get_cur_hash();
		if (max_cnt < hash.count)
			max_cnt = hash.count;
		++it;
	}

	gettimeofday(&end2, &tzp);
	print_time_elapsed("", &start2, &end2, console);

	console->info("Maximum freq: {}", max_cnt);

	console->info("Num distinct elem: {}", cqf.dist_elts());
	console->info("Total num elems: {}", cqf.total_elts());

	// seek to the end of the file and write the k-mer size
	std::ofstream squeakr_file(ds_file, std::ofstream::out |
														 std::ofstream::app | std::ofstream::binary);
	squeakr_file.seekp(0, squeakr_file.end);
	squeakrconfig config;
	config.kmer_size = opts.ksize;
	config.cutoff = opts.cutoff;
	config.contains_counts = opts.contains_counts;
	squeakr_file.write((const char*)&config, sizeof(config));
	squeakr_file.close();

	return 0;
}

