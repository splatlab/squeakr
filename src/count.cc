/*
 * ============================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/06/2016 09:56:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
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
} flush_object;

/*create a multi-prod multi-cons queue for storing the chunk of fastq file.*/
boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files(64);
boost::atomic<int> num_files {0};

/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(flush_object *obj)
{
	CQF<KeyObject>::Iterator it = obj->local_cqf->begin();
	do {
		KeyObject k = *it;
		obj->main_cqf->insert(k);
		++it;
	} while (!it.done());
	obj->local_cqf->reset();
}

/* convert a chunk of the fastq file into kmers */
void reads_to_kmers(chunk &c, flush_object *obj)
{
	auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		std::string read(fs, fe-fs);
		/*cout << read << std::endl;*/

start_read:
		if (read.length() < obj->ksize) // start with the next read if length is smaller than K
			goto next_read;
		{
			__int128_t first = 0;
			__int128_t first_rev = 0;
			__int128_t item = 0;
			for(int i = 0; i < obj->ksize; i++) { //First kmer
				uint8_t curr = Kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i + 1, read.length());
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
			 * If lock can't be accuired in the first attempt then
			 * insert the item in the local QF.
			 */
			KeyObject k(item%obj->main_cqf->range(), 0, 1);
			if (!obj->main_cqf->insert(k)) {
				obj->local_cqf->insert(k);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}

			uint64_t next = (first << 2) & BITMASK(2 * obj->ksize);
			uint64_t next_rev = first_rev >> 2;

			for(uint32_t i = obj->ksize; i < read.length(); i++) { //next kmers
				uint8_t curr = Kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i + 1, read.length());
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
				KeyObject k(item%obj->main_cqf->range(), 0, 1);
				if (!obj->main_cqf->insert(k)) {
					obj->local_cqf->insert(k);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						dump_local_qf_to_main(obj);
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
				reads_to_kmers(c, obj);
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
		dump_local_qf_to_main(obj);
		obj->count = 0;
	}

	return true;
}

/* main method */
int count_main(CountOpts &opts)
{
	int mode = 0;

	if (opts.exact && opts.ksize > 32) {
		std::cout << "Does not support k-mer size > 32 for squeakr-exact." <<
			std::endl;
		return 1;
	}
	
	enum hashmode hash = DEFAULT;
	int num_hash_bits = opts.qbits+8;	// we use 8 bits for remainders in the main QF
	if (opts.exact) {
		num_hash_bits = 2*opts.ksize; // Each base 2 bits.
		hash = INVERTIBLE;
	}

	std::string ser_ext(".cqf");
	std::string log_ext(".log");
	std::string cluster_ext(".cluster");
	std::string freq_ext(".freq");
	struct timeval start1, start2, end1, end2;
	struct timezone tzp;
	uint32_t OVERHEAD_SIZE = 65535;

	std::string filepath(opts.filenames.front());
	auto const pos = filepath.find_last_of('.');
	std::string input_ext = filepath.substr(pos + 1);
	if (input_ext.compare(std::string("fastq")) == 0 ||
			input_ext.compare(std::string("fq")))
		mode = 0;
	else if (input_ext.compare(std::string("gz")))
		mode = 1;
	else if (input_ext.compare(std::string("bz2")))
		mode = 2;
	else {
		std::cout << "Does not support this input file type." << std::endl;
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

	if (opts.prefix.back() != '/') {
		opts.prefix += '/';
	}
	std::string ds_file =      opts.prefix + opts.filename + ser_ext;
	std::string log_file =     opts.prefix + opts.filename + log_ext;
	std::string cluster_file = opts.prefix + opts.filename + cluster_ext;
	std::string freq_file =    opts.prefix + opts.filename + freq_ext;

	// A random large prime number.
	uint32_t seed = 2038074761;
	//Initialize the main  QF
	CQF<KeyObject> cqf(opts.qbits, num_hash_bits, LOCKS_OPTIONAL, hash, seed);
	CQF<KeyObject> *local_cqfs = (CQF<KeyObject>*)calloc(MAX_NUM_THREADS,
																											 sizeof(CQF<KeyObject>));

	boost::thread_group prod_threads;

	for (int i = 0; i < opts.numthreads; i++) {
		local_cqfs[i] = CQF<KeyObject>(QBITS_LOCAL_QF, num_hash_bits,
																	LOCKS_FORBIDDEN, hash, seed);
		flush_object* obj = (flush_object*)malloc(sizeof(flush_object));
		obj->local_cqf = &local_cqfs[i];
		obj->main_cqf = &cqf;
		obj->ksize = opts.ksize;
		obj->exact = opts.exact;
		obj->count = 0;
		prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod,
																							obj));
	}

	std::cout << "Reading from the fastq file and inserting in the QF" <<
		std::endl;
	gettimeofday(&start1, &tzp);
	prod_threads.join_all();
	cqf.serialize(ds_file);
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);

	std::cout << "Calc freq distribution: " << std::endl;
	//ofstream freq_file;
	//freq_file.open(freq_file.c_str());
	uint64_t max_cnt = 0;
	CQF<KeyObject>::Iterator it = cqf.begin();
	gettimeofday(&start2, &tzp);
	do {
		KeyObject k = *it;
		//freq_file << key << " " << count << std::endl;
		if (max_cnt < k.count)
			max_cnt = k.count;
		++it;
	} while (!it.done());
	gettimeofday(&end2, &tzp);
	print_time_elapsed("", &start2, &end2);

	std::cout << "Maximum freq: " << max_cnt << std::endl;
	//freq_file.close();

	std::cout << "Num distinct elem: " << cqf.dist_elts() << std::endl;
	std::cout << "Total num elems: " << cqf.total_elts() << std::endl;

	return 0;
}

