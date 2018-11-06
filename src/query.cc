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
#include "kmer.h"
#include "squeakrconfig.h"

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
int query_main(QueryOpts& opts)
{
	struct timeval start, end;
	struct timezone tzp;
	std::unordered_set<uint64_t> kmers;

	srand(time(NULL));
	spdlog::logger* console = opts.console.get();

	//Initialize the QF
	console->info("Reading kmers into the QF off the disk.");
	CQF<KeyObject> cqf(opts.squeakr_file, FREAD);

	// seek to the end of the file and read the k-mer size
	std::ifstream squeakr_file(opts.squeakr_file, std::ofstream::in);
	squeakr_file.seekg(0, squeakr_file.end);
	uint64_t file_size = squeakr_file.tellg();
	squeakrconfig config;
	squeakr_file.seekg(file_size - sizeof(squeakrconfig));
	squeakr_file.read((char*)&config, sizeof(config));
	squeakr_file.close();
	if (config.version != VERSION) {
		console->error("Squeakr index version is invalid. Expected: {} Available: {}",
									 VERSION, config.version);
		exit(1);
	}
	console->info("kmer size: {}, version: {}", config.kmer_size, config.version);

	if (cqf.is_exact() && config.kmer_size*2 != cqf.keybits()) {
		console->error("K-mer size is not correct.");
		return 1;
	}

	console->info("Parsing query file for {}-mers.", config.kmer_size);
	Kmer::parse_kmers(opts.queryfile.c_str(), config.kmer_size, kmers);
	console->info("Found {} k-mers", kmers.size());

	std::ofstream opfile(opts.output_file.c_str(), std::ofstream::out);

	console->info("Querying kmers in the QF.");
	uint64_t num_not_found = 0;
	gettimeofday(&start, &tzp);
	for (auto it = kmers.begin(); it != kmers.end(); ++it) {
		uint64_t count = cqf.query(KeyObject(*it, 0, 0), 0);
		if (count == 0) {
			console->error(Kmer::int_to_str(*it, config.kmer_size));
			num_not_found++;
		}
		else
			opfile << Kmer::int_to_str(*it, config.kmer_size) << "\t" << count << std::endl;
	}
	gettimeofday(&end, &tzp);
	opfile.close();
	print_time_elapsed("", &start, &end, console);

	console->info("Not found: {}", num_not_found);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

