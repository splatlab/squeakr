/*
 * =====================================================================================
 *
 *       Filename:  kmer_query.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/27/2016 08:48:26 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
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

#include <zlib.h>
#include <bzlib.h>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) - 1ULL)

// A=1, C=0, T=2, G=3
void getRandomKmers(int n, uint64_t range, std::vector<uint64_t>& kmers,
										uint32_t K)
{
	uint64_t kmer;
	for (int j = 0; j < n; j++) {
		kmer = 0;
		for (uint i = 0; i < K; i++) {
			uint8_t c = rand()%4;
			kmer = kmer | c;
			kmer = kmer << 2;
		}
		kmer = kmer >> 2;
		kmers.push_back(kmer%range);
	}
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
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
	CQF<KeyObject> cqf(opts.cqf_file, LOCKS_FORBIDDEN, FREAD);

	if (cqf.is_exact() && opts.ksize != cqf.keybits()) {
		console->error("K-mer size is not correct.");
		return 1;
	}

	console->info("Parsing query file for {}-mers.", opts.ksize);
	Kmer::parse_kmers(opts.queryfile.c_str(), opts.ksize, kmers);
	console->info("Found {} k-mers", kmers.size());

	console->info("Querying kmers in the QF.");
	uint64_t num_not_found = 0;
	gettimeofday(&start, &tzp);
	for (auto it = kmers.begin(); it != kmers.end(); ++it) {
		if (!cqf.query(KeyObject(*it, 0, 0))) {
			num_not_found++;
		}
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end, console);

	console->info("Not find: {}", num_not_found);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

