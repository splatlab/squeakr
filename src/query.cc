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
	CQF<KeyObject> cqf(opts.squeakr_file, LOCKS_FORBIDDEN, FREAD);

	// seek to the end of the file and read the k-mer size
	std::ifstream squeakr_file(opts.squeakr_file, std::ofstream::in);
	squeakr_file.seekg(0, squeakr_file.end);
	uint64_t file_size = squeakr_file.tellg();
	uint64_t kmer_size;
	squeakr_file.seekg(file_size - sizeof(kmer_size));
	squeakr_file.read((char*)&kmer_size, sizeof(kmer_size));
	squeakr_file.close();
	console->info("kmer size: {}", kmer_size);

	if (cqf.is_exact() && kmer_size*2 != cqf.keybits()) {
		console->error("K-mer size is not correct.");
		return 1;
	}

	console->info("Parsing query file for {}-mers.", kmer_size);
	Kmer::parse_kmers(opts.queryfile.c_str(), kmer_size, kmers);
	console->info("Found {} k-mers", kmers.size());

	std::ofstream opfile(opts.output_file.c_str(), std::ofstream::out);

	console->info("Querying kmers in the QF.");
	uint64_t num_not_found = 0;
	gettimeofday(&start, &tzp);
	for (auto it = kmers.begin(); it != kmers.end(); ++it) {
		uint64_t count = cqf.query(KeyObject(*it, 0, 0));
		if (count == 0)
			num_not_found++;
		else
			opfile << Kmer::int_to_str(*it, kmer_size) << "\t" << count << std::endl;
	}
	gettimeofday(&end, &tzp);
	opfile.close();
	print_time_elapsed("", &start, &end, console);

	console->info("Not found: {}", num_not_found);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

