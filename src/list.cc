/*
 * =============================================================================
 *
 *       Filename:  list.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/10/2018 01:58:17 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =============================================================================
 */

#include <iostream>
#include <cstring>
#include <vector>
#include <set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include	<stdlib.h>

#include "clipp.h"
#include "ProgOpts.h"
#include "gqf_cpp.h"
#include "kmer.h"

/* 
 * ===  FUNCTION  ==============================================================
 *         Name:  main
 *  Description:  
 * =============================================================================
 */
int list_main(ListOpts& opts)
{
	struct timeval start, end;
	struct timezone tzp;
	spdlog::logger* console = opts.console.get();

	//Initialize the QF
	console->info("Reading kmers into the QF off the disk.");
	CQF<KeyObject> cqf(opts.squeakr_file, LOCKS_FORBIDDEN, FREAD);

	if (!cqf.is_exact()) {
		console->error("The file is not generated using Squeakr-exact.");
		return 1;
	}
	
	uint64_t kmer_size = cqf.keybits() / 2;
	std::ofstream opfile(opts.output_file.c_str(), std::ofstream::out);

	CQF<KeyObject>::Iterator it = cqf.begin();
	gettimeofday(&start, &tzp);
	do {
		KeyObject k = *it;
		opfile << Kmer::int_to_str(k.key, kmer_size) << "\t" << k.count << std::endl;
		++it;
	} while (!it.done());
	gettimeofday(&end, &tzp);
	opfile.close();
	print_time_elapsed("", &start, &end, console);


	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
