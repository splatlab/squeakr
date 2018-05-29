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
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
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

#include "clipp.h"
#include "gqf_cpp.h"

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int inner_prod_main(int argc, char *argv[])
{
	std::string ds_filea;
	std::string ds_fileb;
	uint64_t inner_prod;
	struct timeval start, end;
	struct timezone tzp;

	using namespace clipp;
	auto cli = (
							required("-a", "--cqf-file-first") & value("cqf-file-first", ds_filea) % "first input CQF file",
							required("-b", "--cqf-file-second") & value("cqf-file-second", ds_fileb) % "second input CQF file",
							option("-h", "--help")  % "show help"
						 );

	auto res = parse(argc, argv, cli);

	if (!res) {
		std::cout << make_man_page(cli, argv[0]) << "\n";
		return 1;
	}


	srand(time(NULL));

	//Initialize the QF
	std::cout << "Mmap the QF from disk" << std::endl;
	CQF<KeyObject> cfa(ds_filea, LOCKS_FORBIDDEN, FREAD);
	CQF<KeyObject> cfb(ds_fileb, LOCKS_FORBIDDEN, FREAD);

	if (cfa.seed() != cfb.seed()) {
		std::cerr << "Input CQFs do not have the same seed." << std::endl;
		return 1;
	}

	std::cout << "Performing inner product querries." << std::endl;

	gettimeofday(&start, &tzp);
	inner_prod = cfa.inner_prod(cfb);
	gettimeofday(&end, &tzp);
	std::cout << "Inner product: " << inner_prod << std::endl;
	print_time_elapsed("", &start, &end);
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

