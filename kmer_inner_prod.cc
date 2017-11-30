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
#include "threadsafe-gqf/gqf.h"

using namespace std;

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(string desc, struct timeval* start, struct timeval* end) 
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	cout << desc << "Total Time Elapsed: " << to_string(time_elapsed) << " seconds" << endl;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{
	QF cfa, cfb;
	
	//if (argc == 2) {
		//string arg_help(argv[1]);
		//if (arg_help.compare("-h") != 0 || arg_help.compare("-help") != 0) {
			//cout << "./squeakr-inner-product [OPTIONS]" << endl
					 //<< "file1: dataset 1 Squeakr representation" << endl
					 //<< "file2: dataset 2 Squeakr representation" << endl;
			//exit(0);
		//}
	//}

	string ds_filea;
	string ds_fileb;
	uint64_t inner_prod;
	struct timeval start, end;
	struct timezone tzp;

	using namespace clipp;
	auto cli = (
							required("-a", "--cqf-file-first") & value("cqf-file-first", ds_filea) % "first input CQF file",
							required("-b", "--cqf-file-second") & value("cqf-file-second", ds_filea) % "second input CQF file",
							option("-h", "--help")  % "show help"
						 );

	auto res = parse(argc, argv, cli);

	if (!res) {
		std::cout << make_man_page(cli, argv[0]) << "\n";
		return 1;
	}


	srand(time(NULL));

	//Initialize the QF
	cout << "Mmap the QF from disk" << endl;
	qf_read(&cfa, ds_filea.c_str());
	qf_read(&cfb, ds_fileb.c_str());

	if (cfa.metadata->seed != cfb.metadata->seed) {
		cerr << "Input CQFs do not have the same seed." << std::endl;
		return 1;
	}

	cout << "Performing inner product querries." << endl;

	gettimeofday(&start, &tzp);
	inner_prod = qf_inner_product(&cfa, &cfb);
	gettimeofday(&end, &tzp);
	cout << "Inner product: " << inner_prod << "\n" << endl;
	print_time_elapsed("", &start, &end);
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

