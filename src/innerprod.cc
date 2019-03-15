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
#include "ProgOpts.h"
#include "gqf_cpp.h"

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int inner_prod_main(InnerProdOpts& opts)
{
	uint64_t inner_prod;
	struct timeval start, end;
	struct timezone tzp;

	srand(time(NULL));
	spdlog::logger* console = opts.console.get();

	//Initialize the QF
	console->info("Mmaping Squeakr files from disk.");
	CQF<KeyObject> cfa(opts.squeakr_filea, FREAD_MODE);
	CQF<KeyObject> cfb(opts.squeakr_fileb, FREAD_MODE);

	if (cfa.seed() != cfb.seed()) {
		console->error("Input CQFs do not have the same seed.");
		return 1;
	}

	console->info("Performing inner product querries.");

	gettimeofday(&start, &tzp);
	inner_prod = cfa.inner_prod(cfb);
	gettimeofday(&end, &tzp);
	console->info("Inner product: {}", inner_prod);
	print_time_elapsed("", &start, &end, console);
	
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */

