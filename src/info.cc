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
#include	<stdlib.h>

#include "clipp.h"
#include "ProgOpts.h"
#include "squeakrconfig.h"
#include "gqf_cpp.h"
#include "kmer.h"

/*
 * ===  FUNCTION  ==============================================================
 *         Name:  main
 *  Description:  
 * =============================================================================
 */
int info_main(InfoOpts& opts)
{
	spdlog::logger* console = opts.console.get();

	//Initialize the QF
	console->info("Reading squeakr file off the disk.");
	CQF<KeyObject> cqf(opts.squeakr_file, FREAD);

	if (!cqf.is_exact()) {
		console->error("The file is not generated using Squeakr-exact.");
		return 1;
	}

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

	std::string contains_counts = config.contains_counts ? "YES" : "NO";
	std::string exact = cqf.is_exact() ? "YES" : "NO";

	console->info("version: {}", config.version);
	console->info("kmer size: {}", config.kmer_size);
	console->info("Cutoff: {}", config.cutoff);
	console->info("Contains Counts: {}", contains_counts);
	console->info("CQF info.....");
	console->info("Is Exact: {}", exact);
	console->info("Total slots: {}", cqf.numslots());
	console->info("Seed: {}", cqf.seed());
	console->info("Keybits: {}", cqf.keybits());
	console->info("Total elements: {}", cqf.total_elts());
	console->info("Distinct elements: {}", cqf.dist_elts());

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
