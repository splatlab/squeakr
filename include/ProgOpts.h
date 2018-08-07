/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#ifndef __PROG_OPTS__
#define __PROG_OPTS__

#include <memory>
#include "spdlog/spdlog.h"

class CountOpts {
	public:
		int exact {0};
		int ksize;
		int cutoff {1};
		int contains_counts{1};
		int qbits {0};
		bool setqbits{false};
		int numthreads{0};
		std::string output_file;
		std::vector<std::string> filenames;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class QueryOpts {
	public:
		std::string squeakr_file;
		std::string queryfile;
		std::string output_file;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class InnerProdOpts {
	public:
		std::string squeakr_filea;
		std::string squeakr_fileb;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class ListOpts {
	public:
		std::string squeakr_file;
		std::string output_file;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

#endif //__MANTIS_PROG_OPTS__
