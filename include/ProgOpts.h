#ifndef __PROG_OPTS__
#define __PROG_OPTS__

#include <memory>
#include "spdlog/spdlog.h"

class CountOpts {
	public:
		int exact {0};
		int ksize;
		int qbits {0};
		int numthreads{0};
		std::string prefix;
		std::string output_dir;
		std::vector<std::string> filenames;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class QueryOpts {
	public:
		std::string cqf_file;
		int ksize {0};
		std::string queryfile;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class InnerProdOpts {
	public:
		std::string cqf_filea;
		std::string cqf_fileb;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

class ListOpts {
	public:
		std::string cqf_file;
		std::string output_file;
		std::shared_ptr<spdlog::logger> console{nullptr};
};

#endif //__MANTIS_PROG_OPTS__
