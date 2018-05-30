#ifndef __PROG_OPTS__
#define __PROG_OPTS__

#include <memory>
#include "spdlog/spdlog.h"

class CountOpts {
 public:
	int exact;
	int ksize;
	int qbits;
	int numthreads{1};
	std::string filename;
  std::string prefix;
	std::vector<std::string> filenames;
  std::shared_ptr<spdlog::logger> console{nullptr};
};

class QueryOpts {
 public:
  std::string cqf_file;
	int ksize;
	int num_query;
	int random;
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
	std::string filename;
};

#endif //__MANTIS_PROG_OPTS__
