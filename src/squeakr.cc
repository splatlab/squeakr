/*
 * =====================================================================================
 *
 *       Filename:  squeakr.cc
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
 * =====================================================================================
 */

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <exception>

#include "ProgOpts.h"
#include "clipp.h"
#include "SqueakrFS.h"

template <typename T>
void explore_options_verbose(T& res) {
  if(res.any_error()) { std::cerr << "error\n"; }

  //aggregated errors
  if(res.unmapped_args_count()) { std::cerr << "error unmapped args count\n"; /* ... */ }
  if(res.any_bad_repeat()) { std::cerr << "error bad repeat \n"; /* ... */ }
  if(res.any_blocked())    { std::cerr << "error blocked \n"; /* ... */ }
  if(res.any_conflict())   { std::cerr << "error conflict\n"; /* ... */ }

  for(const auto& m : res.missing()) { 
    std::cerr << "missing " << m.param() << " after index " << m.after_index() << '\n';
  }

  //per-argument mapping
  for(const auto& m : res) {
    std::cerr << m.index() << ": " << m.arg() << " -> " << m.param();
    std::cerr << " repeat #" << m.repeat();
    if(m.blocked()) std::cerr << " blocked";
    if(m.conflict()) std::cerr << " conflict";
    std::cerr<< '\n';
  }
}

int query_main (QueryOpts& opt);
int count_main (CountOpts& opt);
int inner_prod_main (InnerProdOpts& opt);
int list_main (ListOpts& opt);

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:
 * ============================================================================
 */
int main ( int argc, char *argv[] ) {
  using namespace clipp;
  enum class mode {count, query, inner_prod, list, help};
  mode selected = mode::help;

  auto console = spdlog::stdout_color_mt("squeakr_console");

  CountOpts countopt;
  QueryOpts queryopt;
  InnerProdOpts innerprodopt;
	ListOpts listopt;
  countopt.console = console;
  queryopt.console = console;
  innerprodopt.console = console;
	listopt.console = console;

  auto ensure_file_exists = [](const std::string& s) -> bool {
    bool exists = squeakr::fs::FileExists(s.c_str());
    if (!exists) {
      std::string e = "The required input file " + s + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

  auto ensure_dir_exists = [](const std::string& s) -> bool {
    bool exists = squeakr::fs::DirExists(s.c_str());
    if (!exists) {
      std::string e = "The required input directory " + s + " does not seem to exist.";
      throw std::runtime_error{e};
    }
    return true;
  };

	auto count_mode = (
									command("count").set(selected, mode::count),
									option("-e", "--exact").set(countopt.exact, 1) %
									"squeakr-exact (default is Squeakr approximate)",
									required("-k","--kmer") & value("k-size", countopt.ksize) %
									"length of k-mers to count",
									required("-s","--log-slots") & value("log-slots",
																											 countopt.qbits) % "log of number of slots in the CQF",
									option("-t","--threads") & value("num-threads",
																										 countopt.numthreads) %
									"number of threads to use to count (default = number of hardware threads)",
									required("-p","--prefix") & value("prefix",
																														 countopt.prefix)
									% "output file prefix",
									option("-o","--output-dir") & value("out-dir",
																											countopt.output_dir) %
									"directory where output should be written (default = \"./\")",
									values("files", countopt.filenames) % "list of files to be counted(supported files: fastq and compressed gzip or bzip2 fastq files)"
									//option("-h", "--help")      % "show help"
						 );

	auto query_mode = (
							command("query").set(selected, mode::query),
							required("-f", "--cqf-file") & value("cqf-file",
																									 queryopt.cqf_file) % "input CQF file",
							required("-k","--kmer") & value("k-size", queryopt.ksize) %
							"length of k-mers to query. Must be same the as the size of counted k-mers",
							required("-q","--query-file") & value("query-file",
																								queryopt.queryfile) % "input query file",
							required("-o", "--output-file") & value("output-file",
																									 queryopt.output_file) % "output file"
							//option("-h", "--help")  % "show help"
						 );

	auto inner_prod_mode = (
							command("inner_prod").set(selected, mode::inner_prod),
							required("-a", "--cqf-file-first") & value("cqf-file-first",
																												 innerprodopt.cqf_filea)
							% "first input CQF file",
							required("-b", "--cqf-file-second") & value("cqf-file-second",
																													innerprodopt.cqf_fileb)
							% "second input CQF file"
							//option("-h", "--help")  % "show help"
						 );

	auto list_mode = (
							command("list").set(selected, mode::list),
							required("-f", "--cqf-file") & value("cqf-file",
																									 listopt.cqf_file) % "input CQF file",
							required("-o", "--output-file") & value("output-file",
																									 listopt.output_file) % "output file"
							//option("-h", "--help")  % "show help"
							);

  auto cli = (
							(count_mode | query_mode | inner_prod_mode | list_mode |
							 command("help").set(selected,mode::help) ),
							option("-v", "--version").call([]{std::cout << "version 1.0\n\n";}).doc("show version")
							);

  assert(count_mode.flags_are_prefix_free());
  assert(query_mode.flags_are_prefix_free());
  assert(inner_prod_mode.flags_are_prefix_free());
  assert(list_mode.flags_are_prefix_free());

  decltype(parse(argc, argv, cli)) res;
  try {
    res = parse(argc, argv, cli);
  } catch (std::exception& e) {
		std::cout << "\n\nParsing command line failed with exception: " <<
			e.what() << "\n";
    std::cout << "\n\n";
    std::cout << make_man_page(cli, "squeakr");
    return 1;
  }

  if(res) {
    switch(selected) {
    case mode::count: count_main(countopt);  break;
    case mode::query: query_main(queryopt);  break;
    case mode::inner_prod: inner_prod_main(innerprodopt);  break;
    case mode::list: list_main(listopt);  break;
    case mode::help: std::cout << make_man_page(cli, "squeakr"); break;
    }
  } else {
    auto b = res.begin();
    auto e = res.end();
    if (std::distance(b,e) > 0) {
      if (b->arg() == "count") {
        std::cout << make_man_page(count_mode, "squeakr");
      } else if (b->arg() == "query") {
        std::cout << make_man_page(query_mode, "squeakr");
      } else if (b->arg() == "inner_prod") {
        std::cout << make_man_page(inner_prod_mode, "squeakr");
      } else if (b->arg() == "list") {
        std::cout << make_man_page(list_mode, "squeakr");
      } else {
        std::cout << "There is no command \"" << b->arg() << "\"\n";
        std::cout << usage_lines(cli, "squeakr") << '\n';
      }
    } else {
      std::cout << usage_lines(cli, "squeakr") << '\n';
    }
  }

  return 0;
}
