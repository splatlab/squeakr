/*
 * =====================================================================================
 *
 *       Filename:  util.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/10/2018 01:20:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _UTIL_H_
#define _UTIL_H_

#include <iostream>
#include <algorithm>
#include <memory>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include "reader.h"

#ifdef DEBUG
#define PRINT_DEBUG 1
#else
#define PRINT_DEBUG 0
#endif

#define DEBUG(x) do { \
	if (PRINT_DEBUG) { std::cerr << x << std::endl; } \
} while (0)

#define ERROR(x) do { \
	{ std::cerr << x << std::endl; } \
} while (0)

#define PRINT(x) do { \
	{ std::cout << x << std::endl; } \
} while (0)

struct file_pointer {
	std::unique_ptr<reader> freader{nullptr};
	char* part{nullptr};
	char* part_buffer{nullptr};
	int mode{0};
	uint64_t size{0};
	uint64_t part_filled{0};
};

void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end);
inline bool is_eof(reader &file_reader, int mode);
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos);
static bool fastq_read_parts(int mode, file_pointer *fp);
bool getFileReader(int mode, const char* fastq_file, reader*
									 file_reader);

std::string last_part(std::string str, char c);

#endif
