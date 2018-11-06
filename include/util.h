/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
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

#include "spdlog/spdlog.h"

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

void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end, spdlog::logger* console);

std::string last_part(std::string str, char c);

#endif
