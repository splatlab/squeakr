/*
 * =====================================================================================
 *
 *       Filename:  util.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/10/2018 01:21:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include "util.h"

/* Count distinct items in a sorted list */
uint64_t count_distinct_kmers(std::multiset<__int128_t> kmers)
{
	uint64_t cnt = 0;
	__int128_t curr_kmer = 0;

	for(__int128_t kmer: kmers) {
		if (kmer != curr_kmer) {
			curr_kmer = kmer;
			cnt++;
		}
	}
	return cnt;
}

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end, spdlog::logger* console)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	console->info("{}  Total Time Elapsed: {} seconds", desc,
								std::to_string(time_elapsed));
}

std::string last_part(std::string str, char c) {
	uint64_t found = str.find_last_of(c);
	return str.substr(found + 1);
}

