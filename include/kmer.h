/*
 * =====================================================================================
 *
 *       Filename:  kmer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2016 05:06:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _KMER_H_
#define _KMER_H_

#include <stdio.h>
#include <string>

#include "hashutil.h"
#include "common_types.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

namespace kmercounting {
	class Kmer {
		public:
			static inline char map_int(uint8_t base);
			static inline uint8_t map_base(char base);
			static __in128_t str_to_int(string str);
			static string int_to_str(__int128_t kmer, uint64_t kmer_size);
			static inline int reverse_complement_base(int x);
			static __int128_t reverse_complement(__int128_t kmer, uint64_t kmer_size);
			static inline bool compare_kmers(__int128_t kmer, __int128_t kmer_rev);

		private:
			Kmer();
	};
}
#endif
