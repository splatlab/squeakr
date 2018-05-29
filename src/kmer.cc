#include <fstream>
#include "kmer.h"


/*return the integer representation of the base */
inline char Kmer::map_int(uint8_t base)
{
	switch(base) {
		case DNA_MAP::A: { return 'A'; }
		case DNA_MAP::T: { return 'T'; }
		case DNA_MAP::C: { return 'C'; }
		case DNA_MAP::G: { return 'G'; }
		default:  { return DNA_MAP::G+1; }
	}
}

/*return the integer representation of the base */
inline uint8_t Kmer::map_base(char base)
{
	switch(base) {
		case 'A': { return DNA_MAP::A; }
		case 'T': { return DNA_MAP::T; }
		case 'C': { return DNA_MAP::C; }
		case 'G': { return DNA_MAP::G; }
		default:  { return DNA_MAP::G+1; }
	}
}

/**
 * Converts a string of "ATCG" to a uint64_t
 * where each character is represented by using only two bits
 */
__int128_t str_to_int(std::string str)
{
	__int128_t strint = 0;
	for (auto it = str.begin(); it != str.end(); it++) {
		uint8_t curr = 0;
		switch (*it) {
			case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
		}
		strint = strint | curr;
		strint = strint << 2;
	}
	return strint >> 2;
}

/**
 * Converts a uint64_t to a string of "ACTG"
 * where each character is represented by using only two bits
 */
std::string int_to_str(__int128_t kmer, uint64_t kmer_size)
{
	uint8_t base;
	std::string str;
	for (uint32_t i = kmer_size; i > 0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = Kmer::map_int(base);
		str.push_back(chr);
	}
	return str;
}

/* Return the reverse complement of a base */
inline int Kmer::reverse_complement_base(int x) { return 3 - x; }

/* Calculate the revsese complement of a kmer */
__int128_t Kmer::reverse_complement(__int128_t kmer, uint64_t kmer_size)
{
	__int128_t rc = 0;
	uint8_t base = 0;
	for (uint32_t i = 0; i < kmer_size; i++) {
		base = kmer & 3ULL;
		base = reverse_complement_base(base);
		kmer >>= 2;
		rc |= base;
		rc <<= 2;
	}
	rc >>=2;
	return rc;
}

/* Compare the kmer and its reverse complement and return the result 
 * Return true if the kmer is greater than or equal to its
 * reverse complement. 
 * */
inline bool Kmer::compare_kmers(__int128_t kmer, __int128_t kmer_rev)
{
	return kmer >= kmer_rev;
}

