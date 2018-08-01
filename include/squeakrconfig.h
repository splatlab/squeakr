/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#ifndef _SQUEAKR_CONFIG_H_
#define _SQUEAKR_CONFIG_H_

#include <stdio.h>

#define VERSION 2

typedef struct __attribute__ ((__packed__)) squeakrconfig {
	uint64_t kmer_size;
	uint64_t cutoff;
	uint64_t contains_counts;
	uint64_t endianness{0x0102030405060708ULL};
	uint32_t version{VERSION};
} squeakrconfig;

#endif
