/*
 * =====================================================================================
 *
 *       Filename:  squeakrconfig.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/25/2018 05:59:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _SQUEAKR_CONFIG_H_
#define _SQUEAKR_CONFIG_H_

#include <stdio.h>

#define VERSION 2

typedef struct squeakrconfig {
	uint64_t kmer_size;
	uint32_t version{VERSION};
} squeakrconfig;

#endif
