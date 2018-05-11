/*
 * =====================================================================================
 *
 *       Filename:  chunk.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2016 04:49:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _CHUNK_H_
#define _CHUNK_H_

#include <stdio.h>

namespace kmercounting {
	class chunk {
		public:
			chunk();
			chunk(char *reads, uint32_t size);
			inline char *get_reads();
			inline uint32_t get_size();

		private:
			char *_reads;
			uint32_t _size;
	};

	chunk::chunk()
	{
		_reads = NULL;
		_size = 0;
	}

	chunk::chunk(char *reads, uint32_t size)
	{
		_reads = reads;
		_size = size;
	}

	inline char *chunk::get_reads()
	{
		return _reads;
	}

	inline uint32_t chunk::get_size()
	{
		return _size;
	}
}

#endif
