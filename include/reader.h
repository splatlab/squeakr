/*
 * =====================================================================================
 *
 *       Filename:  reader.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/06/2016 09:56:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <stdio.h>

#include <zlib.h>
#include <bzlib.h>

#ifndef _READER_H_
#define _READER_H_

namespace kmercounting {
	class reader {
		public:
			reader();
			reader(FILE *in, gzFile in_gzip, BZFILE *in_bzip2, int bzerror);

			FILE *in = nullptr;
			gzFile in_gzip = nullptr;
			BZFILE *in_bzip2 = nullptr;
			int bzerror;
	};

	reader::reader()
	{
		in = nullptr;
		in_gzip = nullptr;
		in_bzip2 = nullptr;
		bzerror = 0;
	}

	reader::reader(FILE *_in, gzFile _in_gzip, BZFILE *_in_bzip2, int _bzerror)
	{
		in = _in;
		in_gzip = _in_gzip;
		in_bzip2 = _in_bzip2;
		bzerror = _bzerror;
	}

}

#endif 
