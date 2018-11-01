/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include <iostream>
#include <fstream>
#include <stdio.h>

#include <zlib.h>
#include <bzlib.h>

#ifndef _READER_H_
#define _READER_H_

struct file_pointer;

class reader {
	public:
		reader();
		reader(FILE *in, gzFile in_gzip, BZFILE *in_bzip2, int bzerror);

		bool is_eof(int mode);

		static bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos);
		static bool fastq_read_parts(int mode, file_pointer *fp);
		static bool getFileReader(int mode, const char* fastq_file, reader*
															file_reader);

		FILE *in = nullptr;
		gzFile in_gzip = nullptr;
		BZFILE *in_bzip2 = nullptr;
		int bzerror;
};

struct file_pointer {
	std::unique_ptr<reader> freader{nullptr};
	char* part{nullptr};
	char* part_buffer{nullptr};
	int mode{0};
	uint64_t size{0};
	uint64_t part_filled{0};
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

/* check if it's the end of the file. */
bool reader::is_eof(int mode) {
	if (mode == 0)
		return feof(in) != 0;
	else if (mode == 1)
		return gzeof(in_gzip) != 0;
	else if (mode == 2)
		return bzerror == BZ_STREAM_END;

	return true;
}

/* move the pointer to the end of the next newline. */
bool reader::skip_next_eol(char *part, int64_t &pos, int64_t max_pos) {
	int64_t i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' ||
																								 part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;
	pos = i+1;

	return true;
}

/* read a part of the fastq file. */
bool reader::fastq_read_parts(int mode, file_pointer *fp) {
	char *& _part = (fp->part);
	uint64_t& _size = fp->size;
	char*& part_buffer = (fp->part_buffer);
	uint64_t& part_filled = fp->part_filled;
	reader& file_reader = *(fp->freader.get());

	uint32_t OVERHEAD_SIZE = 65535;
	uint64_t part_size = 1ULL << 25;
	char *part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
	memcpy(part, part_buffer, part_filled);

	if(file_reader.is_eof(mode))
		return false;

	uint64_t readed = 0;

	if (mode == 0)
		readed = fread(part+part_filled, 1, part_size, file_reader.in);
	else if (mode == 1)
		readed = gzread(file_reader.in_gzip, part+part_filled, (int) part_size);
	else if (mode == 2)
		readed = BZ2_bzRead(&file_reader.bzerror, file_reader.in_bzip2,
												part+part_filled, (int) part_size);
	else 
		readed = 0;

	int64_t total_filled = part_filled + readed;
	int64_t i;
	if(part_filled >= OVERHEAD_SIZE)
	{
		std::cout << "Error: Wrong input file!" << std::endl;
		exit(EXIT_FAILURE);
	}
	if(file_reader.is_eof(mode))
	{
		_part = part;
		_size = total_filled;
		part = NULL;
		return true;
	}
	// Looking for a FASTQ record at the end of the area
	{
		int64_t line_start[9];
		int32_t j;
		i = total_filled - OVERHEAD_SIZE / 2;
		for(j = 0; j < 9; ++j)
		{
			if(!skip_next_eol(part, i, total_filled))
				break;
			line_start[j] = i;
		}
		_part = part;
		if(j < 9)
			_size = 0;
		else
		{
			int k;
			for(k = 0; k < 4; ++k)
			{
				if(part[line_start[k]+0] == '@' && part[line_start[k+2]+0] == '+')
				{
					if(part[line_start[k+2]+1] == '\n' || part[line_start[k+2]+1] == '\r')
						break;
					if(line_start[k+1]-line_start[k] == line_start[k+3]-line_start[k+2] &&
						 memcmp(part+line_start[k]+1, part+line_start[k+2]+1,
										line_start[k+3]-line_start[k+2]-1) == 0)
						break;
				}
			}
			if(k == 4)
				_size = 0;
			else
				_size = line_start[k];
		}
	}

	std::copy(_part+_size, _part+total_filled, part_buffer);
	part_filled = total_filled - _size;

	return true;
}

bool reader::getFileReader(int mode, const char* fastq_file, reader*
																	file_reader) {
	uint64_t gzip_buffer_size = 1ULL << 26;
	uint64_t bzip2_buffer_size = 1ULL << 26;

	if (mode == 0) {
		if ((file_reader->in = fopen(fastq_file, "rb")) == NULL)
			return false;
	} else if (mode == 1) {
		if ((file_reader->in_gzip = gzopen(fastq_file, "rb")) == NULL)
			return false;
		gzbuffer(file_reader->in_gzip, gzip_buffer_size);
	} else if (mode == 2) {
		file_reader->in = fopen(fastq_file, "rb");
		if (!file_reader->in)
			return false;
		setvbuf(file_reader->in, NULL, _IOFBF, bzip2_buffer_size);
		if ((file_reader->in_bzip2 = BZ2_bzReadOpen(&file_reader->bzerror,
																								file_reader->in, 0, 0, NULL,
																								0)) == NULL) {
			fclose(file_reader->in);
			return false;
		}
	}
	return true;
}


#endif 
