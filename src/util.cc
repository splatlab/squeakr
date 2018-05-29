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
												timeval* end)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	std::cout << desc << "Total Time Elapsed: " << std::to_string(time_elapsed)
		<< "seconds" << std::endl;
}

/* check if it's the end of the file. */
inline bool is_eof(reader &file_reader, int mode) {
	if (mode == 0)
		return feof(file_reader.in) != 0;
	else if (mode == 1)
		return gzeof(file_reader.in_gzip) != 0;
	else if (mode == 2)
		return file_reader.bzerror == BZ_STREAM_END;

	return true;
}

/* move the pointer to the end of the next newline. */
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos) {
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
static bool fastq_read_parts(int mode, file_pointer *fp) {
	char *& _part = (fp->part);
	uint64_t& _size = fp->size;
	char*& part_buffer = (fp->part_buffer);
	uint64_t& part_filled = fp->part_filled;
	reader& file_reader = *(fp->freader.get());

	uint32_t OVERHEAD_SIZE = 65535;
	uint64_t part_size = 1ULL << 23;
	char *part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
	memcpy(part, part_buffer, part_filled);

	if(is_eof(file_reader, mode))
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
	if(is_eof(file_reader, mode))
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

bool getFileReader(int mode, const char* fastq_file, reader* file_reader) {
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

std::string last_part(std::string str, char c) {
	uint64_t found = str.find_last_of(c);
	return str.substr(found + 1);
}

