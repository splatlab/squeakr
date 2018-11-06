//
// Copyright (C) 2017 Prahsant Pandey, Michael Bender,
// Rob Johnson, Rob Patro
//
// This file is part of Squeakr.
//

#ifndef __SQUEAKR_FILESYSTEM_HPP__
#define __SQUEAKR_FILESYSTEM_HPP__

#include <vector>
#include <string>

namespace squeakr {
	namespace fs {
		std::string GetDir(std::string str);
		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool FileExists(const char* path);
		// Taken from
		// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
		bool DirExists(const char* path);
		void MakeDir(const char* path);
		// Taken from
		// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
		std::vector<std::string> GetFilesExt(const char *dir, const char *ext);
	}
}

#endif //__SQUEAKR_FILESYSTEM_HPP__
