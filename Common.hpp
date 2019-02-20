/**************************************************
*
* Common.hpp
*
* Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
*
* This file is part of GapCloser.
*
* GapCloser is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* GapCloser is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with GapCloser. If not, see <http://www.gnu.org/licenses/>.
*
**************************************************/

#ifndef COMMON_HPP_
#define COMMON_HPP_



#ifndef NULL
#define NULL    0
#endif

typedef unsigned long long Number_t;
typedef unsigned int Len_t;
typedef unsigned char Short_Len_t;

//static const unsigned int MAX_STRING_LEN = 256;
static const unsigned int MAX_STRING_LEN = 2048;
static const unsigned int DATA_LEN = 32;
static const unsigned int LEN_TO_BITLEN = 2;

int Debug = 0;
Len_t maxReadLength=100;

int NNumber = 1;

int sameKmer = 0;	
int sameEnd = 0;

#endif /*COMMON_HPP_*/
