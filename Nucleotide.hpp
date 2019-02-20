/**************************************************
*
* Nucleotide.hpp
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

#ifndef NUCLEOTIDE_HPP_
#define NUCLEOTIDE_HPP_

/* 
define operations on nucleotide level, including
    1) convertion between base and number.
    2) comparison between sequences.
    3) reverse complement of base/number/sequence.
 */

typedef char Nucleotide_t;
static const Nucleotide_t Adenine = 0;
static const Nucleotide_t Cytosine = 1;
static const Nucleotide_t Guanine = 2;
static const Nucleotide_t Thymine = 3;

const int alphabet[128] =
{
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

const char bases[4] = {'A', 'C', 'G', 'T'};


template<typename T>
int compare(T const& leftNumber, T const& rightNumber) {
	
	if (leftNumber < rightNumber)
		return -1;
	else if (leftNumber > rightNumber)
		return 1;
	else
		return 0;
}

template<typename T>
int compare(T* leftNumbers, Len_t leftLen, T* rightNumbers, Len_t rightLen) {
	
	//if same length then compare data
	if (leftLen == rightLen) {
		
		Len_t numLen = leftLen/DATA_LEN;
		if (leftLen%DATA_LEN) numLen++;
		int result = 0;
		for (Len_t i=0; i<numLen; i++) {
			
			if (leftNumbers[i] < rightNumbers[i]) {
				result = -1;
				break;
			}
			else if (leftNumbers[i] > rightNumbers[i]) {
				result = 1;
				break;
			}
		}
		return result;
	}
	else if (leftLen < rightLen) {
		return -1;
	}
	else {
		return 1;
	}
}

//Reverses a word into its Watson-Crick reverse complement
inline Number_t reverseComplement(Number_t number, Len_t len) {
	Number_t numReverse = 0;
	Nucleotide_t nucleotide;

	for (Len_t i=0; i<len; ++i) {
		nucleotide = number & 3;
		numReverse <<= 2;
		numReverse += 3 - nucleotide;
		number >>= 2;
	}
	
	return numReverse;
}

template<typename T>
inline void reverseComplement(T* numbers, T*& numbersReverse, Len_t numLen, Len_t len, Len_t start=0) {
	
	numbersReverse = new T[numLen];
	
	//clear
	for (Len_t i=0; i<numLen; ++i) {
		numbersReverse[i] = 0;
	}
	
	//copy to temp array
	T* numbersTemp = new T[numLen];
	for (Len_t i=0; i<numLen; ++i) {
		numbersTemp[i] =  numbers[i];
	}
	
	Nucleotide_t nucleotide;
	for (int i=len-1; i>=(int)start; --i) {
		Len_t index = i/DATA_LEN;
		Len_t indexReverse = (len-1-i)/DATA_LEN;
		nucleotide = numbersTemp[index] & 3;
		numbersReverse[indexReverse] <<= 2;
		numbersReverse[indexReverse] += 3 - nucleotide;
		numbersTemp[index] >>= 2;
	}
	
	delete [] numbersTemp;
}

inline char* reverseComplement(char const* sequence, Len_t len) {

	char* reverseSeq = new char[len+1];
	for (Len_t i=0; i<len; i++) {
		
		int j=len-i-1;
		switch (toupper(sequence[j])) {
		case 'A':
			reverseSeq[i]='T';
			break;
		case 'C':
			reverseSeq[i]='G';
			break;
		case 'G':
			reverseSeq[i]='C';
			break;
		case 'T':
			reverseSeq[i]='A';
			break;
		case 'N':
			reverseSeq[i]='N';
			break;
		default:
			reverseSeq[i]=toupper(sequence[j]);
			break;
		}
	}
	reverseSeq[len] = '\0';
	return reverseSeq;
}

inline int nucleotideToNumber(Nucleotide_t nucleotide) {
	return alphabet[nucleotide];
/*	
	int nucleotideNum;

	switch (nucleotide) {
	case 'A':
		nucleotideNum = Adenine;
		break;
	case 'C':
		nucleotideNum = Cytosine;
		break;
	case 'G':
		nucleotideNum = Guanine;
		break;
	case 'T':
		nucleotideNum = Thymine;
		break;
	case 'a':
		nucleotideNum = Adenine;
		break;
	case 'c':
		nucleotideNum = Cytosine;
		break;
	case 'g':
		nucleotideNum = Guanine;
		break;
	case 't':
		nucleotideNum = Thymine;
		break;
	default:
		nucleotideNum = Adenine;  //N to A
	}
	return nucleotideNum;
*/
}

inline Nucleotide_t numberToNucleotide(int nucleotideNum) {
	return bases[nucleotideNum];
/*	
	Nucleotide_t nucleotide;
	
	switch (nucleotideNum) {
	case Adenine:
		nucleotide = 'A';
		break;
	case Cytosine:
		nucleotide = 'C';
		break;
	case Guanine:
		nucleotide = 'G';
		break;
	case Thymine:
		nucleotide = 'T';
		break;
	default:
		nucleotide = 'A';
	}
	return nucleotide;
*/
}

inline Number_t sequenceToNumber(char const* sequence, Len_t len) {
	
	Number_t number = 0;
	for (Len_t i=0; i<len; ++i) {
		
		number = number*4 + nucleotideToNumber(sequence[i]);
	}
	return number;
}

inline void numberToSequence(Number_t number, char* sequence, Len_t len) {
	
	
	for (int i=len-1; i>=0; --i) {
		
		sequence[i] = numberToNucleotide(number%4);
		number /= 4;
	}
}

template<typename T>
inline void sequenceToNumbers(char const* sequence, Len_t len, T*& numbers, Len_t& numLen) {
	
	if (len%DATA_LEN) { //len%32
		numLen = len/DATA_LEN+1;
	}
	else {
		numLen = len/DATA_LEN;
	}
	numbers = new T[numLen];
	
	//clear
	for (Len_t i=0; i<numLen; ++i) {
		numbers[i] = 0;
	}
	
	for (Len_t i=0; i<len; ++i) {
		
		numbers[i/DATA_LEN] = numbers[i/DATA_LEN]*4 + nucleotideToNumber(sequence[i]);
	}
}

template<typename T>
inline void numbersToSequence(T* numbers, Len_t numLen, char* sequence, Len_t len) {
	
//	Len_t numLen;
//	if (len%DATA_LEN) {
//		numLen = len/DATA_LEN+1;
//	}
//	else {
//		numLen = len/DATA_LEN;
//	}
	//copy to temp array
	T* numbersTemp = new T[numLen];
	for (Len_t i=0; i<numLen; ++i) {
		numbersTemp[i] =  numbers[i];
	}
	
	for (int i=len-1; i>=0; --i) {
		
		sequence[i] = numberToNucleotide(numbersTemp[i/DATA_LEN]%4);
		numbersTemp[i/DATA_LEN] /= 4;
	}
	
	delete [] numbersTemp;
}

#endif /*NUCLEOTIDE_HPP_*/
