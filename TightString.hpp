/**************************************************
*
* TightString.hpp
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

#ifndef TIGHTSTRING_HPP_
#define TIGHTSTRING_HPP_

#include "Nucleotide.hpp"
#include "ArrayBlock.hpp"


class TightString
{
public:
	typedef unsigned char Codon;
	typedef unsigned int Coordinate;

	// Binary 11000000
	static const Codon MASK_REMAIN_I0 = (((Codon) 3) << 6);
	// Binary 00110000
	static const Codon MASK_REMAIN_I1 = (((Codon) 3) << 4);
	// Binary 00001100
	static const Codon MASK_REMAIN_I2 = (((Codon) 3) << 2);
	// Binary 00000011
	static const Codon MASK_REMAIN_I3 = ((Codon) 3);
	
	// Binary 00111111
	static const Codon MASK_CLEAR_I0 = (Codon)~MASK_REMAIN_I0;
	// Binary 11001111
	static const Codon MASK_CLEAR_I1 = ~MASK_REMAIN_I1;
	// Binary 11110011
	static const Codon MASK_CLEAR_I2 = ~MASK_REMAIN_I2;
	// Binary 11111100
	static const Codon MASK_CLEAR_I3 = ~MASK_REMAIN_I3;
	
	//250*4=1K
	static const Len_t BLOCK_LEN = 250;
	
protected:

	Len_t blockLen;
	ArrayBlock<Codon>* sequence;
	Coordinate length;
public:	
	
//	TightString(Len_t _blockLen=0) {
//		
//		initialize(0, _blockLen);
//	}

	TightString(Coordinate len=0, Len_t _blockLen=0) {
		
		initialize(len, _blockLen);
	}
	
	TightString(Number_t number, Coordinate len, Len_t _blockLen=0) {
		
		initialize(len, _blockLen);
		copyNumber(number, len);
	}
	
	TightString(Number_t* numbers, Len_t numLen, Coordinate len, Len_t _blockLen=0) {
		
		initialize(len, _blockLen);
		copyNumbers(numbers, numLen, len);
	}
	
	TightString(char const* sequence, Coordinate len, Len_t _blockLen=0) {
		
		initialize(len, _blockLen);
		copyString(sequence, len);
	}
	
	TightString(const TightString& tStr, Len_t _blockLen=0) {

		initialize(tStr.length, _blockLen);
		Coordinate arrayLength = getArrayLength();
		Coordinate i;
		for (i=0; i<arrayLength; ++i) {
			
			(*sequence)[i] = (*(tStr.sequence))[i];
		}
	}
	
	TightString(const TightString& lStr, const TightString& rStr, Len_t _blockLen=0) {
		
		initialize(lStr.length + rStr.length, _blockLen);
		concatenate(lStr, rStr);
	}

	~TightString() {
		
		destroy();
	}






	void initialize(Coordinate len, Len_t _blockLen=0) {

		length = len;
		Coordinate arrayLength = getArrayLength();
		
		if (arrayLength<=0) {
			arrayLength = 1;
		}
		
		Len_t sequenceLen;
		if (_blockLen)
			blockLen = _blockLen;
		else
			blockLen = arrayLength;
		
		if (arrayLength%blockLen) {
			sequenceLen = ((arrayLength/blockLen)+1)*blockLen;
		}
		else {
			sequenceLen = arrayLength;
		}
		
		sequence = new ArrayBlock<Codon>(blockLen, sequenceLen);
		clearBlock(0, sequenceLen, *sequence);
	}

	void destroy() {
		delete sequence;
	}
	
	Coordinate getLength() const { return length; }
	
	void setLength(Coordinate _length) { length = _length; }
	
	Coordinate getArrayLength() const {
		
		Coordinate arrayLength = length / 4;
		if (length % 4 > 0)
			++arrayLength;
		return arrayLength;
	}
	
	void copyNumber(Number_t number, Coordinate len, Coordinate start=0) {
		
		for (int i=len-1+start; i>=(int)start; --i) {
			writeNucleotideNumber(number%4, (*(this->sequence))[i/4], i%4);
			number /= 4;
		}
	}
	
	void copyNumbers(Number_t* numbers, Len_t numLen, Coordinate len, Coordinate start=0) {
		
		//copy to temp array
		Number_t* numbersTemp = new Number_t[numLen];
		for (Len_t i=0; i<numLen; ++i) {
			numbersTemp[i] =  numbers[i];
		}
		
		for (int i=len-1,j=i+start; i>=0; --i,--j) {
			
			writeNucleotideNumber(numbersTemp[i/DATA_LEN]%4, (*(this->sequence))[j/4], j%4);
			numbersTemp[i/DATA_LEN] /= 4;
		}
		
		delete [] numbersTemp;
	}
	
	// copy a tradionnal string of A,T,G, and C of length size
	void copyString(char const* sequence, Coordinate len) {

		for (Coordinate i = 0; i < len; i++)
			writeNucleotide(sequence[i], (*(this->sequence))[i/4], i%4);
	}

	// Adds a number into the codon at the desired position (0, 1, 2, or 3);
	void writeNucleotideNumber(Nucleotide_t nucleotide, Codon& codon, Coordinate position) {

		if (position == 0) {
			codon &= MASK_CLEAR_I0;
			codon += nucleotide << 6;
		} else if (position == 1) {
			codon &= MASK_CLEAR_I1;
			codon += nucleotide << 4;
		} else if (position == 2) {
			codon &= MASK_CLEAR_I2;
			codon += nucleotide << 2;
		} else if (position == 3) {
			codon &= MASK_CLEAR_I3;
			codon += nucleotide;
		}
	}

	// Adds a nucleotide into the codon at the desired position (0, 1, 2, or 3);
	void writeNucleotide(Nucleotide_t nucleotide, Codon& codon, Coordinate position) {
		
		int nucleotideNum = nucleotideToNumber(nucleotide);
		writeNucleotideNumber(nucleotideNum, codon, position);
	}
	
	void writeNucleotideAtPosition(Nucleotide_t nucleotide, Coordinate position) {
		if (position >= length)
			return;

		writeNucleotideNumber(nucleotide, (*sequence)[position/4], position%4);
	}
	
	char* readTightString() const {

		if (sequence == NULL || length == 0)
			return NULL;
		
		Coordinate index, index4;
		Codon codon;
		char* string = new char[length+1];

		for (index = 0; index < length/4; index++) {
			index4 = index << 2;
			codon = (*sequence)[index];
			string[index4] = numberToNucleotide((codon & MASK_REMAIN_I0) >> 6);
			string[index4 + 1] = numberToNucleotide((codon & MASK_REMAIN_I1) >> 4);
			string[index4 + 2] = numberToNucleotide((codon & MASK_REMAIN_I2) >> 2);
			string[index4 + 3] = numberToNucleotide(codon & MASK_REMAIN_I3);
		}

		if (length % 4) {
			index4 = index << 2;
			codon = (*sequence)[index];
	
			switch ((length-1) % 4) {
			case 3:
				string[index4 + 3] = numberToNucleotide(codon & MASK_REMAIN_I3);
			case 2:
				string[index4 + 2] = numberToNucleotide((codon & MASK_REMAIN_I2) >> 2);
			case 1:
				string[index4 + 1] = numberToNucleotide((codon & MASK_REMAIN_I1) >> 4);
			case 0:
				string[index4] = numberToNucleotide((codon & MASK_REMAIN_I0) >> 6);
			}
		}

		string[length] = '\0';

		return string;
	}
	
	void readTightString(Number_t& number) const {
		
		number = 0;
		Coordinate arrayLength = getArrayLength();
		Coordinate i;
		for (i=0; i<arrayLength; ++i) {
			
			number |= (*sequence)[i];
			number <<= 8;
		}
		number >>= 8;
		//right alignment
		if (length%4)
			number >>= 8 - (length%4)*2;
	}
	
	void readTightStringFragment(Coordinate start, Coordinate finish, char *string) const {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;

		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;

		for (index = 0; index < inFinish - inStart; index++)
			string[index] = numberToNucleotide(getNucleotide(index + inStart));

		string[inFinish - inStart] = '\0';
	}
	
	void writeTightStringFragment(Coordinate start, Coordinate finish, TightString const& string) {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;
		
		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;
		
		for (index = 0; index < inFinish - inStart; index++) {
			
			writeNucleotideAtPosition(string.getNucleotide(index), index + inStart);
		}
	}
	
	void readTightStringFragment(Coordinate start, Coordinate finish, TightString& string) const {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;

		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;
		
		for (index = 0; index < inFinish - inStart; index++) {
			
			string.writeNucleotideAtPosition(getNucleotide(index + inStart), index);
		}
	}
	
	void readTightStringFragment(Coordinate start, Coordinate finish, Number_t& number) const {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;

		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;
		
		number = 0;
		for (index = 0; index < inFinish - inStart; index++) {
			number = number<<2 | getNucleotide(index + inStart);
//			number <<= 2;
		}
//		number >>= 2;
	}
	
	void readTightStringFragment(Coordinate start, Coordinate finish, Number_t*& numbers, Len_t& numLen) const {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;

		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;
		
		Len_t len = inFinish - inStart;
		if (len%DATA_LEN) {
			numLen = len/DATA_LEN+1;
		}
		else {
			numLen = len/DATA_LEN;
		}
		numbers = new Number_t[numLen];
		
		//clear
		for (Len_t i=0; i<numLen; ++i) {
			numbers[i] = 0;
		}
		
		for (index=0; index<len; index++) {
			
			numbers[index/DATA_LEN] = numbers[index/DATA_LEN]*4 + getNucleotide(index + inStart);
		}
	}
	
	void readTightStringFragmentAtReverse(Coordinate start, Coordinate finish, TightString& string) const {
		Coordinate index;
		Coordinate inStart = start;
		Coordinate inFinish = finish;

		if (start > length)
			inStart = finish;
		else if (finish > length)
			inFinish = length;
		
		for (index = 0; index < inFinish - inStart; index++) {
			
			string.writeNucleotideAtPosition(getInverseNucleotide(index + inStart), string.getLength()-index-1);
		}
	}
	
	Nucleotide_t getNucleotide(Coordinate nucleotideIndex) const {
		
		Codon codon = (*sequence)[nucleotideIndex/4];

		switch (nucleotideIndex % 4) {
		case 3:
			return (codon & MASK_REMAIN_I3);
		case 2:
			return (codon & MASK_REMAIN_I2) >> 2;
		case 1:
			return (codon & MASK_REMAIN_I1) >> 4;
		case 0:
			return (codon & MASK_REMAIN_I0) >> 6;
		}

		return '?';
	}
	
	Nucleotide_t getInverseNucleotide(Coordinate nucleotideIndex) const {
		
		Codon codon = (*sequence)[nucleotideIndex/4];

		switch (nucleotideIndex % 4) {
		case 3:
			return (3 - (codon & MASK_REMAIN_I3));
		case 2:
			return (3 - ((codon & MASK_REMAIN_I2) >> 2));
		case 1:
			return (3 - ((codon & MASK_REMAIN_I1) >> 4));
		case 0:
			return (3 - ((codon & MASK_REMAIN_I0) >> 6));
		}

		return '?';
	}
	
	char getNucleotideChar(Coordinate nucleotideIndex) {
		Codon codon = (*sequence)[nucleotideIndex/4];

		switch (nucleotideIndex % 4) {
		case 3:
			return numberToNucleotide(codon & MASK_REMAIN_I3);
		case 2:
			return numberToNucleotide((codon & MASK_REMAIN_I2) >> 2);
		case 1:
			return numberToNucleotide((codon & MASK_REMAIN_I1) >> 4);
		case 0:
			return numberToNucleotide((codon & MASK_REMAIN_I0) >> 6);
		}

		return '?';
	}
	
	
	char getInverseNucleotideChar(Coordinate nucleotideIndex) {
		Codon codon = (*sequence)[nucleotideIndex/4];

		switch (nucleotideIndex % 4) {
		case 3:
			return numberToNucleotide(3 - (codon & MASK_REMAIN_I3));
		case 2:
			return numberToNucleotide(3 - ((codon & MASK_REMAIN_I2) >> 2));
		case 1:
			return numberToNucleotide(3 - ((codon & MASK_REMAIN_I1) >> 4));
		case 0:
			return numberToNucleotide(3 - ((codon & MASK_REMAIN_I0) >> 6));
		}

		return '?';
	}
	
	TightString concatenate(const TightString& rStr) {
		
		if (rStr.length == 0)
			return *this;
		
		TightString unionStr(length + rStr.length);

		//copy left member
		//unit is codon
		Coordinate arrayLength = getArrayLength();
		Coordinate i;
		for (i=0; i<arrayLength; ++i) {
			
			(*(unionStr.sequence))[i] = (*sequence)[i];
		}
		
		//copy right member
		//unit is bit
		for (i=0; i<rStr.length; ++i) {
			
			unionStr.writeNucleotideAtPosition(rStr.getNucleotide(i), length+i);
		}
		
		return unionStr;
	}
	
	void concatenate(const TightString& lStr, const TightString& rStr) {
		
		//copy left member
		//unit is codon
		Coordinate arrayLength = lStr.getArrayLength();
		Coordinate i;
		for (i=0; i<arrayLength; ++i) {
			
			(*sequence)[i] = (*lStr.sequence)[i];
		}
		
		//copy right member
		//unit is bit
		for (i=0; i<rStr.length; ++i) {
			
			writeNucleotideAtPosition(rStr.getNucleotide(i), lStr.length+i);
		}
	}
	
	void increaseBlock() {
		
		Len_t oldLength = (*sequence).getLength();
		(*sequence).setLength(oldLength+blockLen);
		clearBlock(oldLength, oldLength+blockLen, *sequence);
	}
	
	void clearBlock(Len_t iStart, Len_t iEnd, ArrayBlock<Codon>& sequence) {
		
		for (Len_t i=iStart; i<iEnd; i++) {
			sequence[i] = 0;
		}
	}
	
	void append(TightString const& string, Len_t len) {
		
		Coordinate oldLength = length;
		length += len;
		Coordinate arrayLength = getArrayLength();
		
		//check array length
		while (arrayLength >= (*sequence).getLength()) {
			increaseBlock();
		}
		
		for (Len_t i=0; i<len; ++i) {
			writeNucleotideAtPosition(string.getNucleotide(i), oldLength+i);
		}
	}
	
	void append(Number_t number, Len_t len) {
		
		Coordinate oldLength = length;
		length += len;
		Coordinate arrayLength = getArrayLength();
		
		//check array length
		while (arrayLength >= (*sequence).getLength()) {
			increaseBlock();
		}
		
		for (int i=len-1; i>=0; --i) {
			writeNucleotideAtPosition(number%4, oldLength+i);
			number /= 4;
		}
	}
	
	TightString& operator= (const TightString& tStr) {
		
		if (sequence || length) 
			destroy();
		
		initialize(tStr.length);
		Coordinate arrayLength = getArrayLength();
		Coordinate i;
		for (i=0; i<arrayLength; ++i) {
			
			(*sequence)[i] = (*(tStr.sequence))[i];
		}
		return *this;
	}

	
	TightString operator+ (const TightString& rStr) {
		
		return concatenate(rStr);
	}
	
//	TightString& operator+= (const TightString& rStr) {
//		
//		TightString tStr(*this, rStr);
//	}
	
	Nucleotide_t const operator [] (Coordinate nucleotideIndex) const {
		
		return getNucleotide(nucleotideIndex);
	}
	
	Nucleotide_t operator [] (Coordinate nucleotideIndex) {
		
		return getNucleotide(nucleotideIndex);
	}
	
	int compare(TightString const& rStr) const {
		
		Len_t leftLen = getLength();
		Len_t rightLen = rStr.getLength();
		
		if (leftLen < rightLen) {
			return -1;
		}
		else if (leftLen > rightLen) {
			return 1;
		}
		//if same length then compare data
		else {
			
			int result = 0;
			Coordinate arrayLength = getArrayLength();
			for (Coordinate i=0; i<arrayLength; ++i) {
				
				if ((*sequence)[i] < (*(rStr.sequence))[i]) {
					result = -1;
					break;
				}
				else if ((*sequence)[i] > (*(rStr.sequence))[i]) {
					result = 1;
					break;
				}
			}
			
			return result;
		}
	}
	
	bool operator ==(TightString const& right) {
		return compare(right) == 0;
	}

	bool operator !=(TightString const& right) {
		return compare(right) != 0;
	}

	bool operator <=(TightString const& right) {
		return compare(right) <= 0;
	}

	bool operator <(TightString const& right) {
		return compare(right) < 0;
	}

	bool operator >=(TightString const& right) {
		return compare(right) >= 0;
	}

	bool operator >(TightString const& right) {
		return compare(right) > 0;
	}
	
};

#endif /*TIGHTSTRING_HPP_*/
