/**************************************************
*
* Read.hpp
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

#ifndef READ_HPP_
#define READ_HPP_

#include "Nucleotide.hpp"
#include "TightString.hpp"
using namespace std;


/*
#if defined B3
#define BYTE8_3
#elif defined B4
#define BYTE8_4
#elif defined B5
#define BYTE8_5
#elif defined B6
#define BYTE8_6
#else
//#define BYTE8_3
#define BYTE8_5
#endif
*/

/* data structure and operation on a read */
class Read
{
public:
/*	
#if defined BYTE8_3
	static const Len_t BYTE8_NUM = 3;
	static const Len_t DATA_ARRAY_SIZE = 1;
#elif defined BYTE8_4
	static const Len_t BYTE8_NUM = 4;
	static const Len_t DATA_ARRAY_SIZE = 2;
#elif defined BYTE8_5
	static const Len_t BYTE8_NUM = 5;
	static const Len_t DATA_ARRAY_SIZE = 3;
#elif defined BYTE8_6
	static const Len_t BYTE8_NUM = 6;
	static const Len_t DATA_ARRAY_SIZE = 4;
#endif
*/

	static const Len_t BYTE8_NUM = 6;  // total BYTE8 number of one read
	static const Len_t DATA_ARRAY_SIZE=4;  // BYTE8 number of major data array
	
	static const Len_t DATA_REMAIN_ARRAY_SIZE = 1; //BYTE8 number of extra array
	static const Len_t ARRAY_SIZE = DATA_ARRAY_SIZE + DATA_REMAIN_ARRAY_SIZE;  // BYTE8 number of all data array

	
	static const Len_t BITLEN_DATA_REMAIN = 54;  // bit length of part of read in DATA_REMAIN_ARRAY_SIZE
	static const Len_t BITLEN_LENGTH = 9;  // record read length
	static const Len_t BITLEN_FLAG = 1;


	static const Len_t BITLEN_ID = 34;  // read ID
	static const Len_t BITLEN_DEPTH = 30;  // read depth
	
//	static const Len_t DATA_MAXLEN = (DATA_ARRAY_SIZE*bitsizeof(Number_t) + BITLEN_DATA_REMAIN) / 2;
	static Len_t DATA_MAXLEN;
	
protected:
	
	//the shift length in bit
	static const Len_t SHIFTLEN_DATA_REMAIN = BITLEN_LENGTH+BITLEN_FLAG;
	static const Len_t SHIFTLEN_LENGTH = BITLEN_FLAG;
	static const Len_t SHIFTLEN_ID = BITLEN_DEPTH;
	
	// Binary 11111111 11111111 11111111 11111111 11111111 11111111 11111100 00000000
	static const Number_t MASK_REMAIN_DATA_REMAIN = ~((Number_t)((1 << SHIFTLEN_DATA_REMAIN) - 1));
	// Binary 00000000 00000000 00000000 00000000 00000000 00000000 00000011 11111110
	static const Number_t MASK_REMAIN_LENGTH = ((Number_t)((1 << BITLEN_LENGTH) - 1) << SHIFTLEN_LENGTH);
	// Binary 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000001
	static const Number_t MASK_REMAIN_FLAG = ((Number_t) 1);
	// Binary 11111111 11111111 11111111 11111111 11000000 00000000 00000000 00000000
	static const Number_t MASK_REMAIN_ID = ~((Number_t)((1 << SHIFTLEN_ID) - 1));
	// Binary 00000000 00000000 00000000 00000000 00111111 11111111 11111111 11111111
	static const Number_t MASK_REMAIN_DEPTH = ((Number_t)((1 << BITLEN_DEPTH) - 1));
	
	// Binary 00000000 00000000 00000000 00000000 00000000 00000000 00000011 11111111
	static const Number_t MASK_CLEAR_DATA_REMAIN = ~MASK_REMAIN_DATA_REMAIN;
	// Binary 11111111 11111111 11111111 11111111 11111111 11111111 11111100 00000001
	static const Number_t MASK_CLEAR_LENGTH = ~MASK_REMAIN_LENGTH;
	// Binary 11111111 11111111 11111111 11111111 11111111 11111111 11111111 11111110
	static const Number_t MASK_CLEAR_FLAG = ~MASK_REMAIN_FLAG;
	// Binary 00000000 00000000 00000000 00000000 00111111 11111111 11111111 11111111
	static const Number_t MASK_CLEAR_ID = ~MASK_REMAIN_ID;
	// Binary 11111111 11111111 11111111 11111111 11000000 00000000 00000000 00000000
	static const Number_t MASK_CLEAR_DEPTH = ~MASK_REMAIN_DEPTH;
	
protected:
	
	Number_t data[ARRAY_SIZE];
	Number_t id_depth;  // record read ID and depth
	
public:
	Read() : 
		id_depth(0)
	{
		for (Len_t i=0; i<ARRAY_SIZE; i++) {
			data[i] = 0;
		}
	}

	Read(char* sequence, Len_t len) : 
		id_depth(0)
	{
		initialize(sequence, len);
	}

	Read(Number_t* numbers, Len_t len) : 
		id_depth(0)
	{
		initialize(numbers, len);
	}
	
	Read(Read& read) :
		id_depth(read.id_depth)
	{		
		for (Len_t i=0; i<ARRAY_SIZE; i++) {
			data[i] = read.data[i];
		}
	}
	
	Read& operator= (Read const& read) {

		for (Len_t i=0; i<ARRAY_SIZE; i++) {
			data[i] = read.data[i];
		}
		
		id_depth = read.id_depth;
		
		return *this;
	}
	
	~Read()
	{

	}
	
	void initialize(char* sequence, Len_t len) {
		
		Number_t* numbers;
		Len_t numLen;
		sequenceToNumbers(sequence, len, numbers, numLen);
		setReadData(numbers, len);
		delete [] numbers;
	}
	
	void initialize(Number_t* numbers, Len_t len) {
		
		setReadData(numbers, len);

	}
	
public:
	
	void setReadData(Number_t* readData, Len_t len) {
		
		Len_t arraySize = len/DATA_LEN;
		Len_t remainLen = len%DATA_LEN;
		
		Len_t i=0;
		for (; i<arraySize; i++) {
			data[i] = readData[i];
		}
		
		//set remain sequence
		if (remainLen > 0) {
			if ( arraySize >= DATA_ARRAY_SIZE ) {
				setReadDataRemain(readData[i], remainLen);
			}
			else {
				data[i] = getLeftAlignData(readData[i], remainLen);
			}
		}
		
		setLen(len);
	}
	
	void getReadData(Number_t*& readData, Len_t& arraySize, Len_t len) const {
		
		getReadData(readData, arraySize);
		
		Short_Len_t readLen = getLen();
		//right shift readData for readLen-len bits
		if (readLen > len) {
			for (int i=readLen-1; i>=(int)len; --i) {
				readData[i/DATA_LEN] /= 4;
			}
		}
	}

	/* Copy data in length of 'arraySize' from 'data' to 'readData'. It's client's responsibity 
	    to delete 'readData' when it is done. */
	void getReadData(Number_t*& readData, Len_t& arraySize) const {
		
		Short_Len_t len = getLen();
		
		arraySize = len/DATA_LEN;
		Len_t remainLen = len%DATA_LEN;
		
		if (remainLen == 0) {
			readData = new Number_t[arraySize];
		}
		else {
			readData = new Number_t[arraySize+1];
		}
		
		Len_t i=0;
		for (; i<arraySize; i++) {
			readData[i] = data[i];
		}
		
		if (remainLen > 0) {
			
			readData[i] = getRightAlignData(data[i], remainLen);
			arraySize++;
		}
	}
	
	void getReadDataReverse(Number_t*& readDataReverse, Len_t& arraySize) const {
		
		Number_t* readData;
		getReadData(readData, arraySize);
		
		reverseComplement(readData, readDataReverse, arraySize, getLen());
		delete [] readData;
	}
	
	void getSequence(char* string, Number_t* readData, Len_t len, Len_t start=0) const {
		
		Short_Len_t readLen = getLen();
		//right shift readData for readLen-len-start bits
		if (readLen > len+start) {
			for (int i=readLen-1; i>=(int)(len+start); --i) {
				readData[i/DATA_LEN] /= 4;
			}
		}
		else if (readLen < len+start) {
			len = readLen-start;
		}
		
		for (int i=len+start-1,j=len-1; i>=(int)start; --i,--j) {
			
			string[j] = numberToNucleotide(readData[i/DATA_LEN]%4);
			readData[i/DATA_LEN] /= 4;
		}
		
		string[len] = '\0';
	}
	
	void getSequence(char* string, Len_t len, Len_t start=0) const {
		
		Number_t* readData;
		Len_t arraySize;
		getReadData(readData, arraySize);
		
		getSequence(string, readData, len, start);
		
		delete [] readData;
	}
	
	void getSequenceReverse(char* string, Len_t len, Len_t start=0) const {
		
		Number_t* readDataReverse;
		Len_t arraySize;
		getReadDataReverse(readDataReverse, arraySize);
		
		getSequence(string, readDataReverse, len, start);
		
		delete [] readDataReverse;
	}
	
	void getSequence(TightString& string, Number_t* readData, Len_t len, Len_t start=0) const {
		
		Short_Len_t readLen = getLen();
		//right shift readData for readLen-len-start bits
		if (readLen > len+start) {
			for (int i=readLen-1; i>=(int)(len+start); --i) {
				readData[i/DATA_LEN] /= 4;
			}
		}
		else if (readLen < len+start) {
			len = readLen-start;
		}
		
		for (int i=len+start-1,j=len-1; i>=(int)start; --i,--j) {
			
			string.writeNucleotideAtPosition(readData[i/DATA_LEN]%4, j);
			readData[i/DATA_LEN] /= 4;
		}
	}
	
	void getSequence(TightString& string, Len_t len, Len_t start=0) const {
		
		Number_t* readData;
		Len_t arraySize;
		getReadData(readData, arraySize);
		
		getSequence(string, readData, len, start);
		
		delete [] readData;
	}
	
	void getSequenceReverse(TightString& string, Len_t len, Len_t start=0) const {
		
		Number_t* readDataReverse;
		Len_t arraySize;
		getReadDataReverse(readDataReverse, arraySize);
		
		getSequence(string, readDataReverse, len, start);
		
		delete [] readDataReverse;
	}
	
	void reverse() {
		
		Number_t* readDataReverse;
		Len_t arraySize;
		getReadDataReverse(readDataReverse, arraySize);
		
		setReadData(readDataReverse, getLen());
		
		delete [] readDataReverse;
	}
	
	Number_t getLeftAlignData(Number_t _data, Len_t len) const {
		
		return ( _data << (DATA_LEN - len)*LEN_TO_BITLEN );
	}
	
	Number_t getRightAlignData(Number_t _data, Len_t len) const {
		
		return ( _data >> (DATA_LEN - len)*LEN_TO_BITLEN );
	}
	
	Len_t getArraySize() const {
		
		Len_t arraySize;
		Short_Len_t len = getLen();
		if (len%DATA_LEN) {
			arraySize = len/DATA_LEN+1;
		}
		else {
			arraySize = len/DATA_LEN;
		}
		return arraySize;
	}
	
	void setReadDataRemain(Number_t remain, Len_t remainLen) {
		
		data[DATA_ARRAY_SIZE] &= MASK_CLEAR_DATA_REMAIN;
		data[DATA_ARRAY_SIZE] |= getLeftAlignData(remain, remainLen);
	}
	
	Number_t getReadDataRemain(Len_t remainLen) const {
		
		return getRightAlignData(data[DATA_ARRAY_SIZE], remainLen);
	}
	
	void setLen(Short_Len_t len) {
		
		data[DATA_ARRAY_SIZE] &= MASK_CLEAR_LENGTH;
		data[DATA_ARRAY_SIZE] |= len << SHIFTLEN_LENGTH;
	}
	
	Short_Len_t getLen() const {
		
		return (data[DATA_ARRAY_SIZE] & MASK_REMAIN_LENGTH) >> SHIFTLEN_LENGTH;
	}
	
	bool isFlag() const {
		
		return data[DATA_ARRAY_SIZE] & MASK_REMAIN_FLAG;;
	}
	
	void setFlag() {
		
		data[DATA_ARRAY_SIZE] |= MASK_REMAIN_FLAG;
	}
	
	void clearFlag() {
		
		data[DATA_ARRAY_SIZE] &= MASK_CLEAR_FLAG;
	}
	
	void setID(Number_t id) {
		
		id_depth &= MASK_CLEAR_ID;
		id_depth |= id << SHIFTLEN_ID;
	}
	
	Number_t getID() const {
		
		return id_depth >> SHIFTLEN_ID;
	}
	
	void setDepth(Len_t depth) {
		
		id_depth &= MASK_CLEAR_DEPTH;
		id_depth |= depth;
	}
	
	Len_t getDepth() const {
		
		return id_depth & MASK_REMAIN_DEPTH;
	}
	
	int compare(Number_t* rightNumbers, Len_t rightLen) const {
		
		Number_t* readData;
		Len_t arraySize;
		getReadData(readData, arraySize);
		int result = ::compare(readData, getLen(), rightNumbers, rightLen);
		delete [] readData;
		return result;
	}
	
	int compare(Len_t leftLen, Number_t* rightNumbers, Len_t rightLen) const {
		
		Number_t* readData;
		Len_t arraySize;
		getReadData(readData, arraySize, leftLen);
		int result = ::compare(readData, leftLen, rightNumbers, rightLen);
		delete [] readData;
		return result;
	}
	
	int compare(Read const& read) const {
		
		Number_t const* rightData = read.data;
		int result = 0;
		for (Len_t i=0; i<ARRAY_SIZE; i++) {
			
			if (data[i] < rightData[i]) {
				result = -1;
				break;
			}
			else if (data[i] > rightData[i]) {
				result = 1;
				break;
			}
		}
		return result;
	}
	/*
	int compareReverse(Read const& read) const {
		
		Number_t* leftReadDataReverse;
		getReadDataReverseLeftAlign(leftReadDataReverse);
		
		Number_t* rightReadDataReverse;
		read.getReadDataReverseLeftAlign(rightReadDataReverse);
		
		int result = ::compare(leftReadDataReverse, DATA_LEN*(DATA_ARRAY_SIZE+1), rightReadDataReverse, DATA_LEN*(DATA_ARRAY_SIZE+1));
		
		delete [] leftReadDataReverse;
		delete [] rightReadDataReverse;
		
		return result;
	}*/
	
	bool operator ==(Read const& right) {
		return compare(right) == 0;
	}
	
	bool operator !=(Read const& right) {
		return compare(right) != 0;
	}
	
	bool operator <=(Read const& right) {
		return compare(right) <= 0;
	}
	
	bool operator <(Read const& right) {
		return compare(right) < 0;
	}
	
	bool operator >=(Read const& right) {
		return compare(right) >= 0;
	}
	
	bool operator >(Read const& right) {
		return compare(right) > 0;
	}
};

int compare(Read const& leftRead, Read const& rightRead) {
	
	return leftRead.compare(rightRead);
}
/*
int compareReverse(Read const& leftReadReverse, Read const& rightReadReverse) {
	
	return leftReadReverse.compareReverse(rightReadReverse);
}*/



#endif /*READ_HPP_*/
