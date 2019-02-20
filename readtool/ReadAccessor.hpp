/**************************************************
*
* ReadAccessor.hpp
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

#ifndef READACCESSOR_HPP_
#define READACCESSOR_HPP_

#include "readhash/ReadHash.hpp"
#include "ReadElement.hpp"
#include "readtool/PairInfo.hpp"
#include "sort/MedianOfThreeQuickSorter.hpp"

using namespace data_structure;

/* accessor to reads in readHash */
class ReadAccessor
{
public:
	static const Short_Len_t inAll = 1;
	static const Short_Len_t inFlag = 2;
	static const Short_Len_t inNoFlag = 3;
	
protected:
	ReadHash const& readHash;
	PairInfo const& pairInfo;
	
	Short_Len_t overlapParam;
	
public:
	
	ReadAccessor(ReadHash const& _readHash, PairInfo const& _pairInfo) : 
		readHash(_readHash), 
		pairInfo(_pairInfo), 
		overlapParam(_readHash.getOverlapParam())
	{
		ReadElement::staticInitialize(readHash.getReads());
	}
	
	Number_t getReadsSum() const { return readHash.getReadsSum(); }

	// get reads began with 'sequence'
	void getReadsBeginWith(char const* sequence, LinkedList<ReadElement>& readElements, bool isSetFlag=true, Short_Len_t inWhichSet=inNoFlag) {
		
		Len_t sequenceLen = strlen(sequence);
		if (sequenceLen<overlapParam) return;
		
		Number_t overlap = sequenceToNumber(sequence, overlapParam);
		
		getReadsBeginWith(overlap, readElements, isSetFlag, inWhichSet);
	}
	
	void getReadsBeginWith(TightString const& sequence, LinkedList<ReadElement>& readElements, bool isSetFlag=true, Short_Len_t inWhichSet=inNoFlag) {
		
		Len_t sequenceLen = sequence.getLength();
		if (sequenceLen<overlapParam) return;
		
		Number_t overlap;
		sequence.readTightStringFragment(0, overlapParam, overlap);	//binary fromat stored in overlap
		
		getReadsBeginWith(overlap, readElements, isSetFlag, inWhichSet);
	}
	
	void getReadsBeginWith(Number_t overlap, LinkedList<ReadElement>& readElements, bool isSetFlag=true, Short_Len_t inWhichSet=inNoFlag) {
		
		getReads(overlap, readElements, inWhichSet);

		if (isSetFlag) {
			setReadsFlag(readElements);
		}
	}
	
	void setReadsFlag(LinkedList<ReadElement>& readElements) {
		
		ListElement<ReadElement> const* ptrRead;
		for (ptrRead = readElements.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
			
			const_cast<ReadElement&>(ptrRead->getDatum()).setFlag();
		}
	}
	
	void getReads(Number_t overlap, LinkedList<ReadElement>& readElements, Short_Len_t inWhichSet=inNoFlag) {
		
		switch(inWhichSet) {
			case inAll:
				getReadsInAll(overlap, readElements);
				break;
			case inFlag:
				getReadsInFlag(overlap, readElements);
				break;
			case inNoFlag:
				getReadsInNoFlag(overlap, readElements);
				break;
			default:
				getReadsInNoFlag(overlap, readElements);
				break;
		}
	}
	
	void getReadsInAll(Number_t overlap, LinkedList<ReadElement>& readElements) {
		
		if (readHash.getHash().find(overlap)) {
			Number_t readsSum = getReadsSum();
			ReadsHashElement readsHashElement = readHash.getHash()[overlap];

			// forward
			Number_t start = readsHashElement.serialNumForward.start;			
			Len_t* serialNums = readHash.getSerialNumsForward();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], true);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}

			// reversed
			start = readsHashElement.serialNumReverse.start;
			serialNums = readHash.getSerialNumsReverse();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], false);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}
		}
	}
	
	void getReadsInFlag(Number_t overlap, LinkedList<ReadElement>& readElements) {
		
		if (readHash.getHash().find(overlap)) {
			Number_t readsSum = getReadsSum();
			ReadsHashElement readsHashElement = readHash.getHash()[overlap];
			
			Number_t start = readsHashElement.serialNumForward.start;			
			Len_t* serialNums = readHash.getSerialNumsForward();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], true);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					if (readElement.isFlag())
						readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}
			
			start = readsHashElement.serialNumReverse.start;
			serialNums = readHash.getSerialNumsReverse();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], false);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					if (readElement.isFlag())
						readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}
		}
	}
	
	void getReadsInNoFlag(Number_t overlap, LinkedList<ReadElement>& readElements) {
		
		if (readHash.getHash().find(overlap)) {
			Number_t readsSum = getReadsSum();
			ReadsHashElement readsHashElement = readHash.getHash()[overlap];
			
			Number_t start = readsHashElement.serialNumForward.start;			
			Len_t* serialNums = readHash.getSerialNumsForward();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], true);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					if (!readElement.isFlag())
						readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}
			
			start = readsHashElement.serialNumReverse.start;
			serialNums = readHash.getSerialNumsReverse();
			
			for (Len_t i=start; i<readsSum;) {
				
				ReadElement readElement(&serialNums[i], false);
				
				Number_t* readData;
				Len_t arraySize;
				readElement.getReadData(readData, arraySize);
				
				Number_t readOverlap;
				Len_t readLen = readElement.getLen();
				if (readLen < DATA_LEN) {
					readOverlap = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
				}
				else {
					readOverlap = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
				}
				
				delete [] readData;
				
				if (overlap == readOverlap) {
					if (!readElement.isFlag())
						readElements.append(readElement);
					i += readElement.getDepth();
				}
				else {
					break;
				}
			}
		}
	}
	
	void getReadIdsByPair(LinkedList<Number_t>& pairs, LinkedList<Number_t>& ids) {
		
		Len_t pairsSum = pairs.getCount();
		if (pairsSum>0) {
			Number_t* readsID = new Number_t[pairsSum];
			Read* reads = readHash.getReads();
			Len_t i=0;
			ListElement<Number_t> const* ptrPair;
			for (ptrPair = pairs.getHead(); ptrPair != 0; ptrPair = ptrPair->getNext()) {
				
				Number_t iReadPair = ptrPair->getDatum();
				readsID[i] = reads[iReadPair].getID();
				i++;
			}
			
			MedianOfThreeQuickSorter<Number_t> sorter;
			sorter.sort(readsID, pairsSum);
			
			//append read with different readID
			for (i=0; i<pairsSum-1; i++) {
				
				if (readsID[i] != readsID[i+1]) {
					ids.append(readsID[i]);
				}
			}
			ids.append(readsID[i]);
			delete [] readsID;
		}
	}
	
	void getReadIdsByPair(ReadElement const& readElement, LinkedList<Number_t>* ids) {
		
		Len_t arrayLen = pairInfo.getArrayLen();
		LinkedList<Number_t>* pairs = new LinkedList<Number_t>[arrayLen];
		pairInfo.getPairs(readElement, pairs);
		
		for (Len_t i=0; i<arrayLen; i++) {
			
			getReadIdsByPair(pairs[i], ids[i]);
		}
		
		delete [] pairs;
	}
	
	void getReadIdsByPair(Number_t id, LinkedList<Number_t>* ids) {
		
		getReadIdsByPair(getRead(id), ids);
	}

	// function:	 getReadsByPair
	// describtion: keep only one read element for one readID
	// input: 		 pairs -- reads' sirialNum
	// output:	 readElements -- reads with different readID
	void getReadsByPair(LinkedList<Number_t>& pairs, LinkedList<ReadElement>& readElements) {
		
		Len_t pairsSum = pairs.getCount();
		if (pairsSum>0) {
			Number_t* readsID = new Number_t[pairsSum];
			Read* reads = readHash.getReads();
			Len_t i=0;
			ListElement<Number_t> const* ptrPair;
			for (ptrPair = pairs.getHead(); ptrPair != 0; ptrPair = ptrPair->getNext()) {
				
				Number_t iReadPair = ptrPair->getDatum();
				readsID[i] = reads[iReadPair].getID();
				i++;
			}
			
			MedianOfThreeQuickSorter<Number_t> sorter;
			sorter.sort(readsID, pairsSum);
			
			//append read with different readID
			for (i=0; i<pairsSum-1; i++) {
				
				if (readsID[i] != readsID[i+1]) {
					readElements.append(getRead(readsID[i]));
				}
			}
			readElements.append(getRead(readsID[i]));
			delete [] readsID;
		}
	}

	// function:	 getReadsByPair
	// describtion: find paired-end read elemenst of currently konwn read element in libraries with specific insert size
	// input:		 readElement -- currently konwn read element
	//			 insertSize -- specified insert size
	// output: 	 readElements -- found paired-end read elements
	void getReadsByPair(ReadElement const& readElement, LinkedList<ReadElement>& readElements, Len_t insertSize=0) {
		
		LinkedList<Number_t> pairs;
		pairInfo.getPairs(readElement, pairs, insertSize);
		
		getReadsByPair(pairs, readElements);
	}

	// function:	 getReadsByPair
	// describtion: find paired-end read element of currently konwn read element in all libraries
	// input: 		 readElement -- currently konwn read element
	// output:	 readElements -- found paired-end read element
	void getReadsByPair(ReadElement const& readElement, LinkedList<ReadElement>* readElements) {
		
		Len_t arrayLen = pairInfo.getArrayLen();

		// paired-end reads' serialNum
		LinkedList<Number_t>* pairs = new LinkedList<Number_t>[arrayLen];
		pairInfo.getPairs(readElement, pairs);
		
		for (Len_t i=0; i<arrayLen; i++) {
			
			getReadsByPair(pairs[i], readElements[i]);
		}
		
		delete [] pairs;
	}
	
	void setReadFlag(ReadElement& readElement) {
		
		readElement.setFlag();
	}
	
	ReadElement getReadSetFlag(Number_t id) {
		

		ReadElement readElement = getRead(id);
		setReadFlag(readElement);
		return readElement;
	}
	
	ReadElement getRead(Number_t id) {
		
		Len_t* serialNums = &readHash.getSerialNumsForward()[id];
		ReadElement readElement(&serialNums[0], true);
		
		return readElement;
	}
	
	bool isReadFlag(Number_t id) const {
		
		return readHash.getReads()[readHash.getSerialNumsForward()[id]].isFlag();
	}
	
	
	
	class Iterator {
		ReadAccessor const& readAccessor;
		Read* reads;
		Len_t* serialNums;
		Number_t len;
		Number_t sum;
		Number_t indexStart;
//		Number_t indexEnd;
		Number_t index;
		Number_t count;
	public:
		
		Iterator(ReadAccessor const& _readAccessor) : 
			readAccessor(_readAccessor), 
			reads(readAccessor.readHash.getReads()), 
			serialNums(readAccessor.readHash.getSerialNumsForward())
		{
			len = readAccessor.readHash.getReadsSum();
			sum = len;
			indexStart = 0;
			reset();
		}
		
		Iterator(ReadAccessor const& _readAccessor, Number_t _indexStart) : 
			readAccessor(_readAccessor), 
			reads(readAccessor.readHash.getReads()), 
			serialNums(readAccessor.readHash.getSerialNumsForward()) 
		{
			len = readAccessor.readHash.getReadsSum();
			sum = len;
			
			//get next different read
			while(reads[serialNums[_indexStart]].getID() == reads[serialNums[_indexStart+1]].getID()) {
				
				_indexStart = (_indexStart+1)%len;
			}
			_indexStart = (_indexStart+1)%len;
			
			indexStart = _indexStart;
			reset();
		}
		
		Iterator(ReadAccessor const& _readAccessor, Number_t _indexStart, Number_t _indexEnd) : 
			readAccessor(_readAccessor), 
			reads(readAccessor.readHash.getReads()), 
			serialNums(readAccessor.readHash.getSerialNumsForward())
		{
			len = readAccessor.readHash.getReadsSum();
			sum = _indexEnd - _indexStart;
			
			//get next different read
			while(reads[serialNums[_indexStart]].getID() == reads[serialNums[_indexStart+1]].getID()) {
				
				_indexStart = (_indexStart+1)%len;
			}
			_indexStart = (_indexStart+1)%len;
			
			indexStart = _indexStart;
			reset();
		}
		
		bool isDone() const {
			return count >= sum;
		}
			
		ReadElement operator * () const {
			return const_cast<ReadAccessor&>(readAccessor).getReadSetFlag(index);
		}
		
		void operator ++ () {
			
			//iterate reads which not set flag
			do {
				
				Len_t depth = reads[serialNums[index]].getDepth();
				index = (index + depth)%len;
				count += depth;
				
			} while(!isDone() && readAccessor.isReadFlag(index));
		}
		
		void reset() {
			
			index = indexStart;
			count = 0;
		}
	};
	
	Iterator& newIterator() const {
		return *new Iterator(*this);
	}
	
	Iterator& newIterator(Number_t _indexStart) const {
		return *new Iterator(*this, _indexStart);
	}
	
	Iterator& newIterator(Number_t _indexStart, Number_t _indexEnd) const {
		return *new Iterator(*this, _indexStart, _indexEnd);
	}
	
	friend class Iterator;
};



#endif /*READACCESSOR_HPP_*/
