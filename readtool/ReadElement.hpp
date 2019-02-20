/**************************************************
*
* ReadElement.hpp
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

#ifndef READELEMENT_HPP_
#define READELEMENT_HPP_

class ReadElement
{
protected:
	
	static Read* read;
	
	Len_t len;
	Number_t id;
	Len_t depth;
	Len_t* serialNums;
	
	bool forwardFlag;
	
public:
	
	ReadElement(Len_t* _serialNums, bool _forwardFlag) : 
		len(read[_serialNums[0]].getLen()), 
		id(read[_serialNums[0]].getID()), 
		depth(read[_serialNums[0]].getDepth()), 
		serialNums(_serialNums), 
		forwardFlag(_forwardFlag)
	{
	}
	
	static void staticInitialize(Read* _read) {
		read = _read;
	}
	
	Read* getRawReads() const { return &read[serialNums[0]]; }
	Len_t getLen() const { return len; }
	
	bool isFlag() const {
		
		return read[serialNums[0]].isFlag();
	}
	
	void setFlag() {
		
		if (!isFlag()) {
			for(Len_t i=0; i<depth; i++) {
				read[serialNums[i]].setFlag();
			}
		}
	}
	
	void clearFlag() {
		
		if (isFlag()) {
			for(Len_t i=0; i<depth; i++) {
				read[serialNums[i]].clearFlag();
			}
		}
	}
	
	Number_t getID() const { return id; }
	Len_t getDepth() const { return depth; }
	
	Len_t* getSerialNums() const { return serialNums; }
	
	void getSequence(char* sequence) const {
		
		getSequence(sequence, len, 0);
	}
	
	void getSequence(char* sequence, Len_t len, Len_t start) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getSequence(sequence, len, start);
		}
		else {
			read[serialNums[0]].getSequenceReverse(sequence, len, start);
		}
	}
	
	void getSequenceReverse(char* sequenceReverse) const {
		
		getSequenceReverse(sequenceReverse, len, 0);
	}
	
	void getSequenceReverse(char* sequenceReverse, Len_t len, Len_t start) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getSequenceReverse(sequenceReverse, len, start);
		}
		else {
			read[serialNums[0]].getSequence(sequenceReverse, len, start);
		}
	}
	
	void getSequence(TightString& sequence) const {
		
		getSequence(sequence, len, 0);
	}
	
	void getSequence(TightString& sequence, Len_t len, Len_t start) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getSequence(sequence, len, start);
		}
		else {
			read[serialNums[0]].getSequenceReverse(sequence, len, start);
		}
	}
	
	void getSequenceReverse(TightString& sequenceReverse) const {
		
		getSequenceReverse(sequenceReverse, len, 0);
	}
	
	void getSequenceReverse(TightString& sequenceReverse, Len_t len, Len_t start) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getSequenceReverse(sequenceReverse, len, start);
		}
		else {
			read[serialNums[0]].getSequence(sequenceReverse, len, start);
		}
	}
	
	void getReadData(Number_t*& readData, Len_t& arraySize) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getReadData(readData, arraySize);
		}
		else {
			read[serialNums[0]].getReadDataReverse(readData, arraySize);
		}
	}
	
	void getReadDataReverse(Number_t*& readDataReverse, Len_t& arraySize) const {
		
		if (forwardFlag) {
			read[serialNums[0]].getReadDataReverse(readDataReverse, arraySize);
		}
		else {
			read[serialNums[0]].getReadData(readDataReverse, arraySize);
		}
	}
	
};

Read* ReadElement::read=NULL;



#endif /*READELEMENT_HPP_*/
