/**************************************************
*
* ArrayBlock.hpp
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

#ifndef ARRAYBLOCK_HPP_
#define ARRAYBLOCK_HPP_

#include "base/LinkedList.hpp"

using namespace data_structure;

/*
define operations on block level, including:
  1) get element in specific position in way of function calling and operator usage.
  2) increase block array's length.
  3) append specific length of blocks to current block array.
  4) combine two block arrays into one.
  5) reverse block array.
*/
template<class T>
class ArrayBlock {
protected:
	
	LinkedList<T*> data;
	unsigned int blockLen;		// blcok number
	unsigned int length;		// sequence length
public:
	
	ArrayBlock& operator =(ArrayBlock const&);
	
	ArrayBlock() :
		data(), blockLen(0), length(0) {
	}
	
	ArrayBlock(unsigned int len) :
		data(), blockLen(len), length(len) {
		
		T* block = new T [len];
		data.append(block);
	}
	
	ArrayBlock(unsigned int _blockLen, unsigned int _length) :
		data(), blockLen(_blockLen), length(_length) {
		
		if (length) {
			unsigned int blockSum = (length-1) / blockLen;
			
			for (unsigned int i=0; i <= blockSum; ++i) {
				data.append(new T [blockLen]);
			}
		}
	}
	
//	ArrayBlock(ArrayBlock<T> const& array) :
//		data(), blockLen(array.blockLen), length(array.length) {
//		
//		ListElement<T*> const* ptr;
//		for (ptr = array.getHead(); ptr != 0; ptr = ptr->getNext())
//			data.append(ptr->getDatum());
//	}
	
	~ArrayBlock() {
		
		ListElement<T*> const* ptr;
		for (ptr = data.getHead(); ptr != 0; ptr = ptr->getNext())
			delete [] ptr->getDatum();
		
		data.purge();
	}

	T& getElement(unsigned int position) const {
		
		if (position >= length)
			throw out_of_range ("invalid position");
//			setLength(length+blockLen);
		
		unsigned int blockIndex = position / blockLen;
		ListElement<T*> const* ptr;
		unsigned int i=0;
		for (ptr = data.getHead(); i!=blockIndex && ptr!=0; ++i, ptr = ptr->getNext());
		
		return ptr->getDatum()[position%blockLen];
	}
	
	T const& operator [](unsigned int position) const {
		return getElement(position);
	}
	
	T& operator [](unsigned int position) {
		return getElement(position);
	}
	
	unsigned int getBlockLen() const {
		return blockLen;
	}
	
	unsigned int getLength() const {
		return length;
	}
	
	void setLength(unsigned int newLength) {
		
		if (newLength) {
			if (length) {
				unsigned int newBlockSum = (newLength-1) / blockLen;
				unsigned int oldBlockSum = (length-1) / blockLen;
				if (newBlockSum > oldBlockSum) {
					
					for (unsigned int i=0; i < newBlockSum-oldBlockSum; ++i) {
						data.append(new T [blockLen]);
					}
				}
			}
			else {
				unsigned int newBlockSum = (newLength-1) / blockLen;
				for (unsigned int i=0; i <= newBlockSum; ++i) {
					data.append(new T [blockLen]);
				}
			}
			length = newLength;
		}
	}
	
	void clear(unsigned int start, unsigned int end) {
		
		for (unsigned int i=start; i<end; i++) {
			getElement(i) = 0;
		}
	}
	
	void increase(unsigned int newLength) {
		
		unsigned int oldLength = getLength();
		setLength(newLength);
		clear(oldLength, newLength);
	}
	
	void increaseOneBlock() {
		
		unsigned int oldLength = getLength();
		unsigned int blockLen = getBlockLen();
		setLength(oldLength+blockLen);
		clear(oldLength, oldLength+blockLen);
	}
	
	void append(ArrayBlock<T> const& append, unsigned int appendLen, unsigned int start, unsigned int startOrigin) {
		
		for (unsigned int i=0; i<appendLen; i++) {
			getElement(startOrigin+i) = append[start+i];
		}
	}
	
	void combine(ArrayBlock<T> const& left, unsigned int leftLen, ArrayBlock<T> const& right, unsigned int rightLen) {
				
		unsigned int i=0;
		for (unsigned int j=0; j<leftLen; j++,i++) {
			getElement(i) = left[j];
		}
		for (unsigned int j=0; j<rightLen; j++,i++) {
			getElement(i) = right[j];
		}
	}
	
	void reverse(ArrayBlock<T>& reverse, unsigned int start, unsigned int end) const {
		
		for (unsigned int i=start; i<end; i++) {
			reverse[end-1-i] = getElement(i);
		}
	}
	
};



#endif /*ARRAYBLOCK_HPP_*/
