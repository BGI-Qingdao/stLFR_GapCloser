/**************************************************
*
* Array.hpp
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

#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include <stdexcept>
using namespace std;

namespace data_structure {


/*
template<class T>
class Array {
protected:
	T* data;
	unsigned int length;
public:

	Array& operator =(Array const&);

	Array() :
		data(new T [0]), length(0) {
	}

	Array(unsigned int n) :
		data(new T [n]), length(n) {
	}

	Array(Array<T> const& array) :
		data(new T [array.length]), length(array.length) {
		for (unsigned int i = 0; i < length; ++i)
			data[i] = array.data[i];
	}

	~Array() {
		delete [] data;
	}

	T const& operator [](unsigned int position) const {
		if (position >= length)
			throw out_of_range ("invalid position");
		return data[position];
	}

	T& operator [](unsigned int position) {
		if (position >= length)
			throw out_of_range ("invalid position");
		return data[position];
	}

	T const* getData() const {
		return data;
	}

	T* getData() {
		return data;
	}

	unsigned int getLength() const {
		return length;
	}

	void setLength(unsigned int newLength) {
		T* const newData = new T [newLength];
		unsigned int const min =length < newLength ? length : newLength;
		for (unsigned int i = 0; i < min; ++i)
			newData [i] = data[i];
		delete [] data;
		data = newData;
		length = newLength;
	}
};
*/



template<class T>
class Array {
protected:
	T* data;
	unsigned int base;
	unsigned int length;
public:

	Array& operator =(Array const&);

	Array() :
		data(new T [0]), base(0), length(0) {
	}

	Array(unsigned int n, unsigned int m = 0) :
		data(new T [n]), base(m), length(n) {
	}

	Array(Array<T> const& array) :
		data(new T [array.length]), base(array.base), length(array.length) {
		for (unsigned int i = 0; i < length; ++i)
			data[i] = array.data[i];
	}

	~Array() {
		delete [] data;
	}

	T const& operator [](unsigned int position) const {
		unsigned int const offset = position - base;
		if (offset >= length)
			throw out_of_range ("invalid position");
		return data[offset];
	}

	T& operator [](unsigned int position) {
		unsigned int const offset = position - base;
		if (offset >= length)
			throw out_of_range ("invalid position");
		return data[offset];
	}

	T const* getData() const {
		return data;
	}
	
	T* getData() {
		return data;
	}

	unsigned int getBase() const {
		return base;
	}

	unsigned int getLength() const {
		return length;
	}

	void setBase(unsigned int newBase) {
		base = newBase;
	}

	void setLength(unsigned int newLength) {
		T* const newData = new T [newLength];
		unsigned int const min =length < newLength ? length : newLength;
		for (unsigned int i = 0; i < min; ++i)
			newData [i] = data[i];
		delete [] data;
		data = newData;
		length = newLength;
	}
};



}

#endif /*ARRAY_HPP_*/
