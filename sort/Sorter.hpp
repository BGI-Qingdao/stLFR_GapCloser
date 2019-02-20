/**************************************************
*
* Sorter.hpp
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

#ifndef SORTER_HPP_
#define SORTER_HPP_

#include "base/Array.hpp"

namespace data_structure
{



template<class T> class Sorter {
protected:
	unsigned int n;

	static void swap(T& x, T& y) {
		T const tmp = x;
		x = y;
		y = tmp;
	}

	virtual void doSort(Array<T>&) = 0;
	virtual void doSort(T[]) = 0;
public:
	void sort(Array<T>& array) {
		n = array.getLength();
		if (n > 0) {
			unsigned int const tmp = array.getBase();
			array.setBase(0);
			doSort(array);
			array.setBase(tmp);
		}
	}
	
	void sort(T array[], unsigned int len) {
		n = len;
		if (n > 0) {
			doSort(array);
		}
	}
};



}

#endif /*SORTER_HPP_*/
