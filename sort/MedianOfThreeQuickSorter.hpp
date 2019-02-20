/**************************************************
*
* MedianOfThreeQuickSorter.hpp
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

#ifndef MEDIANOFTHREEQUICKSORTER_HPP_
#define MEDIANOFTHREEQUICKSORTER_HPP_

#include "sort/QuickSorter.hpp"

namespace data_structure
{



template<class T> class MedianOfThreeQuickSorter : public QuickSorter<T> {
protected:
//	unsigned int selectPivot(Array<T>& array, unsigned int left, unsigned int right) {
//		
//		return selectPivot(array.getData(), left, right);
//	}
	
	static void swap(unsigned int& x, unsigned int& y) {
		unsigned int const tmp = x;
		x = y;
		y = tmp;
	}
	
	unsigned int selectPivot(T array[], unsigned int left, unsigned int right) {
		unsigned int middle = (left+right)/2;
		if (array[left] > array[middle])
			swap(left, middle);
		if (array[left] > array[right])
			swap(left, right);
		if (array[middle] > array[right])
			swap(middle, right);
		return middle;
	}
};



}

#endif /*MEDIANOFTHREEQUICKSORTER_HPP_*/
