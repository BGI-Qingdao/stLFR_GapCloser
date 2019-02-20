/**************************************************
*
* QuickSorter.hpp
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

#ifndef QUICKSORTER_HPP_
#define QUICKSORTER_HPP_

#include "sort/ClassesOfSorters.hpp"
#include "sort/StraightInsertionSorter.hpp"

namespace data_structure
{



template<class T> class QuickSorter : public ExchangeSorter<T> {
protected:
	static unsigned int const cutOff = 2;

//	virtual unsigned int selectPivot(Array<T>&, unsigned int, unsigned int) = 0;
	virtual unsigned int selectPivot(T[], unsigned int, unsigned int) = 0;
	
	void doSort(Array<T>& array, unsigned int left, unsigned int right) {
		
		doSort(array.getData(), left, right);
	}
	
	void doSort(T array[], unsigned int left, unsigned int right) {
		
		if (right-left+1 > cutOff) {
			
			unsigned int const p = selectPivot(array, left, right);
			swap(array[p], array[right]);
			T& pivot = array[right];
			unsigned int i = left;
			unsigned int j = right - 1U;
			for (;;)
			{
				while (i < j && array[i] < pivot) ++i;
				while (i < j && array[j] > pivot) --j;
				if (i >= j) break;
				swap(array[i++], array[j--]);
			}
			if (array[i] > pivot)
				swap(array[i], pivot);
			if (left < i)
				doSort(array, left, i-1U);
			if (right > i)
				doSort(array, i+1, right);
		}
	}

	void doSort(Array<T>& array) {
		
		doSort(array.getData());
	}
	
	void doSort(T array[]){
		
		doSort(array, 0, this->n - 1U);
		StraightInsertionSorter<T> s;
		s.sort(array, this->n);
	}

};



}

#endif /*QUICKSORTER_HPP_*/
