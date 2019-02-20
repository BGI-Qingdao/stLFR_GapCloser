/**************************************************
*
* StraightInsertionSorter.hpp
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

#ifndef STRAIGHTINSERTIONSORTER_HPP_
#define STRAIGHTINSERTIONSORTER_HPP_

#include "sort/ClassesOfSorters.hpp"

namespace data_structure
{



template<class T> class StraightInsertionSorter : public InsertionSorter<T> {
protected:
	void doSort(Array<T>& array) {
		doSort(array.getData());
	}
	
	void doSort(T array[]) {
		for (unsigned int i=1; i < this->n; ++i)
			for (unsigned int j=i; j>0 && array[j-1U]>array[j]; --j)
				swap(array[j], array[j-1U]);
	}
};



}

#endif /*STRAIGHTINSERTIONSORTER_HPP_*/
