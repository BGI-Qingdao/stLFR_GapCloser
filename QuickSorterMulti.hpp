/**************************************************
*
* QuickSorterMulti.hpp
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

#ifndef QUICKSORTERMULTI_HPP_
#define QUICKSORTERMULTI_HPP_



template<class T>
class QuickSorterSingle
{
protected:
	
	typedef int (fun_compare_t)(T const& leftKey, T const& rightKey);
	
public:
	
	void doSort(T keys[], fun_compare_t& compare, unsigned int length) {
		
		doQuickSort(keys, compare, 0, length-1);
		doStraightInsertionSort(keys, compare, length);
	}
	
protected:
	
	void doStraightInsertionSort(T keys[], fun_compare_t& compare, unsigned int length) {
		for (unsigned int i=1; i <length; ++i)
			for (unsigned int j=i; j>0 && compare(keys[j-1U],keys[j])>0; --j) {
				swap(keys[j], keys[j-1U]);
			}
	}
	
	void doQuickSort(T keys[], fun_compare_t& compare, unsigned int left, unsigned int right) {
		
		unsigned int const cutOff = 2;
		if (right-left+1 > cutOff) {
			
			unsigned int const p = selectPivot(keys, compare, left, right);
			swap(keys[p], keys[right]);
			
			T& pivot = keys[right];
			unsigned int i = left;
			unsigned int j = right - 1U;
			for (;;)
			{
				while (i < j && compare(keys[i],pivot)<0) ++i;
				while (i < j && compare(keys[j],pivot)>0) --j;
				if (i >= j) break;
				swap(keys[i], keys[j]);
				i++;
				j--;
			}
			if (compare(keys[i],pivot)>0) {
				swap(keys[i], pivot);
			}
			if (left < i)
				doQuickSort(keys, compare, left, i-1U);
			if (right > i)
				doQuickSort(keys, compare, i+1, right);
		}
	}
	
	template<class V>
	void swap(V& x, V& y) {
		V const tmp = x;
		x = y;
		y = tmp;
	}
	
	unsigned int selectPivot(T array[], fun_compare_t& compare, unsigned int left, unsigned int right) {
		unsigned int middle = (left+right)/2;
		if (compare(array[left],array[middle])>0)
			swap(left, middle);
		if (compare(array[left],array[right])>0)
			swap(left, right);
		if (compare(array[middle],array[right])>0)
			swap(middle, right);
		return middle;
	}
};



template<class T, class U>
class QuickSorterMulti
{
protected:
	
	typedef int (fun_compare_t)(T const& leftKey, T const& rightKey);
	
public:
	
	void doSort(T keys[], U values[], fun_compare_t& compare, unsigned int length) {
		
		doQuickSort(keys, values, compare, 0, length-1);
		doStraightInsertionSort(keys, values, compare, length);
	}
	
protected:
	
	void doStraightInsertionSort(T keys[], U values[], fun_compare_t& compare, unsigned int length) {
		for (unsigned int i=1; i <length; ++i)
			for (unsigned int j=i; j>0 && compare(keys[j-1U],keys[j])>0; --j) {
				swap(keys[j], keys[j-1U]);
				swap(values[j], values[j-1U]);
			}
	}
	
	void doQuickSort(T keys[], U values[], fun_compare_t& compare, unsigned int left, unsigned int right) {
		
		unsigned int const cutOff = 2;
		if (right-left+1 > cutOff) {
			
			unsigned int const p = selectPivot(keys, compare, left, right);
			swap(keys[p], keys[right]);
			swap(values[p], values[right]);
			
			T& pivot = keys[right];
			unsigned int i = left;
			unsigned int j = right - 1U;
			for (;;)
			{
				while (i < j && compare(keys[i],pivot)<0) ++i;
				while (i < j && compare(keys[j],pivot)>0) --j;
				if (i >= j) break;
				swap(keys[i], keys[j]);
				swap(values[i], values[j]);
				i++;
				j--;
			}
			if (compare(keys[i],pivot)>0) {
				swap(keys[i], pivot);
				swap(values[i], values[right]);
			}
			if (left < i)
				doQuickSort(keys, values, compare, left, i-1U);
			if (right > i)
				doQuickSort(keys, values, compare, i+1, right);
		}
	}
	
	template<class V>
	void swap(V& x, V& y) {
		V const tmp = x;
		x = y;
		y = tmp;
	}
	
	unsigned int selectPivot(T array[], fun_compare_t& compare, unsigned int left, unsigned int right) {
		unsigned int middle = (left+right)/2;
		if (compare(array[left],array[middle])>0)
			swap(left, middle);
		if (compare(array[left],array[right])>0)
			swap(left, right);
		if (compare(array[middle],array[right])>0)
			swap(middle, right);
		return middle;
	}
	
};





#endif /*QUICKSORTERMULTI_HPP_*/
