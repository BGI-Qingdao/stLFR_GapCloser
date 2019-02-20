/**************************************************
*
* ClassesOfSorters.hpp
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

#ifndef CLASSESOFSORTERS_HPP_
#define CLASSESOFSORTERS_HPP_

#include "sort/Sorter.hpp"

namespace data_structure
{



template<class T> class InsertionSorter : public Sorter<T> {
};

template<class T> class ExchangeSorter : public Sorter<T> {
};

template<class T> class SelectionSorter : public Sorter<T> {
};

template<class T> class MergeSorter : public Sorter<T> {
};

template<class T> class DistributionSorter : public Sorter<T> {
};



}

#endif /*CLASSESOFSORTERS_HPP_*/
