/**************************************************
 *
 * Utils.cpp
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

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <stdexcept>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <set>
#include <algorithm>

#define bitsizeof(T) (sizeof(T)* 8 )

inline void chomp(char* str) {

    size_t size = strlen(str);
    while( size>0 && ((str[size-1] == '\r') || (str[size-1] == '\n')) ) {
        str[size-1] = '\0';
        --size;
    }
}

inline void chomp(std::string& str) {

    int i=str.length()-1;
    while ((i>=0) && ((str[i] == '\r') || (str[i] == '\n')))
        i--;
    if (i < 0) {
        str = "";
    } else {
        str = str.substr(0, i+1);
    }
}

    template< class T>
void SetAdd( const std::set<T>  s1 ,
        const std::set<int> & s2 )
{
    for( const auto & i : s2 )
        s1.insert(i);
}
    template< class T>
std::set<T > SetUnion( const std::set<T> & s1 ,
        const std::set<int> & s2 )
{
    std::set<T> dest1;
    std::set_union(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

    template< class T>
std::set<T > SetDiff( const std::set<T > & s1 ,
        const std::set<T> & s2 )
{
    std::set<T> dest1;
    std::set_difference(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

    template< class T>
std::set<T> SetIntersection( const std::set<T> & s1 ,
        const std::set<T> & s2 )
{
    std::set<T> dest1;
    std::set_intersection(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

#endif /*UTILS_HPP_*/
