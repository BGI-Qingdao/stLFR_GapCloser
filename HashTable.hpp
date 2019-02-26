/**************************************************
 *
 * HashTable.hpp
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

#ifndef HASHTABLE_HPP_
#define HASHTABLE_HPP_

#include <cmath>
#include <stdexcept>
#include "Common.hpp"

Number_t hash1(char c) {
    return abs (c);
}

Number_t hash1(int i) {
    return abs (i);
}

Number_t hash1(Len_t i) {
    return abs ((int)i);
}

//Number_t hash(Number_t i) {
//	return abs (i);
//}

Number_t hash1(Number_t key){
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

Number_t hash1(double d) {
    if (d == 0)
        return 0;
    else {
        int exponent;
        double mantissa = std::frexp(d, &exponent);
        return (Number_t)(2 * fabs(mantissa) - 1) * ~0LLU;
    }
}

static const Len_t shift = 6;
static const Number_t mask = (Number_t)~0U << (bitsizeof(Number_t) - shift);

//static Number_t hash(char const* str, Len_t len) {
//	
//	Number_t result = 0;
//	for (Len_t i=0; i<len; ++i)
//		result = (result & mask) ^ (result << shift) ^ str[i];
//	return result;
//}

//static Number_t hash(char const* str) {
//	
//	return hash(str, (Len_t)strlen(str));
//}

//static Number_t hash(Kmer const& kmer) {
//	
//	return hash((char const*)kmer.getData(), kmer.getDataNum());
//}

static bool isPrime(Number_t num){
    Number_t i, max;
    if(num < 4) return true;
    if(num % 2 == 0) return false;
    max = (Number_t)sqrt((float)num);
    for(i=3;i<max;i+=2){ if(num % i == 0) return false; }
    return true;
}

static Number_t findNextPrime(Number_t num){
    if(num % 2 == 0) num ++;
    while(true){ if(isPrime(num)) return num; num += 2; }
}



template<class T, class U>
struct HashElement
{
    T key;
    U value;
};

template<class T, class U>
class HashTable
{
    protected:

        static const Short_Len_t EMPTY = 0;
        static const Short_Len_t OCCUPIY = 1;
        static const Short_Len_t DELETE = 2;

    protected:

        Number_t length;  // hash table size
        Number_t count;
        Number_t max;
        float loadFactor;
        Number_t conflictSum;

        T* keys;
        U* values;
        Short_Len_t* flags;

        //the space should be released before distroy a object if 'ownFlag' is true
        bool ownFlag;

    public:

        HashTable(Number_t _length = 3, float _loadFactor = 0.75) {

            initialize(_length, _loadFactor);

        }

        HashTable(HashTable const& hash) {

            length = hash.length;
            count = hash.count;
            max = hash.max;
            loadFactor = hash.loadFactor;
            conflictSum = hash.conflictSum;
            keys = hash.keys;
            values = hash.values;
            flags = hash.flags;
            ownFlag = hash.ownFlag;

            const_cast<HashTable&>(hash).ownFlag = false;
        }

        HashTable& operator= (HashTable const& hash) {

            if (&hash!=this) {
                purge();
                length = hash.length;
                count = hash.count;
                max = hash.max;
                loadFactor = hash.loadFactor;
                conflictSum = hash.conflictSum;
                keys = hash.keys;
                values = hash.values;
                flags = hash.flags;
                ownFlag = hash.ownFlag;

                const_cast<HashTable&>(hash).ownFlag = false;
            }
            return *this;
        }

        virtual ~HashTable() {
            purge();
        }

        void purge() {

            if (ownFlag) {

                delete [] keys;
                delete [] values;
                free(flags);
            }
            keys = NULL;
            values = NULL;
            flags = NULL;
            length = 0;
            count = 0;
            max = 0;
            loadFactor = 0;
            conflictSum = 0;
        }

        void reInitialize(Number_t _length = 3, float _loadFactor = 0.75) {
            purge();
            initialize(_length, _loadFactor);
        }

        void initialize(Number_t _length = 3, float _loadFactor = 0.75) {

            length = findNextPrime(_length);
            count = 0;
            max = (Number_t)(length*_loadFactor);
            loadFactor = _loadFactor;
            conflictSum = 0;
            keys = new T[length];
            values = new U[length];
            flags = (Short_Len_t*)calloc((length+3)/4, 1);
            ownFlag = true;
        }

        Number_t c(Number_t i) const {
            return i;
        }

        Number_t getLength() const { return length; }
        Number_t getCount() const { return count; }
        Number_t getMax() const { return max; }
        float getLoadFactor() const { return loadFactor; }
        Number_t getConflictSum() const { return conflictSum; }

        void increaseLength() {

            Number_t newLength = findNextPrime(length*2);

            T* newKeys = new T[newLength];
            U* newValues = new U[newLength];

            if ((newKeys == NULL) || (newValues == NULL))
                throw std::out_of_range ("out of memory");

            Short_Len_t* newFlags = (Short_Len_t*)calloc((newLength+3)/4, 1);

            T* oldKeys = keys;
            U* oldValues = values;
            Short_Len_t* oldFlags = flags;
            Number_t oldLength = length;
            length = newLength;
            count = 0;
            max = (Number_t)(newLength * loadFactor);
            conflictSum = 0;
            keys = newKeys;
            values = newValues;
            flags = newFlags;

            for (Number_t i=0; i<oldLength; ++i) {
                if (isFlag(oldFlags, EMPTY, i)) continue;
                insert(oldKeys[i], oldValues[i]);
            }

            delete [] oldKeys;
            delete [] oldValues;
            free(oldFlags);

        }

        void insert(T& key, U& value) {
            if (count > max)
                increaseLength();
            //			throw out_of_range ("hash table is full");

            Number_t hashValue = hash1(key);
            for (Number_t i=0; i<count+1; ++i) {
                Number_t probe = (hashValue + c(i)) % length;
                if (!isFlag(OCCUPIY, probe)) {
                    setFlag(OCCUPIY, probe);
                    keys[probe] = key;
                    values[probe] = value;
                    ++count;
                    break;
                }
                else if (key == keys[probe]) {
                    values[probe] = value;
                    break;
                }
                else {
                    conflictSum++;
                }
            }
        }

        Number_t findMatch(T const& key) const {
            Number_t hashValue = hash1(key);
            for (Number_t i=0; i<length; ++i) {
                Number_t probe = (hashValue + c(i)) % length;
                if (isFlag(EMPTY, probe))
                    break;
                if (isFlag(OCCUPIY, probe) && key == keys[probe])
                    return probe;
            }
            return length;
        }

        U* find(T const& key) const {
            Number_t offset = findMatch(key);
            if (offset < length)
                return &values[offset];
            else
                return NULL;
        }

        U& getElement(T& key) const {

            U* pValue = find(key);

            if (pValue)
                return *pValue;
            else {
                U value;
                const_cast<HashTable*>(this)->insert(key, value);
                return *find(key);
            }
        }

        U& operator [](T const& key) const {

            return getElement((T&)key);
        }

        U& operator [](T& key) const {

            return getElement(key);
        }

        //	Number_t findUnoccupied(T& key) const {
        //		Number_t hashValue = hash(key);
        //		for (Number_t i=0; i<count+1; ++i) {
        //			Number_t probe = (hashValue + c(i)) % length;
        //			if (!isFlag(OCCUPIY, probe))
        //				return probe;
        //			if (isFlag(OCCUPIY, probe) && key == array[probe].key)
        //				return probe;
        //		}
        //		return length;
        //	}

        //	Number_t findInstance(T& key) const {
        //		Number_t hashValue = hash(key);
        //		for (Number_t i=0; i<length; ++i) {
        //			Number_t probe = (hashValue + c(i)) % length;
        //			if (isFlag(EMPTY, probe))
        //				break;
        //			if (isFlag(OCCUPIY, probe) && key == array[probe].key)
        //				return probe;
        //		}
        //		return length;
        //	}


        //delete an element
        void withdraw(T& key) {
            if (count == 0)
                throw std::out_of_range ("hash table is empty");
            Number_t i = findMatch(key);
            if (i == length)
                throw std::out_of_range ("object not found");
            for (;;)
            {
                Number_t j;
                for (j=(i+1)%length; isFlag(OCCUPIY, j); j=(j+1)%length)
                {
                    Number_t hashValue = hash1(keys[j]);
                    if ((hashValue <= i && i < j) || (i < j && j < hashValue) || (j < hashValue && hashValue <= i))
                        break;
                }
                if (isFlag(EMPTY, j))
                    break;
                keys[i] = keys[j];
                values[i] = values[j];
                i = j;
            }

            setFlag(EMPTY, i);
            --count;
        }



        void setFlag(Short_Len_t flag, Number_t offset) {

            flags[offset/4] &= ~(3 << offset%4*2);
            flags[offset/4] |= flag << offset%4*2;
        }

        bool isFlag(Short_Len_t flag, Number_t offset) const {

            return ((flags[offset/4] >> offset%4*2) & 3) == flag;
        }

        bool isFlag(Short_Len_t* flags, Short_Len_t flag, Number_t offset) const {

            return ((flags[offset/4] >> offset%4*2) & 3) == flag;
        }



        class Iterator {
            HashTable const& hashTable;
            Number_t index;

            public:

            Iterator(HashTable const& _hashTable) : 
                hashTable(_hashTable) 
            {
                reset();
            }

            bool isDone() const {

                return index >= hashTable.length;
            }

            HashElement<T, U> operator * () const {
                HashElement<T, U> element;
                element.key = hashTable.keys[index];
                element.value = hashTable.values[index];
                return element;
            }

            void operator ++ () {

                //iterate hash table get one element
                do {
                    ++index;
                } while(!isDone() && hashTable.isFlag(EMPTY, index));
            }

            void reset() {

                index = 0;
                while(!isDone() && hashTable.isFlag(EMPTY, index)) {
                    ++index;
                }
            }
        };

        Iterator& newIterator() const {
            return *new Iterator(*this);
        }

        friend class Iterator;

};



#endif /*HASHTABLE_HPP_*/
