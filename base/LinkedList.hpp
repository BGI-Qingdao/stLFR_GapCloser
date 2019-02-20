/**************************************************
*
* LinkedList.hpp
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

#ifndef LINKEDLIST_HPP_
#define LINKEDLIST_HPP_

#include <stdexcept>
using namespace std;

namespace data_structure {

template<class T>
class LinkedList;

template<class T>
class ListElement {
	T datum;
	ListElement* next;

	ListElement(T const& _datum, ListElement<T>* _next) :
		datum(_datum), next(_next) {
	}
public:

	T const& getDatum() const {
		return datum;
	}

	ListElement<T> const* getNext() const {
		return next;
	}

	friend class LinkedList<T>;
};

/*
define operations on linked list, including:
  1) get first/last element.
  2) search element.
  3) search and delete element.
  4) insert new element from head/tail or before/after specific element.
*/
template<class T>
class LinkedList {
	ListElement<T>* head;
	ListElement<T>* tail;
	unsigned int count;
public:

	LinkedList() :
		head(0), tail(0), count(0) {
	}

	~LinkedList() {
		purge();
	}

	LinkedList(LinkedList<T> const& linkedList) :
		head(0), tail(0), count(0) {
		ListElement<T> const* ptr;
		for (ptr = linkedList.head; ptr != 0; ptr = ptr->next)
			append(ptr->datum);
	}

	LinkedList<T>& operator =(LinkedList<T> const& linkedList) {
		if (&linkedList != this) {
			purge();
			ListElement<T> const* ptr;
			for (ptr = linkedList.head; ptr != 0; ptr = ptr->next)
				append(ptr->datum);
		}
		return *this;
	}

	void purge() {
		while (head != 0) {
			ListElement<T>* const tmp = head;
			head = head->next;
			delete tmp;
		}
		tail = 0;
		count = 0;
	}

	unsigned int getCount() const {
		return count;
	}
	
	ListElement<T> const* getHead() const {
		return head;
	}

	
	ListElement<T> const* getTail() const {
		return tail;
	}

	bool isEmpty() const {
		return head == 0;
	}

	T const& first() const {
		if (head == 0)
			throw domain_error ("list is empty");
		return head->datum;
	}

	T const& last() const {
		if (tail == 0)
			throw domain_error ("list is empty");
		return tail->datum;
	}

	void prepend(T const& item) {
		ListElement<T>* const tmp = new ListElement<T> (item, head);
		if (head == 0)
			tail = tmp;
		head = tmp;
		++count;
	}

	void append(T const& item) {
		ListElement<T>* const tmp = new ListElement<T> (item, 0);
		if (head == 0)
			head = tmp;
		else
			tail->next = tmp;
		tail = tmp;
		++count;
	}

	void extract(T const& item) {
		ListElement<T>* ptr = head;
		ListElement<T>* prevPtr = 0;
		while (ptr != 0&& ptr->datum != item) {
			prevPtr = ptr;
			ptr = ptr->next;
		}
		if (ptr == 0)
			throw invalid_argument ("item not found");
		if (ptr == head)
			head = ptr->next;
		else
			prevPtr->next = ptr->next;
		if (ptr == tail)
			tail = prevPtr;
		delete ptr;
		--count;
	}

	void myExtract(T const& item) {
		ListElement<T>* ptr = head;
		ListElement<T>* prevPtr = 0;
		while (ptr != 0&& ptr->datum != item) {
			prevPtr = ptr;
			ptr = ptr->next;
		}
		if (ptr == 0)
			return;
		if (ptr == head)
			head = ptr->next;
		else
			prevPtr->next = ptr->next;
		if (ptr == tail)
			tail = prevPtr;
		delete ptr;
		--count;
	}

	bool find(T const& item) {
		ListElement<T>* ptr = head;
		ListElement<T>* prevPtr = 0;
		while (ptr != 0&& ptr->datum != item) {
			prevPtr = ptr;
			ptr = ptr->next;
		}
		if (ptr == 0)
			return false;
		else
			return true;
		
	}

	void insertAfter(ListElement<T> const* arg, T const& item) {
		ListElement<T>* ptr = const_cast<ListElement<T>*> (arg);
		if (ptr == 0)
			throw invalid_argument ("invalid position");
		ListElement<T>* const tmp =new ListElement<T> (item, ptr->next);
		ptr->next = tmp;
		if (tail == ptr)
			tail = tmp;
		++count;
	}

	void insertBefore(ListElement<T> const* arg, T const& item) {
		ListElement<T>* ptr = const_cast<ListElement<T>*> (arg);
		if (ptr == 0)
			throw invalid_argument ("invalid position");
		ListElement<T>* const tmp = new ListElement<T> (item, ptr);
		if (head == ptr)
			head = tmp;
		else {
			ListElement<T>* prevPtr = head;
			while (prevPtr != 0&& prevPtr->next != ptr)
				prevPtr = prevPtr->next;
			if (prevPtr == 0)
				throw invalid_argument ("invalid position");
			prevPtr->next = tmp;
		}
		++count;
	}
};

}

#endif /*LINKEDLIST_HPP_*/
