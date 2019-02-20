/**************************************************
*
* PairInfo.hpp
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

#ifndef PAIRINFO_HPP_
#define PAIRINFO_HPP_



struct PairInfoElement
{
	Len_t insertSize;
	Number_t beginSerialNum;
	Number_t serialSum;
	Len_t variance;
};

class PairInfo
{
	
protected:
	
	ifstream* finp;
	LibInfo* pLibInfo;
	
	Number_t readsSum;
	
	PairInfoElement* pairInfoArray;
	Len_t arrayLen;
	
public:
	
	PairInfo() : 
		finp(NULL), 
		pLibInfo(NULL), 
		readsSum(0), 
		pairInfoArray(NULL), 
		arrayLen(0)
	{
	}
	
	PairInfo(ifstream* _finp, LibInfo* _pLibInfo=NULL) : 
		finp(_finp), 
		pLibInfo(_pLibInfo), 
		readsSum(0), 
		pairInfoArray(NULL), 
		arrayLen(0)
	{
		initialize();
	}
	
	virtual ~PairInfo()
	{
		if (pairInfoArray)
			delete [] pairInfoArray;
	}
	
	Number_t getReadsSum() const { return readsSum; }
	
	PairInfoElement* getArray() const { return pairInfoArray; }
	Len_t getArrayLen() const { return arrayLen; }
	
	void initialize() {
		
		if (pLibInfo) {
			initializePairInfoByLib(*pLibInfo);
		}
		else {
			initializeReadsSum(*finp);
			initializePairInfo(*finp);
		}
	}
	
	void initializePairInfoByLib(LibInfo& libInfo) {
		
		arrayLen = libInfo.getArrayLen();
		pairInfoArray = new PairInfoElement[arrayLen];
		Number_t beginSerialNum = 0;
		
		for (Len_t i=0; i<arrayLen; i++) {
			
			pairInfoArray[i].insertSize = libInfo.getArray()[i].avg_ins;
			pairInfoArray[i].beginSerialNum = beginSerialNum;
			pairInfoArray[i].serialSum = libInfo.getArray()[i].pairReadsSum;
			int maxVariance = libInfo.getArray()[i].max_ins - libInfo.getArray()[i].avg_ins;
			int minVariance = libInfo.getArray()[i].avg_ins - libInfo.getArray()[i].min_ins;
			pairInfoArray[i].variance = maxVariance>minVariance ? maxVariance:minVariance;
			beginSerialNum += libInfo.getArray()[i].pairReadsSum;
		}
	}
	
	void initializeReadsSum(ifstream& fin) {
		
		const char* delim = "#reads_sum\t";
		size_t found;
		string str, str1;
		
		while(!fin.eof()) {
			getline(fin, str);
			found = str.find(delim);
			if (found!=string::npos) {
				
				str1 = str.substr(found+strlen(delim));
				readsSum = atoll(str1.c_str());
				break;
			}
		}
	}
	
	void initializePairInfo(ifstream& fin) {
		
		fin.clear();
		fin.seekg(0, ios_base::beg);
		
		string str, str1;
		while(!fin.eof()) {
			getline(fin, str);
			if ((str.length() == 0) || (str[0] == '#'))
				continue;
			arrayLen++;
		}
		
		pairInfoArray = new PairInfoElement[arrayLen];
		
		fin.clear();
		fin.seekg(0, ios_base::beg);
		Len_t i = 0;
		while(!fin.eof()) {
			getline(fin, str);
			if ((str.length() == 0) || (str[0] == '#'))
				continue;
			const char* delim = "\t";
			size_t found;
			size_t begin=0;
			found = str.find(delim,begin);
			str1 = str.substr(begin,found-begin);
			pairInfoArray[i].insertSize = atoi(str1.c_str());
			begin = found+1;
			found = str.find(delim,begin);
			str1 = str.substr(begin,found-begin);
			pairInfoArray[i].beginSerialNum = atoi(str1.c_str());
			begin = found+1;
			found = str.find(delim,begin);
			str1 = str.substr(begin,found-begin);
			pairInfoArray[i].serialSum = atoi(str1.c_str());
			begin = found+1;
			str1 = str.substr(begin);
			pairInfoArray[i].variance = atoi(str1.c_str());
			i++;
		}
	}
	
	bool checkInsertSize(Number_t serialNum, Len_t insertSize) const {
		
		Number_t beginSerialNum = 0;
		Number_t serialSum = 0;
		for (Len_t i=0; i<arrayLen; i++) {
			
			if (insertSize == pairInfoArray[i].insertSize) {
				
				beginSerialNum = pairInfoArray[i].beginSerialNum;
				serialSum = pairInfoArray[i].serialSum;
				break;
			}
		}
		
		if (beginSerialNum) {
			if ((serialNum<beginSerialNum) || (serialNum>=(beginSerialNum+serialSum)))
				return false;
			else
				return true;
		}
		else
			return false;
	}
	/*
	bool getPair(Number_t iRead, Number_t& iReadPair, Len_t insertSize=0) const {
		
		Number_t serialNum = iRead >> 1;
		//last bit store the reverse flag
		//because the pair sequences not in same strand.
		bool isReverse = !(iRead & 1);
		
		//check the insertSize
		if (insertSize) {
			
			if (!checkInsertSize(serialNum, insertSize))
				return false;
		}
		
		if (serialNum%2==0) {
			iReadPair = ((serialNum-1)<<1)|isReverse;
		}
		else {
			iReadPair = ((serialNum+1)<<1)|isReverse;
		}
		return true;
	}*/
	
	bool getPair(Number_t iRead, Number_t& iReadPair, Len_t insertSize=0) const {
		
		Number_t serialNum = iRead;
		
		//check the insertSize
		if (insertSize) {
			
			if (!checkInsertSize(serialNum, insertSize))
				return false;
		}
		
		if (serialNum%2==0) {
			iReadPair = serialNum+1;
		}
		else {
			iReadPair = serialNum-1;
		}
		return true;
	}
	
	void getPairs(ReadElement const& readElement, LinkedList<Number_t>& pairs, Len_t insertSize=0) const {
		
		Len_t depth = readElement.getDepth();
		for (Len_t i=0; i<depth; i++) {
			Number_t iRead = readElement.getSerialNums()[i];
			Number_t iReadPair;
			if (getPair(iRead, iReadPair, insertSize))
				pairs.append(iReadPair);
		}
	}
	
	/*
	void getPairs(ReadElement const& readElement, LinkedList<Number_t>* pairs) const {
		
		Len_t depth = readElement.getDepth();
		for (Len_t i=0; i<depth; i++) {
			Number_t iRead = readElement.getSerialNums()[i];
			Number_t serialNum = iRead >> 1;
			//last bit store the reverse flag
			//because the pair sequences not in same strand.
			bool isReverse = !(iRead & 1);
			
			for (Len_t i=0; i<arrayLen; i++) {
				
				Number_t beginSerialNum = pairInfoArray[i].beginSerialNum;
				Number_t endSerialNum = beginSerialNum+pairInfoArray[i].serialSum;
				if ((serialNum>=beginSerialNum) && (serialNum<endSerialNum)) {
					
					Number_t iReadPair;
					if (serialNum%2==0) {
						iReadPair = ((serialNum-1)<<1)|isReverse;
					}
					else {
						iReadPair = ((serialNum+1)<<1)|isReverse;
					}
					pairs[i].append(iReadPair);
					break;
				}
			}
		}
	}*/


	// input: readElemtnt -- currently konwn read
	// output: pairs -- found paired-end reads' sirialNum
	void getPairs(ReadElement const& readElement, LinkedList<Number_t>* pairs) const {
		
		Len_t depth = readElement.getDepth();
		for (Len_t i=0; i<depth; i++) {
			Number_t iRead = readElement.getSerialNums()[i];
			Number_t serialNum = iRead;
			
			for (Len_t i=0; i<arrayLen; i++) {
				
				Number_t beginSerialNum = pairInfoArray[i].beginSerialNum;
				Number_t endSerialNum = beginSerialNum+pairInfoArray[i].serialSum;
				if ((serialNum>=beginSerialNum) && (serialNum<endSerialNum)) {
					
					Number_t iReadPair;
					if (serialNum%2==0) {
						iReadPair = serialNum+1;
					}
					else {
						iReadPair = serialNum-1;
					}
					pairs[i].append(iReadPair);
					break;
				}
			}
		}
	}
	
	Len_t getMaxInsertSize() const {
		
		if (arrayLen > 0)
			return pairInfoArray[arrayLen-1].insertSize;
		else 
			return 200;
	}
	
	Len_t getMinInsertSize() const {
		
		if (arrayLen > 0)
			return pairInfoArray[0].insertSize;
		else 
			return 200;
	}
	
};

#endif /*PAIRINFO_HPP_*/
