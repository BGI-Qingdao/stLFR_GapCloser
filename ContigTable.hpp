/**************************************************
*
* ContigTable.hpp
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

#ifndef CONTIGTABLE_HPP_
#define CONTIGTABLE_HPP_

#include "ContigForFill.hpp"

/* initialization of contigs and gaps from 'fin' */
class ContigTable
{
protected:
	
	ifstream& fin;
	
	ContigForFill* contigs;
	Len_t contigsSum;
	
	Len_t endNumLen;
	Len_t mismatchLen;
	
public:
	
	ContigTable(ifstream& _fin, Len_t _endNumLen=10, Len_t _mismatchLen=5) : 
		fin(_fin), 
		contigs(NULL),
		contigsSum(0), 
		endNumLen(_endNumLen), 
		mismatchLen(_mismatchLen)
	{
		initialize();
	}
	
	virtual ~ContigTable()
	{
		freeMemorys();
	}
	
	ContigForFill* getContigs() const { return contigs; }
	Len_t getContigsSum() const { return contigsSum; }
	
	void initialize() {
		
		cout << ">>>>>>>>>>scaffold initializing<<<<<<<<<<" << endl;
		cout << endl;
		time_t total_start_time=time(NULL);
		
		
		
		cout << "counting scaffolds" << endl;
		time_t start_time=time(NULL);
		
		countContigs(fin);
		
		cout << "scaffolds sum: " << contigsSum <<endl;
		cout << "spent time: " << time(NULL)-start_time << "s" << endl;
		cout << "scaffolds counted" << endl;
		cout << endl;
		
		
		
		cout << "allocating memorys" << endl;
		start_time=time(NULL);
		
		allocateMemorys();
		
		cout << "spent time: " << time(NULL)-start_time << "s" << endl;
		cout << "memorys allocated" << endl;
		cout << endl;
		
		
		
		cout << "initializing scaffolds" << endl;
		start_time=time(NULL);
		
		initializeContigs(fin);
		
		cout << "spent time: " << time(NULL)-start_time << "s" << endl;
		cout << "scaffolds initialized" << endl;
		cout << endl;
		
		
		
		cout << "scaffold sum " << contigsSum <<endl;
		cout << "spent total time: " << time(NULL)-total_start_time << "s" <<endl;
		cout << ">>>>>>>>>>scaffold initialization finished<<<<<<<<<<" << endl;
		cout << endl;
	}
	
	void countContigs(ifstream& fin) {
		
		string str;
		while(!fin.eof()) {
			
			getline(fin, str);
			if (str[0] == '>') {
				
				contigsSum++;
			}
		}
	}
	
	void allocateMemorys() {
				
		contigs = new ContigForFill[contigsSum+1];
	}
	
	void freeMemorys() {
		
		delete [] contigs;
	}
	
	void initializeContig(Len_t iContig, string const& strSeq) {
		
		Len_t contigLength = strSeq.length();
		
		Len_t i=0;
		//ignore 'N' string at head
		for (; i < contigLength; i++) {
			if (strSeq[i]!='N') 
				break;
		}
		
		GapInfo gap;
		Len_t seqLen = 0;
		for (; i < contigLength; i++) {
			
			if (strSeq[i]=='N') {
				
				if (seqLen > 0) {
					TightString* seq = new TightString(&strSeq.c_str()[i-seqLen], seqLen);
					contigs[iContig].appendSequence(seq);
					seqLen = 0;
				}
				
				gap.length++;
			}
			else {
				
				if (gap.length > 0) {
					
					Len_t sameFlag = 1;	

					HashTable< Number_t, Len_t > endNumHash;
					for (Len_t j=0; j<mismatchLen; j++) {
						
//						Len_t currEndNumLen = contigLength - i - j;
//						if (currEndNumLen > endNumLen) 
//							currEndNumLen = endNumLen;
						
						if ((int)(contigLength - i - j) >= (int)endNumLen) {
						
							Number_t num = sequenceToNumber(&strSeq.c_str()[i+j], endNumLen);
							
//							endNum.pos = j;
							
							if (!endNumHash.find(num)) {
								endNumHash[num] = j;

								sameFlag = 0;	
							}
//
							else{
								endNumHash[num] = (Len_t)-1;	
								sameKmer++;
							}
//
						}
					}
					if (contigLength - i < endNumLen){
						sameFlag = 0;
					}
					if (sameFlag == 1){
						sameEnd++;
					}
					
					gap.endNumHash = endNumHash;
					contigs[iContig].appendGap(gap);
					gap.length = 0;
				}
				
				seqLen++;
			}
		}
		
		if (seqLen > 0) {
			TightString* seq = new TightString(&strSeq.c_str()[i-seqLen], seqLen);
			contigs[iContig].appendSequence(seq);
			seqLen = 0;
		}
		Len_t oldMismatchLen = mismatchLen;
		Len_t newMismatchLen = mismatchLen;
		
		Len_t addUnique = 1;
		while (addUnique == 1){
			addUnique = 0;
			oldMismatchLen = newMismatchLen;
			
		LinkedList<GapInfo> &gaps = contigs[iContig].getGaps();
		ListElement<GapInfo> const*gapPtr = gaps.getHead();
		LinkedList<TightString*> const& seqs = contigs[iContig].getSequences();
		ListElement<TightString*> const *seqPtr = seqs.getHead();

		for ( ; gapPtr!=0; gapPtr=gapPtr->getNext()){
			seqPtr = seqPtr->getNext();

			GapInfo const&gapInfo = gapPtr->getDatum();
			TightString* const& tSeq = seqPtr->getDatum();
			Len_t tSeqLen = tSeq->getLength();
			Len_t uniqueNum = 0;
			Len_t i=0;
			for ( ; i<oldMismatchLen; i++){
				if ((int)(tSeqLen-i) >= (int)endNumLen){
					Number_t num;
					tSeq->readTightStringFragment(i, i+endNumLen, num);
					if (gapInfo.endNumHash.find(num) && (gapInfo.endNumHash[num]!=(Len_t)(-1))){
						uniqueNum = 1;
						break;
					}
				}
				else{
					uniqueNum = 1;
				}
			}

			if (uniqueNum == 1)
				continue;
			else{
				Len_t mismatchNum = 0;
				while ((int)(tSeqLen-i) >= (int)endNumLen){
					Number_t num;
					tSeq->readTightStringFragment(i, i+endNumLen, num);
					if (!gapInfo.endNumHash.find(num)){
						gapInfo.endNumHash[num] = i;
						mismatchNum++;

						if (i > newMismatchLen){
								newMismatchLen = i;
						}

						addUnique = 1;

						if (mismatchNum == mismatchLen){
						
							break;
						}
					}
					else{
						gapInfo.endNumHash[num] = (Len_t)-1;
						sameKmer++;
					}

					i++;
				}
				
			}
		}
		}
	}
	
	void initializeContigs(ifstream& fin) {
		
		fin.clear();
		fin.seekg(0, ios_base::beg);
		Len_t iContig=0;
		string str, strName, strSeq;
		while(!fin.eof()) {
			
			getline(fin, str);
			chomp(str);
			if (str[0] == '>') {
				
				if (strSeq.length() > 0) {
					iContig++;
					contigs[iContig].setName(strName);
					initializeContig(iContig, strSeq);
					strSeq = "";
				}
				strName = str;
			}
			else {
				
//				chomp(str);
				strSeq += str;
			}
		}
		if (strSeq.length() > 0) {
			iContig++;
			contigs[iContig].setName(strName);
			initializeContig(iContig, strSeq);
			strSeq = "";
		}
	}
	
	
	
	
};



#endif /*CONTIGTABLE_HPP_*/
