/**************************************************
*
* ContigAssembler.hpp
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

#ifndef CONTIGASSEMBLER_HPP_
#define CONTIGASSEMBLER_HPP_

#include "Contig.hpp"
//#include <map>

class ContigAssembler
{
public:
	static const Short_Len_t fixedOverlapMode = 1;
	static const Len_t maxReadsCount = 10000;
	
protected:
	char* outfile;
	ofstream& fout;
//	ofstream foutDG;
	ofstream foutContig;
	ReadAccessor& readAccessor;
	PairInfo const& pairInfo;
	
	Len_t threadSum;
	pthread_t* threadsID;
	
	Len_t numberOfContigs;
	Len_t maxReadLength;
	Len_t reExtendLength;
	
	Short_Len_t overlapMode;
	Short_Len_t overlapParam;
	
	pthread_mutex_t mutexNumberOfContigs;
	pthread_mutex_t mutexOutput;
	pthread_mutex_t mutexContigsPos;
	
public:
	
	ContigAssembler(char* _outfile, ofstream& _fout, ReadAccessor& _readAccessor, PairInfo const& _pairInfo, Len_t _threadSum, Len_t _maxReadLength=35, Short_Len_t _overlapMode=fixedOverlapMode, Short_Len_t _overlapParam=25) : 
		outfile(_outfile), 
		fout(_fout), 
		readAccessor(_readAccessor), 
		pairInfo(_pairInfo), 
		threadSum(_threadSum), 
		threadsID(new pthread_t[_threadSum]), 
		numberOfContigs(0), 
		maxReadLength(_maxReadLength), 
		reExtendLength(_maxReadLength), 
		overlapMode(_overlapMode), 
		overlapParam(_overlapParam)
	{
		pthread_mutex_init(&mutexNumberOfContigs, NULL);
		pthread_mutex_init(&mutexOutput, NULL);
		pthread_mutex_init(&mutexContigsPos, NULL);
//		assemble();
	}
	
	virtual ~ContigAssembler()
	{
		delete [] threadsID;
		pthread_mutex_destroy(&mutexNumberOfContigs);
		pthread_mutex_destroy(&mutexOutput);
		pthread_mutex_destroy(&mutexContigsPos);
	}
	
public:
	
	void assemble() {
		
		cout << ">>>>>>>>>>assembling<<<<<<<<<<" << endl;
		cout << endl;
		time_t total_start_time=time(NULL);
		
//		char dgName[MAX_STRING_LEN];
//		strcpy(dgName, outfile);
//		strcat(dgName,".dg");
//		foutDG.open(dgName);
		
		char contigName[MAX_STRING_LEN];
		strcpy(contigName, outfile);
		strcat(contigName,".contig");
		foutContig.open(contigName);
		
		Contig::readAccessor = &readAccessor;
		
		createThread();
		
//		foutDG.close();
		foutContig.close();
		
		cout << "spend total time " << time(NULL)-total_start_time<<endl;
		cout << ">>>>>>>>>>assembling finished<<<<<<<<<<" << endl;
		cout << endl;
	}
	
	
	
protected:
	
	struct Arg 
	{
		ContigAssembler* pContigAssembler;
		Len_t iThread;
	};
	
	void createThread() {
		
		for (Len_t i=0; i<threadSum; i++) {
			Arg* pArg = new Arg;
			pArg->pContigAssembler = this;
			pArg->iThread = i;
			pthread_create(&threadsID[i], NULL, (void *(*)(void *))threadFunc, (void*)pArg);
		}
		
		for (Len_t i=0; i<threadSum; i++) {
			pthread_join(threadsID[i], NULL);
		}
	}
	
	static void threadFunc(Arg* pArg) {
		
		ContigAssembler& contigAssembler = *(pArg->pContigAssembler);
		Len_t iThread = pArg->iThread;
		delete pArg;
		
		contigAssembler.assembleInThread(iThread); 
	}
	
	
	
	void assembleInThread(Len_t iThread) { //assembly by contig.
		
		pthread_mutex_lock(&mutexNumberOfContigs);
		numberOfContigs++;
		Len_t iContig = numberOfContigs;
		cout << "constructing " << numberOfContigs << " contig" <<endl;
		pthread_mutex_unlock(&mutexNumberOfContigs);
		
		ReadAccessor::Iterator& pReadElement = readAccessor.newIterator(readAccessor.getReadsSum()*iThread/threadSum);
		while (!pReadElement.isDone()) {
			
			ReadElement readElement = *pReadElement;
			Len_t readLength = readElement.getLen();
			
			
			
			////////////////////////////////unique contig assembly////////////////////////////////
			
			Contig uniqueContigForward(readLength);
			
			readElement.getSequence(uniqueContigForward.getTightString());
			
			Len_t depth = readElement.getDepth();
			for(Len_t i=0; i<readLength; i++) {
				uniqueContigForward.getDepths()[i] += depth;
			}
			
			//forward unique contig assembly
			assembleUniqueContig(uniqueContigForward, maxReadLength*2);
			Len_t uniqueContigForwardLen = uniqueContigForward.getLength();
			
			if (uniqueContigForwardLen < maxReadLength*2) {
				++pReadElement;
				continue;
			}
			
			
			
			Contig uniqueContigReverse(readLength);
			
			readElement.getSequenceReverse(uniqueContigReverse.getTightString());
			
			for(Len_t i=0; i<readLength; i++) {
				uniqueContigReverse.getDepths()[i] += depth;
			}
			
			//reverse unique contig assembly
			assembleUniqueContig(uniqueContigReverse, maxReadLength*2);
			Len_t uniqueContigReverseLen = uniqueContigReverse.getLength();
			
			if (uniqueContigReverseLen < maxReadLength*2) {
				++pReadElement;
				continue;
			}
			
			
			
			////////////////////////////////contig assembly by pair-end////////////////////////////////			
			Contig contigReverse(readLength, iContig);
			
			readElement.getSequenceReverse(contigReverse.getTightString());
			
			for(Len_t i=0; i<readLength; i++) {
				contigReverse.getDepths()[i] += depth;
			}
			
			//reverse contig assembly by pair-end, strand flag set to 1
			Number_t endOtherContigPosReverse = 0;
			assembleContigByPE(contigReverse, endOtherContigPosReverse, 1);
			Len_t contigReverseLen = contigReverse.getLength();
			
			if (contigReverseLen < overlapParam) {
				++pReadElement;
				continue;
			}
			
			
			
			Contig contigForward(0, iContig);
			contigReverse.reverse(contigForward, 0, contigReverseLen);
			
			//forward contig assembly by pair-end, strand flag set to 0
			Number_t endOtherContigPosForward = 0;
			assembleContigByPE(contigForward, endOtherContigPosForward, 0);
			Len_t contigForwardLen = contigForward.getLength();
			
			if (contigForwardLen < maxReadLength) {
				
				++pReadElement;
				continue;
			}
			
			
//			contigForward.setEndOtherContigPosReverse(endOtherContigPosReverse);
//			contigForward.setEndOtherContigPosForward(endOtherContigPosForward);
			
			
			pthread_mutex_lock(&mutexOutput);
			contigForward.outputContigFastAFormat(fout);
			contigForward.outputContig(foutContig);
			pthread_mutex_unlock(&mutexOutput);
			
			++pReadElement;
			if (!pReadElement.isDone()) {
				pthread_mutex_lock(&mutexNumberOfContigs);
				numberOfContigs++;
				iContig = numberOfContigs;
				cout << "constructing " << numberOfContigs << " contig" <<endl;
				pthread_mutex_unlock(&mutexNumberOfContigs);
			}
		}
		delete &pReadElement;
	}
	
	void assembleContigByPE(Contig& contig, Number_t& endOtherContigPos, Short_Len_t forwardReverseFlag, Len_t maxContigLength=0) {
		
		switch(overlapMode) {
		
			case fixedOverlapMode:
				fixedOverlapByPE(contig, endOtherContigPos, forwardReverseFlag, maxContigLength);
				break;
			default:
				fixedOverlapByPE(contig, endOtherContigPos, forwardReverseFlag, maxContigLength);
				break;
		}
	}
	
	void getContigPos(LinkedList<Number_t>& ids, LinkedList<Len_t>& contigPos, Contig& contig) {
		
		
		ListElement<Number_t> const* ptrId;
		for (ptrId = ids.getHead(); ptrId != 0; ptrId = ptrId->getNext()) {
			
			Number_t id = ptrId->getDatum();
			
			LinkedList<Len_t> pos;
			contig.getContigPos(id, pos);
			
			ListElement<Len_t> const* ptr;
			for (ptr = pos.getHead(); ptr != 0; ptr = ptr->getNext()) {
				contigPos.append(ptr->getDatum());
			}
		}
	}
	
	void getContigPosByPair(ReadElement const& readElement, LinkedList<Len_t>* contigPos, Contig& contig) {
		
		Len_t arrayLen = pairInfo.getArrayLen();
		LinkedList<Number_t>* ids = new LinkedList<Number_t>[arrayLen];
		readAccessor.getReadIdsByPair(readElement, ids);
		
		for (Len_t i=0; i<arrayLen; i++) {
			
			getContigPos(ids[i], contigPos[i], contig);
		}
		
		delete [] ids;
	}
	
//	void getContigPos(LinkedList<ReadElement>& reads, LinkedList<Len_t>& contigPos, Contig& contig) {
//		
//		ListElement<ReadElement> const* ptrRead;
//		for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
//			
//			ReadElement const& readElement = ptrRead->getDatum();
//			
//			LinkedList<Len_t> pos;
//			contig.getContigPos(readElement, pos);
//			
//			ListElement<Len_t> const* ptr;
//			for (ptr = pos.getHead(); ptr != 0; ptr = ptr->getNext()) {
//				contigPos.append(ptr->getDatum());
//			}
//		}
//	}
//	
//	void getContigPosByPair(ReadElement const& readElement, LinkedList<Len_t>& contigPos, Contig& contig, Len_t insertSize=0) {
//		
//		LinkedList<ReadElement> reads;
//		readAccessor.getReadsByPair(readElement, reads, insertSize);
//		
//		getContigPos(reads, contigPos, contig);
//	}
//	
//	void getContigPosByPair(ReadElement const& readElement, LinkedList<Len_t>* contigPos, Contig& contig) {
//		
//		Len_t arrayLen = pairInfo.getArrayLen();
//		LinkedList<ReadElement>* reads = new LinkedList<ReadElement>[arrayLen];
//		readAccessor.getReadsByPair(readElement, reads);
//		
//		for (Len_t i=0; i<arrayLen; i++) {
//			
//			getContigPos(reads[i], contigPos[i], contig);
//		}
//		
//		delete [] reads;
//	}
	
	struct Repeat 
	{
		Len_t pos;
		Len_t len;
	};
	
	bool selfCircleCheckByPE(
			LinkedList<ReadElement> const& reads, 
			LinkedList<Repeat>*& pOldRepeats, 
			Contig& contig
			) {
		
		TightString& tStrContig = contig.getTightString();
		ArrayBlock< LinkedList<Number_t> >& readPositions = contig.getReadPositions();
		
		LinkedList<Repeat>& oldRepeats = *pOldRepeats;
		LinkedList<Repeat>& newRepeats = *new LinkedList<Repeat>;
		
		LinkedList<Len_t> pos;
		ListElement<ReadElement> const* ptrRead;
		for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
			
			ReadElement const& readElement = ptrRead->getDatum();
			contig.getContigPos(readElement, pos);
		}
		
		
		
		bool result = false;
		if (pos.getCount()>0) {
			
			//compare repeat area, if repeat length great than max insert size,
			//then we find a self circle
			if (oldRepeats.getCount()>0) {
				
				ListElement<Repeat> const* ptr;
				for (ptr = oldRepeats.getHead(); ptr != 0; ptr = ptr->getNext()) {
					
					Repeat& repeat = const_cast<Repeat&>(ptr->getDatum());
					
					//get next repeat position in this repeat area
					Len_t repeatPosNext=repeat.pos+repeat.len;
					for (; repeatPosNext<tStrContig.getLength(); repeatPosNext++) {
						
						if (readPositions[repeatPosNext].getCount())
							break;
					}
					
					if (repeatPosNext<tStrContig.getLength()) {
						
						ListElement<Number_t> const* ptr;
						for (ptr = readPositions[repeatPosNext].getHead(); ptr != 0; ptr = ptr->getNext()) {
							
							for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
								
								ReadElement const& readElement = ptrRead->getDatum();
								if (readElement.getID() == ptr->getDatum()) {
									
									repeat.len = repeatPosNext-repeat.pos+1;
									if (repeat.len > pairInfo.getMaxInsertSize()) {
										delete pOldRepeats;
										pOldRepeats = &newRepeats;
										return true;
									}
									else {
										newRepeats.append(repeat);
										goto break1;
									}
								}
							
							}
						}
						break1:;
					}
				}
			}
			//first meet repeat area
			else {
				
				ListElement<Len_t> const* ptr;
				for (ptr = pos.getHead(); ptr != 0; ptr = ptr->getNext()) {
					
					Repeat repeat;
					repeat.pos = ptr->getDatum();
					repeat.len = 1;
					newRepeats.append(repeat);
				}
			}
		}
		
		delete pOldRepeats;
		pOldRepeats = &newRepeats;
		return result;
	}

	// function:     matrixProcess
	// description: deal with sequence matrix consisting of contig and reads
	// input:          strMatrix -- sequence matrix consisting of contig and reads
	//                   depthsMatrix -- depth matrix consisting of contig and reads' depth
	//                   xLength -- maximum consensus length
	//                   yLength -- reads number plusing contig number
	//                   cutoffScore -- cutoff for dominating base
	//                   errorCutoff -- maximum number of mismatch
	// output:        tStrFinal -- consensus sequence
	//                   depthsFinal -- depth of consensus
	//                   usedFlag -- flag to indicate whether a read is used or not, 1 for used
	//                   replaceFlag -- flag to indicate whether to replace contig with consensus, 1 to replace
	void matrixProcess(
			char** strMatrix, 
			Len_t** depthsMatrix, 
			TightString& tStrFinal, 
			Len_t* depthsFinal, 
			LinkedList<Len_t> &usedFlag,		
			int& branchPosition, 
			Len_t &xLength, 	
			Len_t yLength, 
			int &replaceFlag,		
			float cutoffScore=1.0, 
			Len_t errorCutoff=2
			) {
		
		Len_t errorSum = 0;
		int branchBegin = -1;
		Len_t i,j;

		Len_t leftReadNum = yLength - 1;


		Len_t *errorArr = new Len_t[yLength];
		for (Len_t m=0; m<yLength; m++){
			errorArr[m] = 0;
		}
		
		for(j=0; j<xLength; ) {
			//statistic
			Len_t bitScore[4]={0,0,0,0};
			for(i=0;i<yLength;i++) {
				
				if (strlen(strMatrix[i])>j)
					bitScore[nucleotideToNumber(strMatrix[i][j])] += depthsMatrix[i][j];
			}
			Len_t scoreSum=0;
			for(i=0;i<4;i++) {
				scoreSum+=bitScore[i];
			}
			//get nucleotide which score above the cutoff score
			int nucleotide = -1;
			for(i=0;i<4;i++) {
				if (bitScore[i]/(float)scoreSum >= cutoffScore) {
					nucleotide = i;
					break;
				}
			}
			
			//count error sum
			if ((nucleotide>=0) && (bitScore[nucleotide]<scoreSum)) {
				errorSum++;
				if (branchBegin<0)
					branchBegin = j;

				// check if there is read containing > errorCutoff low frequency bases
				if (strlen(strMatrix[0])>j && bitScore[nucleotideToNumber(strMatrix[0][j])]/(float)scoreSum >= cutoffScore){
					xLength = 0;

					usedFlag.purge();		
					
					Len_t deleteNum = 0;
					
					for (Len_t m=1; m<yLength; m++){		
						Len_t remainLen = strlen(strMatrix[m]);
						if (depthsMatrix[m][0] == 0){
							usedFlag.append(0);
							
						}
						else if (remainLen>j && nucleotideToNumber(strMatrix[m][j])!=nucleotide){
							errorArr[m]++;	

							// delete read containing > errorCutoff low frequency bases
							if (errorArr[m] > errorCutoff){	
								for (Len_t n=0; n<remainLen; n++){
									depthsMatrix[m][n] = 0;
								}

								leftReadNum--;

								usedFlag.append(0);		

								deleteNum++;
							}
							else{
								if (remainLen > xLength)
									xLength = remainLen;
								
								usedFlag.append(1);
							}
							
						}
						else {
							if (remainLen > xLength)
								xLength = remainLen;
							usedFlag.append(1);		
						}
						
					}

					if (leftReadNum == 0){  // no read left
						branchPosition = -1;		
						replaceFlag = 0;
						break;
					}
					else{
						if (deleteNum > 0){	
							errorSum = 0;	
							branchBegin = -1;	
							branchPosition = -1;		
							j = 0;

							for (Len_t m=0; m<yLength; m++){
								errorArr[m] = 0;
							}
							
							continue;
						}
					}
					
				}
			}
			
			if ((nucleotide>=0) && (errorSum<=errorCutoff)) {  // no branch
				tStrFinal.writeNucleotideAtPosition(nucleotide,j);
				depthsFinal[j] = bitScore[nucleotide];
			}
			//can not find this nucleotide, it's a branch
			else {
				branchPosition = j;	
				break;
			}

			j++;		
		}

		delete [] errorArr;		
	}


	// function:     branchProcess
	// description: check if there is branch when calling consensus and deal with consensus accordingly
	// input:          contig -- contig on the left side of contig
	//                   reads -- found reads
	//                   offset -- seed's position in contig
	//                   cutoffScore -- cutoff for dominating base
	//                   errorCutoff -- maximum number of position where mismatch appeared in consensus
	//                   branchProcessFlag -- flag to indicate whether process branch or not, true for processing
	// output:        usedFlag -- flag to indicate whether a read is used or not, 1 for used
	void branchProcess(
			Contig& contig, 
			LinkedList<ReadElement> const& reads, 
			LinkedList<Len_t> &usedFlag, 		
			Len_t offset, 
			int& branchPosition, 
			float cutoffScore=1.0, 
			Len_t errorCutoff=2, 
			bool branchProcessFlag=true
			) {
		
		TightString& tStrContig = contig.getTightString();
		ArrayBlock<Len_t>& contigDepths = contig.getDepths();
		
		branchPosition = -1;
		Len_t readsSum = reads.getCount();
		Len_t i,j;
		

		char** strMatrix = new char*[readsSum+1];
		Len_t** depthsMatrix = new Len_t*[readsSum+1];
		
		Len_t maxLengthRemain = 0;
		
		// store reads' remain sequences
		i=1;
		ListElement<ReadElement> const* ptrRead;
		for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
			ReadElement const& readElement = ptrRead->getDatum();
			Len_t readLength = readElement.getLen();
			Len_t readLengthRemain = readLength - overlapParam;
			
			strMatrix[i] = new char[readLengthRemain+1];
			readElement.getSequence(strMatrix[i], readLengthRemain, overlapParam);
			
			depthsMatrix[i] = new Len_t[readLengthRemain];
			Len_t depth = 0;
			if (readElement.getDepth() < maxReadsCount){
				depth = readElement.getDepth();
			}
			
			for(j=0;j<readLengthRemain;j++) {
//				depthsMatrix[i][j] = readElement.getDepth();
				depthsMatrix[i][j] = depth;
			}
			
			if (maxLengthRemain<readLengthRemain) maxLengthRemain=readLengthRemain;
			
			i++;
		}

		
		//store contig's remain sequence
		Len_t contigLengthRemain = tStrContig.getLength() - (offset+overlapParam);

		if (contigLengthRemain>maxLengthRemain)
			contigLengthRemain = maxLengthRemain;
		strMatrix[0] = new char[contigLengthRemain+1];
		tStrContig.readTightStringFragment(offset+overlapParam, offset+overlapParam+contigLengthRemain, strMatrix[0]);
		depthsMatrix[0] = new Len_t[contigLengthRemain];
		for(j=0;j<contigLengthRemain;j++) {
			depthsMatrix[0][j] = contigDepths[offset+overlapParam+j];
		}
		
		TightString tStrFinal(maxLengthRemain);
		Len_t* depthsFinal = new Len_t[maxLengthRemain];

		int replaceFlag = 1;

		matrixProcess(strMatrix, depthsMatrix, tStrFinal, depthsFinal, usedFlag, branchPosition, maxLengthRemain, readsSum+1, replaceFlag, cutoffScore, errorCutoff);
		
		for(i=0;i<readsSum+1;i++) {
			delete [] strMatrix[i];
			delete [] depthsMatrix[i];
		}
		
		delete [] strMatrix;
		delete [] depthsMatrix;
		
		//branch process
		if (branchPosition>=0) {
			
			if (branchProcessFlag) {
				
				if ((Len_t)branchPosition>contigLengthRemain) {
					
					tStrContig.writeTightStringFragment(offset+overlapParam, tStrContig.getLength(), tStrFinal);
					
					Len_t appendLen = branchPosition - contigLengthRemain;
					TightString tStrAppend(appendLen);
					tStrFinal.readTightStringFragment(contigLengthRemain, branchPosition, tStrAppend);				
					contig.append(tStrAppend, appendLen);
				}
				else {
					
					tStrContig.writeTightStringFragment(offset+overlapParam, offset+overlapParam+branchPosition, tStrFinal);
				}
			
				//set depths
				for(j=0; j<(Len_t)branchPosition; j++) {
					contigDepths[offset+overlapParam+j] = depthsFinal[j];
				}
			}
		}
		else if (replaceFlag == 1){		
			
			if (maxLengthRemain>contigLengthRemain) {
				tStrContig.writeTightStringFragment(offset+overlapParam, offset+overlapParam+contigLengthRemain, tStrFinal);
				Len_t appendLen = maxLengthRemain - contigLengthRemain;
				TightString tStrAppend(appendLen);
				tStrFinal.readTightStringFragment(contigLengthRemain, maxLengthRemain, tStrAppend);
				contig.append(tStrAppend, appendLen);
			}
			else{
				tStrContig.writeTightStringFragment(offset+overlapParam, offset+overlapParam+maxLengthRemain, tStrFinal);
			}
			
			//set depths
			for(j=0; j<maxLengthRemain; j++) {
				contigDepths[offset+overlapParam+j] = depthsFinal[j];
			}
		}
		
		delete [] depthsFinal;
	}
	
	void fixedOverlapByPE(Contig& contig, Number_t& endOtherContigPos, Short_Len_t forwardReverseFlag, Len_t maxContigLength) {
		
		TightString& tStrContig = contig.getTightString();
		ArrayBlock< LinkedList<Number_t> >& readPositions = contig.getReadPositions();
		ArrayBlock<Len_t>& contigDepths = contig.getDepths();
		
		LinkedList<Repeat>* pRepeats = new LinkedList<Repeat>;
		
		int originLength = contig.getLength()-reExtendLength*2;
		if (originLength<0) originLength = 0;
		
		//clean depth and read positions
		for (Len_t i=originLength; i<contig.getLength(); i++) {
			
			ListElement<Number_t> const* ptr;
			for (ptr = readPositions[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
				
				ReadElement const& readElement = readAccessor.getRead(ptr->getDatum());
				Len_t depth = readElement.getDepth();
				Len_t len = readElement.getLen();
				len = i+len>contig.getLength() ? contig.getLength()-i : len;
				for (Len_t j=0; j<len; j++) {
					contigDepths[i+j] -= depth;
				}
			}
			
			readPositions[i].purge();
		}
		
		Len_t offset = originLength;
		while (true) {
			
			//get overlap part by param
			TightString tStrOverlap(overlapParam);
			tStrContig.readTightStringFragment(offset, offset+overlapParam, tStrOverlap);
			
			LinkedList<ReadElement> reads;
			readAccessor.getReadsBeginWith(tStrOverlap,reads,true,ReadAccessor::inAll);
			
			if ( (reads.getCount() > 0) && (reads.getCount() < maxReadsCount) ) {
				/*
				//TODO for test
				if ((iContig==3 || iContig==32) && offset>=0) {
					char* seq = tStrContig.readTightString();
					cout << seq <<endl;
					delete [] seq;
					
					cout << "offset:" << offset << endl;
					Len_t iRead = 0;
					ListElement<ReadElement> const* ptrRead;
					for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
						
						ReadElement const& readElement = ptrRead->getDatum();
						
						Len_t readLength = readElement.getLen();
						TightString tStrRead(readLength);
						readElement.getSequence(tStrRead);
						char* seq = tStrRead.readTightString();
						cout << iRead << "+:";
						cout << seq << endl;
						delete [] seq;
						
						TightString tStrReadReverse(readLength);
						tStrRead.readTightStringFragmentAtReverse(0,readLength,tStrReadReverse);
						seq = tStrReadReverse.readTightString();
						cout << iRead << "-:";
						cout << seq << endl;
						delete [] seq;
						
						iRead++;
					}
					cout << endl << endl;
				}*/
				
				//use different insert size to resolve the branch
				Len_t arrayLen = pairInfo.getArrayLen();
				Len_t maxPairsCount = 1;
				Len_t iMaxPair = 0;
				LinkedList<ReadElement> eligibleReads;
				ArrayBlock< LinkedList<ReadElement> > maxPairEligibleReads(arrayLen);
				ListElement<ReadElement> const* ptrRead;
				for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
					
					ReadElement const& readElement = ptrRead->getDatum();
					
					if (readElement.getDepth() >= maxReadsCount)
						continue;
					
					LinkedList<Len_t>* contigPos = new LinkedList<Len_t>[arrayLen];
					getContigPosByPair(readElement, contigPos, contig);
					
					Len_t pairsCount = 0;
					for (Len_t i=0; i<arrayLen; i++) {
						
						PairInfoElement const& pairInfoElement = pairInfo.getArray()[i];
						
						ListElement<Len_t> const* ptr;
						for (ptr = contigPos[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
							
							Len_t pos = ptr->getDatum();
							//check the actual insert size
							Len_t actualInsertSize = offset - pos + readElement.getLen();
							if ( ( actualInsertSize >= (pairInfoElement.insertSize-pairInfoElement.variance) ) && 
									( actualInsertSize <= (pairInfoElement.insertSize+pairInfoElement.variance) ) ) {
								
								pairsCount++;
								maxPairEligibleReads[i].append(readElement);
								if (i>iMaxPair) iMaxPair = i;
								break;
							}
						}
					}
					
					delete [] contigPos;
					
					if (pairsCount > maxPairsCount) {
						maxPairsCount = pairsCount;
						eligibleReads.purge();
						eligibleReads.append(readElement);
					}
					else if (pairsCount == maxPairsCount) {
						eligibleReads.append(readElement);
					}
				}
				
				
				
				if (eligibleReads.getCount() > 0) {
					
					if (selfCircleCheckByPE(eligibleReads, pRepeats, contig))						
						break;
					
					int branchPosition;
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					
					if (branchPosition >= 0) {
						
						tStrContig.setLength(offset+overlapParam+branchPosition);
					}
					
					//save read positions and set depth on this contig
					ListElement<ReadElement> const* ptrRead;
					for (ptrRead = eligibleReads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
						
						ReadElement const& readElement = ptrRead->getDatum();
						
						readPositions[offset].append(readElement.getID());
						
						contig.appendContigPos(readElement, offset);
						
						//set depths
						Len_t depth = readElement.getDepth();
						for(Len_t j=0; j<overlapParam; j++) {
							contigDepths[offset+j] += depth;
						}
					}
				}
				
				
				
				else if (reads.getCount() == 1) {
					
					if (selfCircleCheckByPE(reads, pRepeats, contig))						
						break;
					
					int branchPosition;
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					
					//save read positions and set depth on this contig
					ListElement<ReadElement> const* ptrRead;
					for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
						
						ReadElement const& readElement = ptrRead->getDatum();
						
						readPositions[offset].append(readElement.getID());
						
						contig.appendContigPos(readElement, offset);
						
						//set depths
						Len_t depth = readElement.getDepth();
						for(Len_t j=0; j<overlapParam; j++) {
							contigDepths[offset+j] += depth;
						}
					}
				}

				
				
				else if (reads.getCount() > 1) {
					
					int branchPosition;
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					
					if (branchPosition < 0) {
						
						if (selfCircleCheckByPE(reads, pRepeats, contig))						
							break;
						
						//save read positions and set depth on this contig
						ListElement<ReadElement> const* ptrRead;
						for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
							
							ReadElement const& readElement = ptrRead->getDatum();
							
							readPositions[offset].append(readElement.getID());
							
							contig.appendContigPos(readElement, offset);
							
							//set depths
							Len_t depth = readElement.getDepth();
							for(Len_t j=0; j<overlapParam; j++) {
								contigDepths[offset+j] += depth;
							}
						}
						
					}
					
				}
				
				
			}
			
			if (maxContigLength && tStrContig.getLength() >= maxContigLength) {
				break;
			}
			
			if ( (offset+overlapParam) >= tStrContig.getLength() ) {
				break;
			}
			else
				offset++;
		}
		
		delete pRepeats;
	}
	
	
	
	void assembleUniqueContig(Contig& contig, Len_t maxContigLength) {
		
		switch(overlapMode) {
		
			case fixedOverlapMode:
				fixedOverlapUnique(contig, maxContigLength);
				break;
			default:
				fixedOverlapUnique(contig, maxContigLength);
				break;
		}
	}
	
	void fixedOverlapUnique(Contig& contig, Len_t maxContigLength) {
		
		TightString& tStrContig = contig.getTightString();
		ArrayBlock<Len_t>& contigDepths = contig.getDepths();
		
		Len_t offset = 0;
		while (true) {
			
			//get overlap part by param
			TightString tStrOverlap(overlapParam);
			tStrContig.readTightStringFragment(offset, offset+overlapParam, tStrOverlap);
			
			LinkedList<ReadElement> reads;
			readAccessor.getReadsBeginWith(tStrOverlap,reads,true,ReadAccessor::inAll);
			
			if ( (reads.getCount() > 0) && (reads.getCount() < maxReadsCount) ) {
				
				int branchPosition;
				float cutoffScore=0.8;
				Len_t errorCutoff=2;
				
				if (branchPosition >= 0) {					
					return;
				}
				else {
					
					//set depth on this contig
					ListElement<ReadElement> const* ptrRead;
					for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
						
						ReadElement const& readElement = ptrRead->getDatum();
						
						//set depths
						Len_t depth = readElement.getDepth();
						for(Len_t j=0; j<overlapParam; j++) {
							contigDepths[offset+j] += depth;
						}
					}
				}
			}
			
			if (tStrContig.getLength() >= maxContigLength) {
				return;
			}
			
			if ( (offset+overlapParam) >= tStrContig.getLength() ) {
				return;
			}
			else
				offset++;
		}
	}

};

#endif /*CONTIGASSEMBLER_HPP_*/
