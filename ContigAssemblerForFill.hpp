/**************************************************
*
* ContigAssemblerForFill.hpp
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

#ifndef CONTIGASSEMBLERFORFILL_HPP_
#define CONTIGASSEMBLERFORFILL_HPP_

#include "ContigAssembler.hpp"
#include "ContigTable.hpp"

/*
  Fill gaps in scaffolds in parallel, one thread handling one scaffold.
*/
class ContigAssemblerForFill : public ContigAssembler
{
protected:
	
	ofstream foutFill;
	ContigTable const& contigTable;
	float deviation;
	Len_t endNumLen;
	Len_t mismatchLen;
	
	Len_t minSearchLen;
	Len_t mapReadsRange;
	Len_t actualGapSum;
	Len_t extendGapSum;
	Len_t actualGapCount;
	Len_t finishGapCount;
	
public:
	
	ContigAssemblerForFill(char* _outfile, ofstream& _fout, ReadAccessor& _readAccessor, PairInfo const& _pairInfo, ContigTable const& _contigTable, Len_t _threadSum, float _deviation=0.5, Len_t _endNumLen=10, Len_t _mismatchLen=5, Len_t _maxReadLength=35, Short_Len_t _overlapMode=fixedOverlapMode, Short_Len_t _overlapParam=25) : 

		ContigAssembler(_outfile, _fout, _readAccessor, _pairInfo, _threadSum, _maxReadLength, _overlapMode, _overlapParam), 
		contigTable(_contigTable), 
		deviation(_deviation), 
		endNumLen(_endNumLen), 
		mismatchLen(_mismatchLen), 
		minSearchLen(1000), 
		mapReadsRange(1000), 
		actualGapSum(0), 
		extendGapSum(0), 
		actualGapCount(0), 
		finishGapCount(0)
	{
		reExtendLength = _overlapParam;
	}
	
public:
	
	void assemble() {
		
		cout << ">>>>>>>>>>assembling<<<<<<<<<<" << endl;
		cout << endl;
		time_t total_start_time=time(NULL);
		
		char fillName[MAX_STRING_LEN];
		strcpy(fillName, outfile);
		strcat(fillName,".fill");
		foutFill.open(fillName);
		
		Contig::readAccessor = &readAccessor;
		
		createThread();
		//TODO
//		assembleInThread(0);
		
		foutFill.close();

		cout << "---------------------------------" << endl;
		cout << "actual gap sum: " << actualGapSum << endl;
		cout << "extend gap sum: " << extendGapSum << endl;
		cout << "actual gap count: " << actualGapCount << endl;
		cout << "finish gap count: " << finishGapCount << endl;
		
		cout << "spent total time: " << time(NULL)-total_start_time << "s"<< endl;
		cout << ">>>>>>>>>>assembling finished<<<<<<<<<<" << endl;
		cout << endl;
	}
	
	
	
protected:
	
	struct Arg 
	{
		ContigAssemblerForFill* pContigAssemblerForFill;
		Len_t iThread;
	};
	
	void createThread() {
		
		for (Len_t i=0; i<threadSum; i++) {
			Arg* pArg = new Arg;
			pArg->pContigAssemblerForFill = this;
			pArg->iThread = i;
			pthread_create(&threadsID[i], NULL, (void *(*)(void *))threadFunc, (void*)pArg);
		}
		
		for (Len_t i=0; i<threadSum; i++) {
			pthread_join(threadsID[i], NULL);
		}
	}
	
	static void threadFunc(Arg* pArg) {
		
		ContigAssemblerForFill& contigAssemblerForFill = *(pArg->pContigAssemblerForFill);
		Len_t iThread = pArg->iThread;
		delete pArg;
		
		contigAssemblerForFill.assembleInThread(iThread);
	}
	
	struct GapResultInfo
	{
		Len_t start;
		Len_t end;
		Len_t flag;
		Len_t quality;
		ExtendInfo extendInfo;
		
		GapResultInfo() : 
			start(0), 
			end(0), 
			flag(0), 
			quality(0)
		{
		}
	};

		
	/* assembly in one thread */
	void assembleInThread(Len_t iThread) {
		
		Len_t actualGapSumInThread=0;
		Len_t extendGapSumInThread=0;
		
		Len_t actualGapCountInThread=0;
		Len_t finishGapCountInThread=0;
		
		for (Len_t i=1; i<=contigTable.getContigsSum(); i++) {

			ContigForFill& contig = contigTable.getContigs()[i];

			// get an unused scaffold
			bool isUsed;
			pthread_mutex_lock(&mutexNumberOfContigs);
			isUsed = contig.isUsed();
			if (!isUsed) contig.setUsedFlag();
			pthread_mutex_unlock(&mutexNumberOfContigs);
			if (isUsed) continue;
			
			pthread_mutex_lock(&mutexNumberOfContigs);
			numberOfContigs++;
			cout << "constructing " << numberOfContigs << " scaffold" <<endl;
			pthread_mutex_unlock(&mutexNumberOfContigs);

			
			actualGapSumInThread += contig.getGapSum();
			actualGapCountInThread += contig.getGapCount();
			
			LinkedList<TightString*> const& sequences = contig.getSequences();
			LinkedList<GapInfo> const& gapsList = contig.getGaps();
						
			Len_t contigCount = sequences.getCount();
			Len_t gapCount = gapsList.getCount();
			Contig* contigs = new Contig[contigCount];
			GapInfo* gaps = new GapInfo[gapCount];

			// map reads to contigs
			Len_t j=0;
			ListElement<TightString*> const* ptr;
			for (ptr = sequences.getHead(); ptr != 0; ptr = ptr->getNext()) {
				
				TightString const& seq = *(ptr->getDatum());
				Contig contig(seq);
				if (contig.getLength() > mapReadsRange*2) {
					contig.mapReads(0, mapReadsRange, overlapParam);
					contig.mapReads(contig.getLength()-mapReadsRange, mapReadsRange, overlapParam);
				}
				else {
					contig.mapReads(0, contig.getLength(), overlapParam);
				}
				contigs[j] = contig;
				j++;
			}

			// get gaps
			j=0;
			ListElement<GapInfo> const* ptrGap;
			for (ptrGap = gapsList.getHead(); ptrGap != 0; ptrGap = ptrGap->getNext()) {
				
				GapInfo const& gap = ptrGap->getDatum();
				gaps[j] = gap;
				j++;
			}
			
			Contig* contigsResult = new Contig[contigCount];
			GapInfo* gapsResult = new GapInfo[gapCount];
			
			// fill gaps in forward direction
			fillGaps(contigs, gaps, contigsResult, gapsResult, contigCount, gapCount);
			
			delete [] contigs;
			delete [] gaps;
			Contig* contigsReverse = new Contig[contigCount];
			GapInfo* gapsReverse = new GapInfo[gapCount];
			
			//reverse contigs and gaps
			for (j=0; j<gapCount; j++) {
				
				Contig const& contig = contigsResult[j];
				GapInfo const& gap = gapsResult[j];
				
				Contig contigReverse;
				contig.reverse(contigReverse, 0, contig.getLength());
				
				contigsReverse[contigCount-j-1] = contigReverse;
				
				GapInfo gapReverse;
				gapReverse.length = gap.length;
				gapReverse.extendInfo.leftLen = gap.extendInfo.rightLen;
				gapReverse.extendInfo.rightLen = gap.extendInfo.leftLen;
				gapReverse.isFilled = gap.isFilled;
				gapReverse.quality = gap.quality;
				
				if (!gapReverse.isFilled) {
					
					Len_t sameFlag = 1;	
					for (Len_t j=0; j<mismatchLen; j++) {
						
						if ((int)(contigReverse.getLength() - j) >= (int)endNumLen) {
							
							Number_t numReverse;
							contigReverse.getTightString().readTightStringFragment(j, j+endNumLen, numReverse);
							
							if (!gapReverse.endNumHash.find(numReverse)) {
								gapReverse.endNumHash[numReverse] = j;

								sameFlag = 0;	
							}
							else {
								gapReverse.endNumHash[numReverse] = (Len_t)-1;	
								sameKmer++;
							}
						}
					}

					if (contigReverse.getLength() < endNumLen) {
						sameFlag = 0;
					}
					if (sameFlag == 1) {
						sameEnd++;
					}
				}
				
				gapsReverse[gapCount-j-1] = gapReverse;
			}
			
			for (; j<contigCount; j++) {
				
				Contig const& contig = contigsResult[j];
				Contig contigReverse;
				contig.reverse(contigReverse, 0, contig.getLength());
				contigsReverse[contigCount-j-1] = contigReverse;
			}

			// add unique anchor of next contig to gap releated hash
			for (j=0; j<gapCount; j++) {
				Contig const& contigReverse = contigsReverse[j+1];
				GapInfo & gapReverse = gapsReverse[j];

				if (gapReverse.isFilled) {
					continue;
				}

				Len_t uniqueNum = 0;
				Len_t k=0;
				for ( ; k<mismatchLen; k++) {
					if ((int)(contigReverse.getLength() - k) >= (int)endNumLen) {
						Number_t numReverse;
						contigReverse.getTightString().readTightStringFragment(k, k+endNumLen, numReverse);
							
						if (gapReverse.endNumHash.find(numReverse) && (gapReverse.endNumHash[numReverse] != (Len_t)-1)) {
							uniqueNum = 1;
							break;
						}
					}
					else {
						uniqueNum = 1;
					}
				}

				if (uniqueNum == 1) { // at least an unique anchor was found
					continue;
				}
				else {
					Len_t mismatchNum = 0;
					while ((int)(contigReverse.getLength() - k) >= (int)endNumLen) {
						Number_t numReverse;
						contigReverse.getTightString().readTightStringFragment(k, k+endNumLen, numReverse);
							
						if (!gapReverse.endNumHash.find(numReverse)) {  // an unique anchor was found
							gapReverse.endNumHash[numReverse] = k;							
							break;			
						}
						else {
							gapReverse.endNumHash[numReverse] = (Len_t)-1;
							sameKmer++;
						}

						k++;
					}
				}
			}
								
			delete [] contigsResult;
			delete [] gapsResult;
			contigsResult = new Contig[contigCount];
			gapsResult = new GapInfo[gapCount];
			
			//fill gaps in reversed direction
			fillGaps(contigsReverse, gapsReverse, contigsResult, gapsResult, contigCount, gapCount);

			delete [] contigsReverse;
			delete [] gapsReverse;
			
			//get sequence
			string finalSeqReverse;
			for (j=0; j<gapCount; j++) {
				
				Contig const& contig = contigsResult[j];
//				GapInfo const& gap = gapsResult[j];
				GapInfo & gap = gapsResult[j];
				
				char* strSeq = contig.getTightString().readTightString();
				if (gap.length<0) {
					if ((int)(contig.getLength()+gap.length) > 0)
						finalSeqReverse.append(strSeq, contig.getLength()+gap.length);
				}
				else {
					finalSeqReverse.append(strSeq, contig.getLength());
				}
				delete [] strSeq;
				
				//fill 'N' string
				for (int j=0; j<gap.length; j++) {
					finalSeqReverse += 'N';
				}
				
				if ((!gap.isFilled) && (gap.length==0)) {
					if (NNumber > 0){
						string NSeq(NNumber, 'N');
						finalSeqReverse += NSeq;
						gap.length = NNumber;
					}
				}
			}

			// get last contig sequence
			for (; j<contigCount; j++) {
				
				Contig const& contig = contigsResult[j];
				char* strSeq = contig.getTightString().readTightString();
				finalSeqReverse.append(strSeq, contig.getLength());
				delete [] strSeq;
			}
			
			char* finalSeq = reverseComplement(finalSeqReverse.c_str(),finalSeqReverse.size());
			
			//get extend info
			GapResultInfo gapResultInfo;
			LinkedList<GapResultInfo> gapResultsInfo;
			Len_t pos = 0;
			for (j=0; j<gapCount; j++) {
				
				Contig const& contig = contigsResult[j];
				GapInfo const& gap = gapsResult[j];
				
				pos += contig.getLength();
//				if (gap.length<0) {
//					pos += gap.length;
//				}
				
				gapResultInfo.end = finalSeqReverse.size() - (pos-gap.extendInfo.leftLen);
				gapResultInfo.start = gapResultInfo.end - (gap.extendInfo.leftLen+gap.length+gap.extendInfo.rightLen);
				gapResultInfo.extendInfo.leftLen = gap.extendInfo.rightLen;
				gapResultInfo.extendInfo.rightLen = gap.extendInfo.leftLen;
				gapResultInfo.flag = gap.isFilled;
				gapResultInfo.quality = gap.quality;
				
				gapResultsInfo.prepend(gapResultInfo);
				
				extendGapSumInThread += gap.extendInfo.leftLen+gap.extendInfo.rightLen;
				if (gap.isFilled) {
					finishGapCountInThread++;
				}
				
				pos += gap.length;
			}
			
			delete [] contigsResult;
			delete [] gapsResult;
			
			pthread_mutex_lock(&mutexOutput);
			//output sequence, fasta format
			fout << contig.getName() << endl;
//			fout << finalSeq << endl;
			
//			Len_t j;
			Len_t finalSeqLen = strlen(finalSeq);
			for (j=0; j<finalSeqLen;) {
				
				fout << finalSeq[j];
				j++;
				if (j%100==0)
					fout << endl;
			}
			if (j%100)
				fout << endl;
			
			//output extend info
			foutFill << contig.getName() << endl;

			ptrGap = gapsList.getHead();
			ListElement<GapResultInfo> const* ptrGapResult;
			for (ptrGapResult = gapResultsInfo.getHead(); ptrGapResult != 0; ptrGapResult = ptrGapResult->getNext()) {
				
				gapResultInfo = ptrGapResult->getDatum();
//				foutFill << gapResultInfo.start << "\t" << gapResultInfo.end << "\t" << gapResultInfo.extendInfo.leftLen << "\t" <<gapResultInfo.extendInfo.rightLen << "\t" <<gapResultInfo.flag << "\t" <<gapResultInfo.quality << endl;
				foutFill << gapResultInfo.start << "\t" << gapResultInfo.end << "\t" << gapResultInfo.extendInfo.leftLen << "\t" <<gapResultInfo.extendInfo.rightLen << "\t" <<gapResultInfo.flag << "\t" <<gapResultInfo.quality << "\t" <<ptrGap->getDatum().length << "\t" << gapResultInfo.end - gapResultInfo.start << endl;

				ptrGap = ptrGap->getNext();
			}
			
			pthread_mutex_unlock(&mutexOutput);
			
			delete [] finalSeq;
		}
		
		pthread_mutex_lock(&mutexNumberOfContigs);
		actualGapSum += actualGapSumInThread;
		extendGapSum += extendGapSumInThread;
		actualGapCount += actualGapCountInThread;
		finishGapCount += finishGapCountInThread;
		pthread_mutex_unlock(&mutexNumberOfContigs);
	}

	// fucntion: 	 fillGaps
	// description: fill gaps in scaffold one by one in direction of left to right
	// input:		 contigs -- contigs before gap filling
	//			 gaps -- gaps before gap filling
	//			 contigCount -- contig number
	//			 gapCount -- gap number
	// output:	 contigsResult -- contigs after gap filling
	//			 gapsResult -- gaps after gap filling
	void fillGaps(
			Contig const* contigs, 
			GapInfo const* gaps, 
			Contig* contigsResult, 
			GapInfo* gapsResult, 
			Len_t contigCount, 
			Len_t gapCount
			) {
		
		Len_t i;
		for (i=0; i<gapCount; i++) {
			
			Contig const& contig = contigs[i];
			GapInfo const& gap = gaps[i];
			Len_t j = i< contigCount-1 ? i+1: i;
			Contig const& nextContig = contigs[j];

			Len_t notNCtgLen = 0;		
			if (!gap.isFilled ){		
				
				Contig contigFill;
				
				
				if ( (contig.getLength() < mapReadsRange) && i>0 ) {  // append previous contig at the beginning of current contig
					
					//calculate the length
					int sum=0;
					sum += contig.getLength();
					Len_t redundantSum = 0;
					bool appendGap = false;
					int j;

					notNCtgLen += contig.getLength();		
					Len_t foundNotFilledGap = 0;			
					
					for(j=i-1; j>=0; j--) {

						if (!gapsResult[j].isFilled){
							foundNotFilledGap = 1;
						}
						
						sum += gapsResult[j].length;
						if (sum>=(int)mapReadsRange) {
							redundantSum = sum - mapReadsRange;
							appendGap = true;
							break;
						}
						if (foundNotFilledGap == 0){
							if (gapsResult[j].length < 0){
								if ((int)contigsResult[j].getLength() + gapsResult[j].length > 0){
									notNCtgLen += contigsResult[j].getLength() + gapsResult[j].length;
								}
							}
							else{
							
								notNCtgLen += contigsResult[j].getLength();
							}
						}
						sum += contigsResult[j].getLength();
						if (sum>=(int)mapReadsRange) {
							redundantSum = sum - mapReadsRange;
							appendGap = false;
							break;
						}
					}
					if (j<0) j=0;
					
					//then append
					if (appendGap) {
						if (!gapsResult[j].isFilled) {
							TightString tStrAppend(gapsResult[j].length-redundantSum);
							contigFill.append(tStrAppend, gapsResult[j].length-redundantSum);
							redundantSum = 0;
						}
						j++;
					}
					
					for (; j<(int)i; j++) {
						
						if (gapsResult[j].length<0) {
							if ((int)(contigsResult[j].getLength()+gapsResult[j].length-redundantSum) > 0)
								contigFill.append(contigsResult[j], redundantSum, contigsResult[j].getLength()+gapsResult[j].length-redundantSum);
						}
						else {
							contigFill.append(contigsResult[j], redundantSum, contigsResult[j].getLength()-redundantSum);
						}
						
						if (redundantSum) redundantSum = 0;
						
						if (!gapsResult[j].isFilled) {
							TightString tStrAppend(gapsResult[j].length);
							contigFill.append(tStrAppend, gapsResult[j].length);
						}
					}
					contigFill.append(contig);
				}
				else {

//					contigFill.append(contig);
					if (contig.getLength() < mapReadsRange) {
						contigFill.append(contig);

						notNCtgLen = contig.getLength();		
					}
					else {
						contigFill.append(contig, contig.getLength()-mapReadsRange, mapReadsRange);

						notNCtgLen = mapReadsRange;	
					}
				}
				if (notNCtgLen > mapReadsRange) {
					notNCtgLen = mapReadsRange;
				}
				
				if (notNCtgLen >= overlapParam) {
				
					Contig gapContig;
					GapInfo gapResult;
					fillGap(contigFill, gap, gapContig, gapResult, notNCtgLen, nextContig);
				
					Contig contigResult(contig,gapContig);
					contigsResult[i] = contigResult;
				
					gapResult.quality = gap.quality;
					for (Len_t j=0; j<gapContig.getLength(); j++) {
						gapResult.quality += gapContig.getQuality()[j];
					}
				
					gapResult.extendInfo.leftLen = gapContig.getLength();
					gapResult.extendInfo.rightLen = gap.extendInfo.rightLen;
					gapsResult[i] = gapResult;

				}
			}
			else if (gap.isFilled || notNCtgLen < overlapParam) {		
				contigsResult[i] = contig;
				
				GapInfo gapResult;
				gapResult.length = gap.length;
				gapResult.extendInfo = gap.extendInfo;
				gapResult.isFilled = gap.isFilled;
				gapResult.quality = gap.quality;
				gapsResult[i] = gapResult;
			}
		}

		// add the last contig
		for (; i<contigCount; i++) {
			
			contigsResult[i] = contigs[i];
		}
	}

	// function:     fillGap
	// description: fill one gap with selected mode. Currenty there is only one mode provided.
	// input:		 contig -- contig on the left side of gap
	//			 gap -- gap to fill
	//			 realContigLen -- length of none N sequence in contig on the left side of gap
	//			 nextcontig -- contig on the right side of gap
	// output:	 gapContig -- extend sequence in gap
	//			 gapResult -- filled gap result
	void fillGap(Contig& contig, GapInfo const& gap, Contig& gapContig, GapInfo& gapResult, Len_t realContigLen, Contig const& nextcontig) {	
		switch(overlapMode) {
		
			case fixedOverlapMode:
				fixedOverlapForFillByPE(contig, gap, gapContig, gapResult, realContigLen, nextcontig);
				break;
			default:
				fixedOverlapForFillByPE(contig, gap, gapContig, gapResult, realContigLen, nextcontig);	
				break;
		}
	}

	// function:     checkBranch
	// description: use all reads that could be gathered to solve branch
	// input:          branchPosition -- branch position in extend sequence
	//                   originLength -- length of contig region whose related information would not be changed during extension, 
	//                                          it is also used as the original start position of extension
	//                   offset -- seed's position in contig
	//                   arraySize -- length of contig region with gathered reads information
	//                   arrayPtr -- pointer to current position of found reads array
	//                   offsetStop -- furthest offset in contig before which reads have beed gathered
	//                   contig -- contig on the left side of gap
	//                   foundFlag -- flag to indicate whether reads are gathered, 1 for yes
	//                   foundReads -- gathered reads
	//                   eligibleFlags -- flag to indicate whether gathered reads have paired-end support, 1 for yes
	//                   usedFlags -- flag to indicate whether reads are used, 1 for yes
	//                   allInsertSizes -- insert size of gathered reads
	//                   preBranchPositions -- branch position before current one
	//                   unreliableBranchPos -- branch position where error might happen with high probability
	//                   errorCutoff -- maximum number of mismatch
	//                   isEligible -- 1 if current reads have paired-end support
	void checkBranch(int branchPosition, int originLength, Len_t offset, Len_t arraySize, Len_t arrayPtr, Len_t& offsetStop, Contig& contig, 
						Len_t *& foundFlag, LinkedList<ReadElement> *& foundReads, LinkedList<Len_t> *& eligibleFlags, LinkedList<Len_t> *& usedFlags, 
						LinkedList<Len_t> **& allInsertSizes, LinkedList<Len_t>& preBranchPositions, LinkedList<Len_t> &unreliableBranchPos, Len_t errorCutoff, bool isEligible){

		Len_t branchPositionInCtg = offset + overlapParam + branchPosition;
		TightString& tStrContig = contig.getTightString();
		ArrayBlock< LinkedList<Number_t> >& readPositions = contig.getReadPositions();
		ArrayBlock<Len_t>& contigDepths = contig.getDepths();
		Len_t localPtr;
		Len_t backoff;
		Len_t remainLen = offsetStop - offset;  // length of region with gathered reads before branch position
		
		if ((int)remainLen > branchPosition)
			remainLen = branchPosition;
		
		if (offset >= originLength + maxReadLength - (overlapParam + branchPosition) - 1){  
			// the branch position's distance to original start position of extension is large enough
			// so that the furthest reads could be retrieved are those who end at branch position
			backoff = maxReadLength - (overlapParam + branchPosition) - 1;
		}
		else{
			// the furthest reads could be retrieved are those who start at the original start position of extension
			backoff = offset - originLength;
		}

		
		if (arrayPtr >= backoff){
			localPtr = arrayPtr - backoff;
		}
		else{
			localPtr = arraySize - (backoff - arrayPtr);
		}

		Len_t backupLocalPtr = localPtr;		

		char* subCtgSeq = new char[contig.getLength()-(offset-backoff)+1];  // contig region starting from the positon where furthest reads start
		contig.getTightString().readTightStringFragment(offset-backoff, contig.getLength(), subCtgSeq);	

		char *readSeq = new char[maxReadLength+1];	

		Len_t beforeDepths[4] = {0, 0, 0, 0};
		Len_t totalDepths[4] = {0, 0, 0, 0};
		Len_t eligibleDepths[4] = {0, 0, 0, 0};

		LinkedList<Len_t> *localUsedFlag = new LinkedList<Len_t>[backoff+1];

		LinkedList<Len_t>* insertSizes4Bases = new LinkedList<Len_t>[4];	

		// deal with gathered reads starting before current branch position
		Len_t stopPos = backoff + remainLen;
		for (Len_t m=0; m<=stopPos; m++) {
			if (foundFlag[localPtr] == 1) {

				ListElement<ReadElement> const *readPtr = foundReads[localPtr].getHead();
				ListElement<Len_t> const *eligiblePtr = eligibleFlags[localPtr].getHead();

				LinkedList<Len_t> *preInsertSizes = allInsertSizes[localPtr];		
				Len_t preInsSubcript = 0;		
				for (; readPtr!=0; readPtr=readPtr->getNext()) {
					ReadElement const& readElement = readPtr->getDatum();

					if (readElement.getDepth() < maxReadsCount ) {
						if (m <= backoff || eligiblePtr->getDatum() > 0 || foundReads[localPtr].getCount() == 1) {
							int readLen = readElement.getLen();
							int endPos = backoff - m + overlapParam + branchPosition;
							
							readElement.getSequence(readSeq, readLen, 0);	
							Len_t used = 1;

							if (readLen > endPos) {  // read goes across branch position
								Len_t errorSum = 0;
								// check the alignment between read and contig
								for (int pos=overlapParam; pos<endPos; pos++) {
									if (subCtgSeq[m+pos] != readSeq[pos]){
										errorSum++;
										if (errorSum > errorCutoff) {
											used = 0;
											break;
										}
									}
								}
							}

							if (m <= backoff)
								localUsedFlag[m].append(used);		
							
							
							if (used == 1 && readLen > endPos){
								Nucleotide_t cmpBase = nucleotideToNumber(readSeq[endPos]);

								if ((isEligible || m > backoff) && eligiblePtr->getDatum() > 0){
									eligibleDepths[(int)cmpBase]++;

									// make statistics of insert size information
									for (ListElement<Len_t> const* ptrIns=preInsertSizes[preInsSubcript].getHead(); ptrIns!=0; ptrIns=ptrIns->getNext()){
										if (!insertSizes4Bases[(int)cmpBase].find(ptrIns->getDatum())){
											insertSizes4Bases[(int)cmpBase].append(ptrIns->getDatum());
												
										}
									}
								}

								if (m <= backoff)
									beforeDepths[(int)cmpBase]++;

								totalDepths[(int)cmpBase]++;

							}
						}

						eligiblePtr = eligiblePtr->getNext();

					}

					preInsSubcript++;
					

				}
			}

			localPtr ++;
			localPtr = localPtr%arraySize;
		}
		
		// find reads starting from contig region without gathered reads
		TightString tStrOverlap(overlapParam);
		ListElement<ReadElement> const* ptrRead;
		Len_t arrayLen = pairInfo.getArrayLen();
		for (Len_t m=remainLen; m<=branchPosition; m++){

			Len_t currentPos = offset+m+1;
			tStrContig.readTightStringFragment(currentPos, currentPos+overlapParam, tStrOverlap);

			LinkedList<ReadElement>& extendReads = foundReads[localPtr];
			LinkedList<Len_t> *&currentInsertSizes = allInsertSizes[localPtr];
			if (foundFlag[localPtr] == 1){
				for (Len_t i=0; i<extendReads.getCount(); i++){
					currentInsertSizes[i].purge();
				}
				delete [] currentInsertSizes;
			}
			extendReads.purge();
			currentInsertSizes = NULL;
			readAccessor.getReadsBeginWith(tStrOverlap,extendReads,true,ReadAccessor::inAll);

			LinkedList<Len_t> &eligibleFlag = eligibleFlags[localPtr];
			eligibleFlag.purge();

			foundFlag[localPtr] = 0;

			LinkedList<Len_t> &usedFlag = usedFlags[localPtr];
			usedFlag.purge();

			if ( (extendReads.getCount() > 0) && (extendReads.getCount() < maxReadsCount) ) {
				
				foundFlag[localPtr] = 1;

				currentInsertSizes = new LinkedList<Len_t>[extendReads.getCount()];	

				Len_t insSubscript = 0;

				for (ptrRead = extendReads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
	
					ReadElement const& readElement = ptrRead->getDatum();
	
					if (readElement.getDepth() >= maxReadsCount){
						usedFlag.append(0);
						insSubscript++;
						continue;
					}

					LinkedList<Len_t>* contigPos = new LinkedList<Len_t>[arrayLen];
					getContigPosByPair(readElement, contigPos, contig);
		
					Len_t pairsCount = 0;
					Len_t PECount = 0;		

					LinkedList<Len_t> insertSizes;	
	
					for (Len_t i=0; i<arrayLen; i++) {	//process rank by rank
						Len_t pairsCountFlag = 1;		
						PairInfoElement const& pairInfoElement = pairInfo.getArray()[i];
		
						ListElement<Len_t> const* ptr;
						for (ptr = contigPos[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
			
							Len_t pos = ptr->getDatum();
							//check the actual insert size
							Len_t actualInsertSize = currentPos - pos + readElement.getLen();
							if ( ( actualInsertSize >= (pairInfoElement.insertSize-pairInfoElement.variance) ) && 
								( actualInsertSize <= (pairInfoElement.insertSize+pairInfoElement.variance) ) ) {
								if (pairsCountFlag == 1){
									pairsCount++;
									pairsCountFlag = 0;
								}
								PECount++;		

								currentInsertSizes[insSubscript].append(pairInfoElement.insertSize);

								Len_t existsFlag = 0;
								for (ListElement<Len_t> const* ptrIns=insertSizes.getHead(); ptrIns!=0; ptrIns=ptrIns->getNext()){
									if (ptrIns->getDatum() == pairInfoElement.insertSize){
										existsFlag = 1;
										break;
									}
								}
								if (existsFlag == 0)
									insertSizes.append(pairInfoElement.insertSize);
							}
							
						}
					}

					for (Len_t n=0; n<arrayLen; n++)
						contigPos[n].purge();
					delete [] contigPos;

					eligibleFlag.append(pairsCount);		
					usedFlag.append(1);	
					insSubscript++;


					if (pairsCount > 0 || extendReads.getCount() == 1){

						int readLen = readElement.getLen();
						
						int endPos = overlapParam+branchPosition-1-m;
						readElement.getSequence(readSeq, readLen, 0);	

					
						Len_t used = 1;

						if (readLen > endPos){  // read goes across branch position
							Len_t errorSum = 0;	
							// check alignment between read and contig
							for (int pos=overlapParam; pos<endPos; pos++){
								if (subCtgSeq[backoff+m+1+pos] != readSeq[pos]){
									errorSum++;
									if (errorSum > errorCutoff){
										
										used = 0;
										break;
									}
								}
							}
						}

						if (used == 1 && readLen > endPos){

							Nucleotide_t cmpBase = nucleotideToNumber(readSeq[endPos]);
							if (pairsCount > 0){
									eligibleDepths[(int)cmpBase]++;

									for (ListElement<Len_t> const* ptrIns=insertSizes.getHead(); ptrIns!=0; ptrIns=ptrIns->getNext()){
										if (!insertSizes4Bases[cmpBase].find(ptrIns->getDatum()))
											insertSizes4Bases[cmpBase].append(ptrIns->getDatum());
									}

							}
							totalDepths[(int)cmpBase]++;
						}
					}
				}
			}

			localPtr++;
			localPtr = localPtr%arraySize;
		}

		// deal with branch position
		int bestBase = -1;
		int secondBestBase = -1;
		Len_t maxEligibleDepth = 0;
		Len_t secondEligibleDepth = 0;
		Len_t maxDiffInsNum = 0;
		Len_t secondDiffInsNum = 0;
		for (Len_t m=0; m<4; m++){

			if (insertSizes4Bases[m].getCount() > maxDiffInsNum) {	
				secondBestBase = bestBase;
				secondDiffInsNum = maxDiffInsNum;

				bestBase = (int)m;

				maxDiffInsNum = insertSizes4Bases[m].getCount();	
			}

			else if (insertSizes4Bases[m].getCount() > secondDiffInsNum) {	
				secondBestBase = (int)m;

				secondDiffInsNum = insertSizes4Bases[m].getCount();	
			}

			insertSizes4Bases[m].purge();		
		}

		delete [] insertSizes4Bases;			

		if (maxDiffInsNum < secondDiffInsNum+2) {  // the dominance of best base is not large enough
			unreliableBranchPos.append(branchPositionInCtg);
		}

		if (maxDiffInsNum == secondDiffInsNum) {
			// select the base having most paired-end reads support
			bestBase = -1;
			secondBestBase = -1;
			for (Len_t m=0; m<4; m++) {
				if (eligibleDepths[m] > maxEligibleDepth) {
					secondBestBase = bestBase;
					secondEligibleDepth = maxEligibleDepth;
				
					bestBase = (int)m;
					maxEligibleDepth = eligibleDepths[m];
				}
				else if (eligibleDepths[m] >secondEligibleDepth) {
					secondBestBase = (int)m;
					secondEligibleDepth = eligibleDepths[m];
				}
				
			}
							
		}

		int preBase = tStrContig.getNucleotide(branchPositionInCtg);

		if (maxDiffInsNum == secondDiffInsNum && maxEligibleDepth == secondEligibleDepth) {

			if (isEligible) {
				if (secondBestBase == preBase)	 // keep the base appeared in contig
					bestBase = secondBestBase;
			}
			else {  // select the base having most read support
				Len_t maxDepth = 0;
				Len_t secondDepth = 0;
				for (Len_t m=0; m<4; m++) {
					if (totalDepths[m] > maxDepth) {
						secondDepth = maxDepth;
						secondBestBase = bestBase;

						maxDepth = totalDepths[m];
						bestBase = m;
					}
					else if (totalDepths[m] > secondDepth) {
						secondBestBase = m;
						secondDepth = totalDepths[m];
					}
				}
			}
				
		}

		if (preBase != bestBase) {
			// replace the previous base
			tStrContig.writeNucleotideAtPosition(bestBase, branchPositionInCtg);
			tStrContig.setLength(offset+overlapParam+branchPosition+1);
			offsetStop = offset + branchPosition;
			
			// update reads aligned to contig
			localPtr = backupLocalPtr;
			for (Len_t m=0; m<=backoff; m++) {
				if (foundFlag[localPtr] == 1) {
					ListElement<Len_t> const *usedPtr = localUsedFlag[m].getHead();
					ListElement<ReadElement> const *readPtr = foundReads[localPtr].getHead();
					for (; readPtr!=0; readPtr=readPtr->getNext()){
						ReadElement const& readElement = readPtr->getDatum();
						
						if (readElement.getDepth() < maxReadsCount ){
							int readLen = readElement.getLen();
							int endPos = readLen - 1 -(int)m - (maxReadLength - overlapParam -backoff-branchPosition-1);
							if (usedPtr->getDatum() == 1 && readLen > endPos){
								Len_t currentPosInCtg = offset-backoff+m;
								readPositions[currentPosInCtg].append(readElement.getID());
		
								contig.appendContigPos(readElement, currentPosInCtg);
							}
							usedPtr = usedPtr->getNext();
						}						
					}
				}

				localPtr++;
				localPtr = localPtr%arraySize;
			}

			// update all branch positions information of contig
			ListElement<Len_t> const* branchPtr = preBranchPositions.getHead();
			while (branchPtr != 0){
				if (branchPtr->getDatum() > branchPositionInCtg){
					Len_t const& position = branchPtr->getDatum();
					branchPtr = branchPtr->getNext();
					preBranchPositions.extract(position);
				}
				else{
					branchPtr = branchPtr->getNext();
				}
			}

			// update all unreliable branch positions information of contig
			branchPtr = unreliableBranchPos.getHead();
			while (branchPtr != 0){
				if (branchPtr->getDatum() > branchPositionInCtg){
					Len_t const& position = branchPtr->getDatum();
					branchPtr = branchPtr->getNext();
					unreliableBranchPos.extract(position);
				}
				else{
					branchPtr = branchPtr->getNext();
				}
			}
			
		}
		else{
			offsetStop = offset + branchPosition + 1;
		}
		
		for (Len_t m=0; m<backoff+1; m++){
			localUsedFlag[m].purge();
		}

		delete [] localUsedFlag;
		delete [] subCtgSeq;
		delete [] readSeq;
			
		for (Len_t pos=offset-backoff; pos<branchPositionInCtg; pos++){
			if (contigDepths[pos] < beforeDepths[bestBase])
				contigDepths[pos] = beforeDepths[bestBase];
		}

		contigDepths[branchPositionInCtg] = beforeDepths[bestBase];

	}

	// function:     fixedOverlapForFillByPE
	// description: iterately extend contig toward gap.
	//                   1) get reads starting with anchor fetched from end of contig.
	//                   2) extend contig with different strategy  according to the chacteristic of found reads
	//                   3) check overlap between extend contig and next contig
	// input:          contig -- contig on the left side of gap
	//                   gap -- gap to fill
	//                   realContigLen -- length of none N sequence in contig on the left side of gap
	//                   nextContig -- contig on the right side of gap
	//output:         gapContig -- extend sequence in gap
	//                   gapResult -- filled gap result
	void fixedOverlapForFillByPE(Contig& contig, GapInfo const& gap, Contig& gapContig, GapInfo& gapResult, Len_t realContigLen, Contig const& nextContig) {	
		
		Len_t contigLenOld = contig.getLength();
		Len_t contigLen = contig.getLength();
		Len_t nextCtgLen = nextContig.getLength();
		
		int start = gap.length - (int)(gap.length*deviation) + (int)contigLenOld;
		int end = gap.length + (int)(gap.length*deviation) + (int)contigLenOld;
		if ( end-start < (int)minSearchLen ) {
			
			start = contigLenOld - overlapParam;
			
			if (start<0) 
				start = 0;
			end = start + minSearchLen;
		}

		
		TightString& tStrContig = contig.getTightString();
		ArrayBlock< LinkedList<Number_t> >& readPositions = contig.getReadPositions();
		ArrayBlock<Len_t>& contigDepths = contig.getDepths();
		ArrayBlock<Len_t>& contigQuality = contig.getQuality();
		
		LinkedList<Repeat>* pRepeats = new LinkedList<Repeat>;
		

		int realStart = realContigLen - maxReadLength + 1;	
		if (realStart < 0) realStart = 0;						
		int originLength = contigLenOld - realContigLen + realStart;  // length of contig region whose related information would not be changed during extension
//		if (originLength<0) originLength = 0;

              // length of contig region whin which no branch is allowed
		Len_t forbiddenBranchLen = contigLenOld - originLength;	

		Len_t compensatedDepths[forbiddenBranchLen];
		for (Len_t i=0; i<forbiddenBranchLen; i++){
			compensatedDepths[i] = 0;
		}
		
		//clean depth and read positions
		for (Len_t i=originLength; i<contig.getLength(); i++) {
			
			ListElement<Number_t> const* ptr;
			for (ptr = readPositions[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
				
				ReadElement const& readElement = readAccessor.getRead(ptr->getDatum());
				Len_t depth = readElement.getDepth();
				Len_t len = readElement.getLen();
				len = i+len>contig.getLength() ? contig.getLength()-i : len;
				for (Len_t j=0; j<len; j++) {
					if (contigDepths[i+j] < depth){
					
						contigDepths[i+j] = 0;
						
					}
					else{
						contigDepths[i+j] -= depth;
					}
				}
			}
			
			readPositions[i].purge();	//clear mapped reads with start position i
		}

		// make sure bases within this contig region have enough depth to be considered as dominated bases
		for (Len_t i=originLength; i<contigLenOld; i++){
			if (contigDepths[i] < 10){
				
				compensatedDepths[i-originLength] = 10 -contigDepths[i];
				contigDepths[i] = 10;
			}
		}


		Len_t arraySize = (maxReadLength - overlapParam + 1)*2;
		LinkedList<ReadElement> *foundReads = new LinkedList<ReadElement>[arraySize];
		LinkedList<Len_t> *eligibleFlags = new LinkedList<Len_t>[arraySize];
		LinkedList<Len_t> **allInsertSizes = new LinkedList<Len_t>*[arraySize];	
		LinkedList<Len_t> *usedFlags = new LinkedList<Len_t>[arraySize];

		Len_t *foundFlag = new Len_t[arraySize];
		Len_t arrayPtr = 0;
		for (arrayPtr=0; arrayPtr<arraySize; arrayPtr++){
			foundFlag[arrayPtr] = 0;

			allInsertSizes[arrayPtr] = NULL;
		}	

		arrayPtr = 0;

		LinkedList<Len_t> preBranchPositions;		
		LinkedList<Len_t> unreliableBranchPos;		

		Len_t offset = originLength;
		Len_t offsetStop = offset;		
		while (true) {
			if (offsetStop < offset){
				offsetStop = offset;
			}

			Len_t extendFlag = 0;		

			//get seed sequence from contig
			TightString tStrOverlap(overlapParam);
			tStrContig.readTightStringFragment(offset, offset+overlapParam, tStrOverlap);	//in binary format
			
			LinkedList<ReadElement> &reads = foundReads[arrayPtr];		

			LinkedList<Len_t> *&currentInsertSizes = allInsertSizes[arrayPtr];

			LinkedList<Len_t> &eligibleFlag = eligibleFlags[arrayPtr];
			LinkedList<Len_t> &usedFlag = usedFlags[arrayPtr];


			if (offsetStop == offset){
			
				if (foundFlag[arrayPtr] == 1){
					for (Len_t i=0; i<reads.getCount(); i++){
						currentInsertSizes[i].purge();
					}
					delete [] currentInsertSizes;
				}
				currentInsertSizes = NULL;

				reads.purge();		
				readAccessor.getReadsBeginWith(tStrOverlap,reads,true,ReadAccessor::inAll);

				foundFlag[arrayPtr] = 0;
				eligibleFlag.purge();				
				usedFlag.purge();
			}
			
			
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

				Len_t maxPairsCount = 0;		
				Len_t iMaxPair = 0;	// the maximum rank of supporting paired-end reads
				LinkedList<ReadElement> eligibleReads;
			
				if (offsetStop == offset){
					foundFlag[arrayPtr] = 1;

					currentInsertSizes = new LinkedList<Len_t>[reads.getCount()];	

					Len_t insSubscript = 0;		
					
					//use different insert size to resolve the branch
					Len_t arrayLen = pairInfo.getArrayLen();
					ArrayBlock< LinkedList<ReadElement> > maxPairEligibleReads(arrayLen);
					ListElement<ReadElement> const* ptrRead;
					for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
						
						ReadElement const& readElement = ptrRead->getDatum();
						
						if (readElement.getDepth() >= maxReadsCount) {
							insSubscript++;		
							usedFlag.append(0);
							continue;
						}
						LinkedList<Len_t>* contigPos = new LinkedList<Len_t>[arrayLen];

						getContigPosByPair(readElement, contigPos, contig);
						
						Len_t pairsCount = 0;
						Len_t PECount = 0;		
						
						for (Len_t i=0; i<arrayLen; i++) {	//process rank by rank
							Len_t pairsCountFlag = 1;		
							PairInfoElement const& pairInfoElement = pairInfo.getArray()[i];
							
							ListElement<Len_t> const* ptr;
							for (ptr = contigPos[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
								
								Len_t pos = ptr->getDatum();
								//check the actual insert size
								Len_t actualInsertSize = offset - pos + readElement.getLen();
								if ( ( actualInsertSize >= (pairInfoElement.insertSize-pairInfoElement.variance) ) && 
										( actualInsertSize <= (pairInfoElement.insertSize+pairInfoElement.variance) ) ) {
									if (pairsCountFlag == 1) {
										pairsCount++;
										pairsCountFlag = 0;
									}
									PECount++;		

									maxPairEligibleReads[i].append(readElement);
									if (i>iMaxPair) iMaxPair = i;

									currentInsertSizes[insSubscript].append(pairInfoElement.insertSize);	
										
								}
							}
						}
						
						for (Len_t m=0; m<arrayLen; m++)
							contigPos[m].purge();
						delete [] contigPos;

						eligibleFlag.append(pairsCount);		
						usedFlag.append(1);			
						
						if (pairsCount > maxPairsCount) {
							maxPairsCount = pairsCount;
							eligibleReads.purge();
							eligibleReads.append(readElement);
						}
						else if (maxPairsCount > 0 && pairsCount == maxPairsCount) {
							eligibleReads.append(readElement);
						}

						insSubscript++;	
					}
				}
				else {
					Len_t insSubscript = 0;
					Len_t pairsCount = 0;
					ListElement<ReadElement> const* ptrRead;
					ListElement<Len_t> const* eligiblePtr = eligibleFlag.getHead();
					for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {
					
						ReadElement const& readElement = ptrRead->getDatum();
					
						if (readElement.getDepth() >= maxReadsCount) {
							insSubscript++;		
							continue;
						}

						pairsCount = eligiblePtr->getDatum();
						if (pairsCount > maxPairsCount) {
							maxPairsCount = pairsCount;
							eligibleReads.purge();
							eligibleReads.append(readElement);
						}
						else if (maxPairsCount > 0 && pairsCount == maxPairsCount) {
							eligibleReads.append(readElement);
						}
						insSubscript++;
						eligiblePtr = eligiblePtr->getNext();
					}
				}
							
				if (eligibleReads.getCount() > 0) {	

					
//					if (selfCircleCheckByPE(eligibleReads, pRepeats, contig))						
//						break;
					int branchPosition;
			
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					
					branchProcess(contig, eligibleReads,  usedFlag, offset, branchPosition, cutoffScore, errorCutoff);
										
					if (branchPosition >= 0) {  // found branch

						extendFlag = 1;
						
						Len_t branchPositionInCtg = offset + overlapParam + branchPosition;

						if ((contig.getLength() > branchPositionInCtg) && (!preBranchPositions.find(branchPositionInCtg)) && (branchPositionInCtg>=contigLenOld)){	

							preBranchPositions.append(branchPositionInCtg);

//							if (stopPtr == arraySize){
//								stopPtr = arrayPtr;
//							}

							checkBranch(branchPosition, originLength, offset, arraySize, arrayPtr, offsetStop, contig, foundFlag, 
										foundReads, eligibleFlags, usedFlags, allInsertSizes, preBranchPositions, unreliableBranchPos, errorCutoff, true);
						
						}
						
					}
					else{
						//save read positions and set depth on this contig
						ListElement<ReadElement> const* ptrRead;
						ListElement<Len_t> const* ptrUsed=usedFlag.getHead();
						for (ptrRead = eligibleReads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {

							if (ptrUsed->getDatum() == 1){

								extendFlag = 1;
						
								ReadElement const& readElement = ptrRead->getDatum();
						
								readPositions[offset].append(readElement.getID());
						
								contig.appendContigPos(readElement, offset);
						
								//set depths
								Len_t depth = readElement.getDepth();
								for(Len_t j=0; j<overlapParam; j++) {
									contigDepths[offset+j] += depth;
								}
							}

							ptrUsed = ptrUsed->getNext();
						}
					}
					//set quality
					if (branchPosition >= 0) {
						
						for(Len_t j=0; j<(Len_t)overlapParam+branchPosition; j++) {
							contigQuality[offset+j] = 1;
						}
					}
					else {
						
						for(Len_t j=offset; j<contig.getLength(); j++) {
							contigQuality[j] = 1;
						}
					}
				}			
				
				else if (reads.getCount() == 1) {	
					
//					if (selfCircleCheckByPE(reads, pRepeats, contig))						
//						break;
					
					int branchPosition;
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					branchProcess(contig, reads, usedFlag, offset, branchPosition, cutoffScore, errorCutoff);

					
					if (branchPosition >= 0){
						extendFlag = 1;
						Len_t branchPositionInCtg = offset + overlapParam + branchPosition;
						if ((contig.getLength() > branchPositionInCtg) && (!preBranchPositions.find(branchPositionInCtg)) 
							&& (branchPositionInCtg>=contigLenOld) && (contigQuality[branchPositionInCtg] == 0)){
							preBranchPositions.append(branchPositionInCtg);	


							checkBranch(branchPosition, originLength, offset, arraySize, arrayPtr, offsetStop, contig, foundFlag, 
										foundReads, eligibleFlags, usedFlags, allInsertSizes, preBranchPositions, unreliableBranchPos, errorCutoff, false);
				
						}

					}
					else{
						//save read positions and set depth on this contig
						ListElement<ReadElement> const* ptrRead;
						ListElement<Len_t> const* ptrUsed=usedFlag.getHead();
						for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {

							if (ptrUsed->getDatum() == 1){

								extendFlag = 1;		
						
								ReadElement const& readElement = ptrRead->getDatum();
						
								readPositions[offset].append(readElement.getID());
						
								contig.appendContigPos(readElement, offset);
						
								//set depths
								Len_t depth = readElement.getDepth();
								for(Len_t j=0; j<overlapParam; j++) {
									contigDepths[offset+j] += depth;
								}
							}

							ptrUsed = ptrUsed->getNext();
						}
					}
				}

				
				
				else if ( (reads.getCount() > 1) && (reads.getCount() < maxReadsCount/10) ) {
					
					int branchPosition;
					float cutoffScore=0.8;
					Len_t errorCutoff=2;
					branchProcess(contig, reads, usedFlag, offset, branchPosition, cutoffScore, errorCutoff, false);
//					branchProcess(contig, reads, offset, branchPosition, cutoffScore, errorCutoff, false);

					
					
					if (branchPosition < 0) {	// no branch
						
//						if (selfCircleCheckByPE(reads, pRepeats, contig))						
//							break;
						
						//save read positions and set depth on this contig
						ListElement<ReadElement> const* ptrRead;
						ListElement<Len_t> const* ptrUsed=usedFlag.getHead();
						for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {

							if (ptrUsed->getDatum() == 1){

								extendFlag = 1;		

								ReadElement const& readElement = ptrRead->getDatum();
							
								readPositions[offset].append(readElement.getID());
							
								contig.appendContigPos(readElement, offset);
							
								//set depths
								Len_t depth = readElement.getDepth();
								for(Len_t j=0; j<overlapParam; j++) {
									contigDepths[offset+j] += depth;
								}
							}
							ptrUsed = ptrUsed->getNext();
						}
						
					}
					
				}


				
			}

			// check whether there is overlap between extend sequence and next contig
			if (extendFlag == 1) {
				Len_t contigLenPrevious = contigLen;
				contigLen = tStrContig.getLength();
				if (start && contigLenPrevious >= start+endNumLen) {
					
					if (nextCtgLen >= endNumLen){	
						//get seed from contig's extend part
						for (Len_t pos=contigLenPrevious-endNumLen; pos<(contigLen-endNumLen); pos++) {
						
							Number_t endNum;
							tStrContig.readTightStringFragment(pos, pos+endNumLen, endNum);

							Len_t endNumPosInNextCtg = 0;	
							if (gap.endNumHash.find(endNum) && ((endNumPosInNextCtg=gap.endNumHash[endNum])!=(Len_t)-1)) {	//seed was found in next contig	

								Len_t ctg1InFrontofCtg2 = (pos >= endNumPosInNextCtg)?1:0;
								Len_t ctg1CompareLen = (ctg1InFrontofCtg2==1)?(contigLen - pos + endNumPosInNextCtg):contigLen;
								Len_t ctg2CompareLen = (ctg1InFrontofCtg2==1)?nextCtgLen:(nextCtgLen-(endNumPosInNextCtg-pos));
								Len_t compareLen = ctg1CompareLen < ctg2CompareLen ? ctg1CompareLen : ctg2CompareLen;
								TightString const& nexttStrContig = nextContig.getTightString();
								Len_t errorSum = 0;
								Len_t errorCutoff = 2;
								Len_t posInCtg1 = (pos>=endNumPosInNextCtg)?(pos-endNumPosInNextCtg):0;
								Len_t posInCtg2 = (pos>=endNumPosInNextCtg)?0:(endNumPosInNextCtg-pos);

								// calculate mismatch number
								for (Len_t i=0; i<compareLen; ++i) {
									if (tStrContig[posInCtg1] != nexttStrContig[posInCtg2]){
										errorSum++;
										if (errorSum > errorCutoff) {
											break;
										}
									}

									posInCtg1++;
									posInCtg2++;
								}

								if (errorSum > errorCutoff) // too many mismatch
									continue;
								
							
								if (pos > contigLenOld) {
									gapContig.append(contig, contigLenOld, pos-contigLenOld);
									gapResult.length = -gap.endNumHash[endNum];
									if (((int)pos+gapResult.length) < 0) gapResult.length = -pos;
								}
								else {
									gapResult.length = (int)pos-(int)contigLenOld-(int)gap.endNumHash[endNum];
									if (((int)contigLenOld+gapResult.length) < 0) gapResult.length = -contigLenOld;
								}
								gapResult.isFilled = true;
								break;
							}
						}
					
						if (gapResult.isFilled)
							break;
					}
					else if (nextCtgLen > 0) {
						// use whole sequece of next contig as seed
						Number_t nextCtgSeq;
						nextContig.getTightString().readTightStringFragment(0, nextCtgLen, nextCtgSeq);

						for (int pos=contigLenPrevious-nextCtgLen; pos<(int)(contigLen-nextCtgLen); pos++){
							Number_t endNum;
							tStrContig.readTightStringFragment(pos, pos+nextCtgLen, endNum);
							if (endNum == nextCtgSeq){
								if (pos > (int)contigLenOld) {
									gapContig.append(contig, contigLenOld, pos-contigLenOld);
									gapResult.length = 0;
									
								}
								else {
									gapResult.length = pos-contigLenOld;
									if ((int)(contigLenOld+gapResult.length) < 0) gapResult.length = -contigLenOld;
								}
								gapResult.isFilled = true;
								break;
							}
						}

						if (gapResult.isFilled)
							break;
					}
							
						
				}
			}	
			
			if (end && (int)contigLen >= end) {  // extend sequence's length was beyond the allowed range
				break;
			}
			
			if ( (offset+overlapParam) >= contigLen ) { // no more seed can be fetched from contig's end
				break;
			}
			else{
				offset++;

				arrayPtr ++;
				arrayPtr = arrayPtr%arraySize;		
					
			}

		}



		
		if (!gapResult.isFilled) {

			Len_t branchNum = unreliableBranchPos.getCount();
			Len_t *stopBranchPositions = NULL;
			if (branchNum > 0 ){
				stopBranchPositions = new Len_t[branchNum];
				Len_t i=0;
				for (ListElement<Len_t> const* ptr=unreliableBranchPos.getHead(); ptr!=0; ptr=ptr->getNext()){
					stopBranchPositions[i] = ptr->getDatum();
					i++;
				}

				for (i=0; i<branchNum-1; i++){
					for (Len_t j=i+1; j<branchNum; j++){
						if (stopBranchPositions[i] > stopBranchPositions[j]){
							Len_t temp = stopBranchPositions[i];
							stopBranchPositions[i] = stopBranchPositions[j];
							stopBranchPositions[j] = temp;
						}
					}
				}
			}
	
			int extendLen = tStrContig.getLength() > contigLenOld ? tStrContig.getLength() - contigLenOld : 0;

			if (branchNum > 0){
				Len_t stopPosition = stopBranchPositions[0];
				for (Len_t i=1; i<branchNum; i++){
					if (stopBranchPositions[i] - stopPosition < overlapParam){
						break;
					}

					stopPosition = stopBranchPositions[i];
				}

				stopPosition = stopPosition > tStrContig.getLength() ? tStrContig.getLength() : stopPosition;
				extendLen = stopPosition > contigLenOld ? stopPosition - contigLenOld : 0;

				delete [] stopBranchPositions;
			}
			
			if (extendLen > gap.length) {
				gapContig.append(contig, contigLenOld, gap.length);
				gapResult.length = 0;
			}
			else if (extendLen > 0) {
				gapContig.append(contig, contigLenOld, extendLen);
				gapResult.length = gap.length - extendLen;
			}
			else {
				gapResult.length = gap.length;
			}
		}

		delete pRepeats;
		for (Len_t i=0; i<arraySize; i++){
			if (foundFlag[i] == 1){
				for (Len_t j=0; j<foundReads[i].getCount(); j++){
					allInsertSizes[i][j].purge();
				}
				
				delete [] allInsertSizes[i];
				allInsertSizes[i] = NULL;
				
				
			}
			foundReads[i].purge();
			eligibleFlags[i].purge();
			usedFlags[i].purge();
		}
		delete [] foundReads;
		delete [] eligibleFlags;
		delete [] foundFlag;
		delete [] usedFlags;
		delete [] allInsertSizes;

		foundReads = NULL;
		eligibleFlags = NULL;
		foundFlag = NULL;
		usedFlags = NULL;
		allInsertSizes = NULL;
	}
	
	
};

#endif /*CONTIGASSEMBLERFORFILL_HPP_*/
