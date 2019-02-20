/**************************************************
*
* ContigForFill.hpp
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

#ifndef CONTIGFORFILL_HPP_
#define CONTIGFORFILL_HPP_

//struct EndNum
//{
//	Len_t len;
//	Len_t pos;
//	
//	EndNum() : 
//		len(0), 
//		pos(0)
//	{
//	}
//};

struct ExtendInfo
{
	int leftStart;	// left start position
	int leftLen;	// left extend length
	int rightStart;
	int rightLen;
	
	ExtendInfo() : 
		leftStart(0), 
		leftLen(0), 
		rightStart(0), 
		rightLen(0)
	{	
	}
};

struct GapInfo
{
	int length;
	ExtendInfo extendInfo;
	bool isFilled;

	// key: kmers coming from next contig's beginning
	// value: positions of kmers
	HashTable< Number_t, Len_t > endNumHash;
	
	int quality;
	
	GapInfo() : 
		length(0), 
		isFilled(false), 
		quality(0)
	{
	}
	
	GapInfo(GapInfo const& gap) {
		
		length = gap.length;
		extendInfo = gap.extendInfo;
		isFilled = gap.isFilled;
		endNumHash = gap.endNumHash;
		quality = gap.quality;
	}
};

/* general contigs and gaps information */
class ContigForFill
{
protected:
	string name;
	LinkedList<TightString*> sequences;
	LinkedList<GapInfo> gaps;
	
	Len_t seqSum;
	Len_t gapSum;
	
	bool usedFlag;
	
public:
	
	ContigForFill() : 
		seqSum(0), 
		gapSum(0), 
		usedFlag(false)
	{
	}
	
	~ContigForFill() {
		
		ListElement<TightString*> const* ptr;
		for (ptr = sequences.getHead(); ptr != 0; ptr = ptr->getNext()) {
			
			delete ptr->getDatum();
		}
	}
	
	void setName(string& _name) { name = _name; }
	string const& getName() const { return name; }
	
	Len_t getSeqCount() const { return sequences.getCount(); }
	Len_t getGapCount() const { return gaps.getCount(); }
	
	Len_t getSeqSum() const { return seqSum; }
	Len_t getGapSum() const { return gapSum; }
	
	LinkedList<TightString*> const& getSequences() const { return sequences; }
//	LinkedList<GapInfo> const& getGaps() const { return gaps; }
	LinkedList<GapInfo> & getGaps() { return gaps; }	
	
	void appendSequence(TightString* seq) { 
		sequences.append(seq);
		seqSum += (*seq).getLength();
	}
	void appendGap(GapInfo& gap) {
		gaps.append(gap);
		gapSum += gap.length;
	}
	
	bool isUsed() const { return usedFlag; }
	void setUsedFlag() { usedFlag = true; }
	void clearUsedFlag() { usedFlag = false; }
	
};

#endif /*CONTIGFORFILL_HPP_*/
