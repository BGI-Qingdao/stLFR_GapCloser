/**************************************************
 *
 * Contig.hpp
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

#ifndef CONTIG_HPP_
#define CONTIG_HPP_

#include <sstream>
#include <set>
#include "Common.hpp"
#include "GlobalAccesser.hpp"
#include "base/LinkedList.hpp"
#include "ArrayBlock.hpp"
#include "readtool/ReadAccessor.hpp"
#include "TightString.hpp"

struct BarcodeInfo
{
    std::set<int> barcodes;
    // the first position that one barcode mapped.
    //  this is used to mask some area's barcodes.
    //  Key : barcode id
    //  Value : position in 0 base ;
    std::map<int , int >  m_data ;

    void UpdateBarcode(const ReadElement & read , int position )
    {
        auto bs = read.getBarodes() ;
        for( int i : bs )
        {
            UpdateBarcode(i , position );
        }
    }

    void UpdateBarcode(int barcode , int position)
    {
        if( barcode < 1 )
            return ;
        if ( m_data.find( barcode ) == m_data.end() )
        {
            m_data[barcode] = position ;
        }
        else
        {
            if( m_data.at(barcode) > position )
                m_data[barcode] = position ;
        }
        barcodes.insert(barcode);
    }
    /*
     * Erase all barcode that at [start , contig_end )
     */
    void Erase( int start )
    {
        for(const auto & pair : m_data )
        {
            if( pair.second >= start )
                barcodes.erase(pair.first);
        }
        auto itr = m_data.begin();
        while (itr != m_data.end()) 
        {
            if (barcodes.find(itr->first) == barcodes.end() )
            {
                auto toErase = itr;
                ++itr;
                m_data.erase(toErase);
            } else {
                ++itr;
            }
        }
    }
    void purge()
    {
        barcodes.clear() ;
        m_data.clear();
    }
};


/*
   define operations on contig level, including:
   1) map reads to contig.
   2) append one contig to another contig.
   3) combine two contig into one.
   4) reverse contig.
   5) output contig sequence.
   */
class Contig
{	
    public:
        static const Len_t maxReadsCount = 10000;

    protected:
        // contig sequence
        TightString* tightString;

        // subscript: position in contig
        // content: list of read IDs  started in the opsition
        ArrayBlock< LinkedList<Number_t> >* readPositions;

        /* delete depth code */
        /*
        // coverage depth
        ArrayBlock<Len_t>* depths;
        */

        /* delete quality code */
        /*
        // for one position, its quality is 1 if reads used to cover it have paired-end support
        ArrayBlock<Len_t>* quality;
        */

        // the barcodes that mapped into this contig.
        //std::set<int> barcodes ;
        BarcodeInfo barcode_info ;
        // key: read ID
        // value: list of start positions in contig
        HashTable< Number_t, LinkedList<Len_t> > readPositionsHash;
        //	map< Number_t, LinkedList<Len_t> > readPositionsHash;

        Len_t contigNum;
        Number_t endOtherContigPosReverse;
        Number_t endOtherContigPosForward;

        // block length in tight string format
        Len_t blockLenTStr;

        // block number of contig
        Len_t contigBlockLen;

        // It's the object's responsibility to free the resource if 'ownFlag' is true
        bool ownFlag;

    public:
        //250*4=1K
        Contig(Len_t length = 0, Len_t _contigNum = 0, Len_t _blockLenTStr = 250, Len_t _contigBlockLen = 1000) : 
            tightString(NULL), 
            readPositions(NULL), 
            /* delete depth code */
            /*
               depths(NULL), 
               */
            /* delete quality code */
            /*
               quality(NULL), 
               */
            contigNum(_contigNum), 
            endOtherContigPosReverse(0), 
            endOtherContigPosForward(0), 
            blockLenTStr(_blockLenTStr), 
            contigBlockLen(_contigBlockLen), 
            ownFlag(true)
    {
        allocateMemorys(length);
    }

        Contig(TightString const& tStr, Len_t _contigNum = 0, Len_t _blockLenTStr = 250, Len_t _contigBlockLen = 1000) : 
            tightString(NULL), 
            readPositions(NULL), 
            /* delete depth code */
            /*
               depths(NULL), 
               */
            /* delete quality code */
            /*
               quality(NULL), 
               */
            contigNum(_contigNum), 
            endOtherContigPosReverse(0), 
            endOtherContigPosForward(0), 
            blockLenTStr(_blockLenTStr), 
            contigBlockLen(_contigBlockLen), 
            ownFlag(true)
    {
        allocateMemorys(tStr.getLength());
        tStr.readTightStringFragment(0, tStr.getLength(), *tightString);
    }

        Contig(Contig const& leftContig, Contig const& rightContig, Len_t _contigNum = 0, Len_t _blockLenTStr = 250, Len_t _contigBlockLen = 1000) : 
            tightString(NULL), 
            readPositions(NULL), 
            /* delete depth code */
            /*
               depths(NULL), 
               */
            /* delete quality code */
            /*
               quality(NULL), 
               */
            contigNum(_contigNum), 
            endOtherContigPosReverse(0), 
            endOtherContigPosForward(0), 
            blockLenTStr(_blockLenTStr), 
            contigBlockLen(_contigBlockLen), 
            ownFlag(true)
    {
        combine(leftContig, rightContig);
    }

        Contig(Contig const& contig) {

            tightString = contig.tightString;
            readPositions = contig.readPositions;
            /* delete depth code */
            /*
               depths = contig.depths;
               */
            /* delete quality code */
            /*
               quality = contig.quality;
               */
            readPositionsHash = contig.readPositionsHash;
            contigNum = contig.contigNum;
            endOtherContigPosReverse = contig.endOtherContigPosReverse;
            endOtherContigPosForward = contig.endOtherContigPosForward;
            blockLenTStr = contig.blockLenTStr;
            contigBlockLen = contig.contigBlockLen;
            ownFlag = contig.ownFlag;
            barcode_info = contig.barcode_info ;
            const_cast<Contig&>(contig).ownFlag = false;
        }

        Contig& operator= (Contig const& contig) {

            if (&contig!=this) {
                purge();
                tightString = contig.tightString;
                readPositions = contig.readPositions;
                /* delete depth code */
                /*
                   depths = contig.depths;
                   */
                /* delete quality code */
                /*
                   quality = contig.quality;
                   */
                readPositionsHash = contig.readPositionsHash;
                contigNum = contig.contigNum;
                endOtherContigPosReverse = contig.endOtherContigPosReverse;
                endOtherContigPosForward = contig.endOtherContigPosForward;
                blockLenTStr = contig.blockLenTStr;
                contigBlockLen = contig.contigBlockLen;
                ownFlag = contig.ownFlag;
                barcode_info = contig.barcode_info;
                const_cast<Contig&>(contig).ownFlag = false;
            }
            return *this;
        }

        ~Contig() {

            purge();
        }

        void allocateMemorys(Len_t length) {

            tightString = new TightString(length, blockLenTStr);

            Len_t blockLen;
            if (length%contigBlockLen) {
                blockLen = ((length/contigBlockLen)+1)*contigBlockLen;
            }
            else {
                blockLen = length;
            }

            readPositions = new ArrayBlock< LinkedList<Number_t> >(contigBlockLen, blockLen);
            /* delete depth code */
            /*
               depths = new ArrayBlock<Len_t>(contigBlockLen, blockLen);
               (*depths).clear(0, blockLen);
               */
            /* delete quality code */
            /*
               quality = new ArrayBlock<Len_t>(contigBlockLen, blockLen);
               (*quality).clear(0, blockLen);
               */
        }

        void purge() {

            if (ownFlag) {
                if (tightString) delete tightString;
                if (readPositions) delete readPositions;
                /* delete depth code */
                /*
                   if (depths) delete depths;
                   */
                /* delete quality code */
                /*
                   if (quality) delete quality;
                   */
            }

            tightString = NULL;
            readPositions = NULL;
            /* delete depth code */
            /*
               depths = NULL;
               */
            /* delete quality code */
            /*
               quality = NULL;
               */
            readPositionsHash.purge();
            barcode_info.purge();
            //		readPositionsHash.clear();
        }

        operator Number_t() const { return getLength(); }

        Len_t getLength() const { return tightString->getLength(); }
        void setLength(Len_t length) { tightString->setLength(length); }

        TightString const& getTightString() const { return *tightString; }
        TightString& getTightString() { return *tightString; }
        void setTightString(TightString* _tightString) { tightString = _tightString; }

        ArrayBlock< LinkedList<Number_t> > const& getReadPositions() const { return *readPositions; }
        ArrayBlock< LinkedList<Number_t> >& getReadPositions() { return *readPositions; }
        void setReadPositions(ArrayBlock< LinkedList<Number_t> >* _readPositions) { readPositions = _readPositions; }

        void getReadPositionsReverse(ArrayBlock< LinkedList<Number_t> >& readPositionsReverse) {

            reverseReadPosition(*readPositions, readPositionsReverse, 0, getLength());
        }

        /* delete depth code */
        /*
           ArrayBlock<Len_t> const& getDepths() const { return *depths; }
           ArrayBlock<Len_t>& getDepths() { return *depths; }
           void setDepths(ArrayBlock<Len_t>* _depths) { depths = _depths; }
           */
        /* delete quality code */
        /*
           ArrayBlock<Len_t> const& getQuality() const { return *quality; }
           ArrayBlock<Len_t>& getQuality() { return *quality; }
           void setQuality(ArrayBlock<Len_t>* _quality) { quality = _quality; }
           */
        const std::set<int> & getBarodes() const { return barcode_info.barcodes ; }
        BarcodeInfo & getBarodeInfo() { return barcode_info ; }

        const BarcodeInfo & getBarodeInfo() const { return barcode_info ; }

        HashTable< Number_t, LinkedList<Len_t> >& getReadPositionsHash() { return readPositionsHash; }
        //	map< Number_t, LinkedList<Len_t> >& getReadPositionsHash() { return readPositionsHash; }

        Len_t getContigNum() const { return contigNum; }
        void setContigNum(Len_t _contigNum) { contigNum = _contigNum; }

        Number_t getEndOtherContigPosReverse() const { return endOtherContigPosReverse; }
        void setEndOtherContigPosReverse(Number_t _endOtherContigPosReverse) { endOtherContigPosReverse = _endOtherContigPosReverse; }

        Number_t getEndOtherContigPosForward() const { return endOtherContigPosForward; }
        void setEndOtherContigPosForward(Number_t _endOtherContigPosForward) { endOtherContigPosForward = _endOtherContigPosForward; }

        /* output coverage and quality information */
        void outputContig(std::ofstream& fout) {

            outputContigFastAFormat(fout);

            Len_t length = getLength();
            Len_t pos;
            for (pos=0; pos<length; ++pos) {

                if ((*readPositions)[pos].getCount()>0) {

                    fout << pos << ":";

                    ListElement<Number_t> const* ptr;
                    for (ptr = (*readPositions)[pos].getHead(); ptr != 0; ptr = ptr->getNext()) {

                        ReadElement const& readElement =
                            GlobalAccesser::the_read_accessor->getRead(ptr->getDatum());
                        Len_t depth = readElement.getDepth();
                        for (Len_t j=0; j<depth; j++) {
                            Number_t iRead = readElement.getSerialNums()[j];
                            //						iRead >>= 1;
                            fout << iRead << ",";
                        }
                    }

                    fout << " ";
                }
            }
            fout << std::endl;



            ArrayBlock< LinkedList<Number_t> > readPositionsReverse(contigBlockLen, length);
            reverseReadPosition(*readPositions, readPositionsReverse, 0, length);
            for (pos=0; pos<length; ++pos) {

                if (readPositionsReverse[pos].getCount()>0) {

                    fout << pos << ":";

                    ListElement<Number_t> const* ptr;
                    for (ptr = readPositionsReverse[pos].getHead(); ptr != 0; ptr = ptr->getNext()) {

                        ReadElement const& readElement =
                            GlobalAccesser::the_read_accessor->getRead(ptr->getDatum());
                        Len_t depth = readElement.getDepth();
                        for (Len_t j=0; j<depth; j++) {
                            Number_t iRead = readElement.getSerialNums()[j];
                            //						iRead >>= 1;
                            fout << iRead << ",";
                        }
                    }

                    fout << " ";
                }
            }
            fout << std::endl;



            /* delete depth code */
            /*
               for (pos=0; pos<length; ++pos) {
               fout << (*depths)[pos] << ',';
               }
               */
            fout << std::endl;

        }

        void outputContigFastAFormat(std::ofstream& fout, int reverseExtendLen, int forwardExtendLen) {

            char* seq = (*tightString).readTightString();
            fout << ">contigs_" << contigNum << "_length_" << getLength() 
                << "_reverse_" << reverseExtendLen 
                << "_forward_" << forwardExtendLen 
                << std::endl;
            fout << seq << std::endl;
            delete [] seq;
        }

        void outputContigFastAFormat(std::ofstream& fout) {

            char* seq = (*tightString).readTightString();
            fout << ">contigs_" << contigNum << "_length_" << getLength() 
                << "_reverse_" << (endOtherContigPosReverse>>33) << "_" << ((endOtherContigPosReverse & 0x00000001FFFFFFFFLLU) >> 2) << "_" << ((endOtherContigPosReverse&2) >> 1) 
                << "_forward_" << (endOtherContigPosForward>>33) << "_" << ((endOtherContigPosForward & 0x00000001FFFFFFFFLLU) >> 2) << "_" << ((endOtherContigPosForward&2) >> 1) 
                << std::endl;
            fout << seq << std::endl;
            delete [] seq;
        }

        void outputDGFormat(std::ofstream& fout) {

            Len_t length = getLength();
            std::ostringstream outputReadsID, outputPos;
            Len_t readsIDSum = 0;
            for (Len_t pos=0; pos<length; ++pos) {

                ListElement<Number_t> const* ptr;
                for (ptr = (*readPositions)[pos].getHead(); ptr != 0; ptr = ptr->getNext()) {

                    ReadElement const& readElement =
                            GlobalAccesser::the_read_accessor->getRead(ptr->getDatum());
                    Len_t depth = readElement.getDepth();
                    for (Len_t j=0; j<depth; j++) {
                        Number_t iRead = readElement.getSerialNums()[j];
                        outputReadsID << iRead << "\t";
                        outputPos << pos << "\t";
                        readsIDSum++;
                    }
                }
            }

            fout << ">contig" << contigNum << "\t" << length << "\t" << readsIDSum << std::endl;
            fout << outputReadsID.str() << std::endl;
            fout << outputPos.str() << std::endl;
        }

        // function:	 mapReads
        // description: map reads to specific region of contig and record aligned information
        // input:		 start -- start position in contig
        //			 len -- length of region for mapping
        //			 overlapParam -- length of seed used to find reads
        //			 _mismatch: mismatch number allowed beyond seed
        // output:	 none
        void mapReads(int start = 0, Len_t len = 0, Short_Len_t overlapParam = 25, Short_Len_t _mismatch = 2) {

            if (len<overlapParam) return;
            TightString& tStrContig = *(this->tightString);
            ArrayBlock< LinkedList<Number_t> >& readPositions = *(this->readPositions);
            /* delete depth code */
            /*
               ArrayBlock<Len_t>& depths = *(this->depths);
               */
            if (start < 0)
                start = 0;
            if ( (len == 0) || (len > getLength()-start) )
                len = getLength()-start;
            Len_t end = start+len;
            for (Len_t i=start; i<=end-overlapParam; ++i) {

                //get overlap part by param
                TightString tStrOverlap(overlapParam);
                tStrContig.readTightStringFragment(i, i+overlapParam, tStrOverlap);

                LinkedList<ReadElement> reads;
                GlobalAccesser::the_read_accessor->getReadsBeginWith
                    (tStrOverlap,reads,false,ReadAccessor::inAll);

                if (reads.getCount() > maxReadsCount)
                    continue;

                Short_Len_t mismatch = _mismatch;
                if (reads.getCount() > maxReadsCount/10) {
                    mismatch = 0;
                }
                else if (reads.getCount() > maxReadsCount/100) {
                    mismatch = 1;
                }

                //save read positions and set depth on this contig
                ListElement<ReadElement> const* ptrRead;
                for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) {

                    ReadElement const& readElement = ptrRead->getDatum();
                    Len_t readLength = readElement.getLen();

                    //check alignment
                    TightString tStrRead(readLength);
                    readElement.getSequence(tStrRead);

                    Len_t errorSum = 0;
                    int readLengthRemain = readLength - overlapParam;
                    int contigLengthRemain = end - overlapParam - i;
                    if (readLengthRemain > contigLengthRemain+mismatch) {
                        continue;
                    }
                    Len_t lengthRemain = readLengthRemain < contigLengthRemain ? readLengthRemain : contigLengthRemain;
                    for(Len_t j=0; j<lengthRemain; j++) {

                        if (tStrRead[j+overlapParam] != tStrContig[i+j+overlapParam]) {
                            errorSum++;
                            if (errorSum>mismatch)
                                break;
                        }
                    }

                    if (errorSum<=mismatch) {

                        readPositions[i].append(readElement.getID());

                        appendContigPos(readElement, i);
                        appendBarcodes(readElement,i);
                        /* delete depth code */
                        /*
                        //set depths
                        Len_t depth = readElement.getDepth();
                        for(Len_t j=0; j<lengthRemain+overlapParam; j++) {
                        depths[i+j] += depth;
                        }
                        */
                    }
                }

            }
        }

        // function: 	 getContigPos
        // description: find read's start positions in contig
        // input:		 id -- read ID
        // output:	 pos -- list of positions
        void getContigPos(Number_t id, LinkedList<Len_t>& pos) const {

            if (readPositionsHash.find(id)) {
                //		if (readPositionsHash.count(id)) {
                pos = readPositionsHash[id];
            }
            }

            // function: 	 getContigPos
            // description: find read's start positions in contig
            // input:		 readElement -- read element
            // output:	 pos -- list of positions
            void getContigPos(ReadElement const& readElement, LinkedList<Len_t>& pos) {

                if (readPositionsHash.find(readElement.getID())) {
                    //		if (readPositionsHash.count(readElement.getID())) {
                    pos = readPositionsHash[readElement.getID()];
                }
                }

                void appendBarcodes(ReadElement const& readElement, int pos )
                {
                    barcode_info.UpdateBarcode(readElement , pos );
                    //SetAdd(barcodes,readElement.getBarodes());
                }

                void appendBarcodes()
                {
                    appendBarcodes(0, getLength());
                }

                void appendBarcodes( Len_t start, Len_t len)
                {
                    ArrayBlock< LinkedList<Number_t> >& readPositions = *(this->readPositions);
                    Len_t end = start+len;
                    for (Len_t i=start; i<end; ++i) {
                        ListElement<Number_t> const* ptr;
                        for (ptr = readPositions[i].getHead(); ptr != 0; ptr = ptr->getNext()) {
                            appendBarcodes(
                                    GlobalAccesser::the_read_accessor->getRead(ptr->getDatum()) 
                                    , i );
                        }
                    }
                }

                // function: 	 appendContigPos
                // description: append read's start position in contig to list
                // input:		 readElement -- read element
                // 			 pos -- start position
                // output:	 none
                void appendContigPos(ReadElement const& readElement, Len_t pos) {

                    readPositionsHash[readElement.getID()].append(pos);
                }

                // function: 	 appendContigPos
                // description: append read's start position in whole contig to list
                // input:		 none
                // output:	 none
                void appendContigPos() {

                    appendContigPos(0, getLength());
                }

                // function: 	 appendContigPos
                // description: append read's start position in specific contig region to list
                // input:		 start -- start position in contig
                //			 len -- region length
                // output:	 none
                void appendContigPos(Len_t start, Len_t len) {

                    ArrayBlock< LinkedList<Number_t> >& readPositions = *(this->readPositions);
                    Len_t end = start+len;
                    for (Len_t i=start; i<end; ++i) {

                        ListElement<Number_t> const* ptr;
                        for (ptr = readPositions[i].getHead(); ptr != 0; ptr = ptr->getNext()) {

                            readPositionsHash[ptr->getDatum()].append(i);
                        }
                    }
                }

                // function: 	 append
                // description: append another whole contig to current contig
                // input:		 contigAppend -- another contig in Contig format
                // output:	 none
                void append(Contig const& contigAppend) {

                    append(contigAppend, 0, contigAppend.getLength());
                }

                // function: 	 append
                // description: append another contig's specific region to current contig
                // input:		 contigAppend -- another contig in Contig format
                //			 start -- start position of another contig
                //			 appendLen -- append length of another contig
                // output:	 none
                void append(Contig const& contigAppend, Len_t start, Len_t appendLen) {

                    TightString tStrAppend(appendLen);
                    contigAppend.getTightString().readTightStringFragment(start, start+appendLen, tStrAppend);
                    Len_t startOrigin = 0;

                    if (!tightString) {
                        allocateMemorys(appendLen);
                        tStrAppend.readTightStringFragment(0, appendLen, *tightString);
                    }
                    else {
                        startOrigin = (*tightString).getLength();
                        (*tightString).append(tStrAppend, appendLen);
                        Len_t newLength = (*tightString).getLength();

                        //			(*readPositions).increase(newLength);
                        (*readPositions).setLength(newLength);
                        /* delete depth code */
                        /*
                           (*depths).increase(newLength);
                           */

                        /* delete quality code */
                        /*
                           (*quality).increase(newLength);
                           */
                    }

                    (*readPositions).append(contigAppend.getReadPositions(), appendLen, start, startOrigin);
                    /* delete depth code */
                    /*
                       (*depths).append(contigAppend.getDepths(), appendLen, start, startOrigin);
                       */

                    /* delete quality code */
                    /*
                       (*quality).append(contigAppend.getQuality(), appendLen, start, startOrigin);
                       */
                    appendContigPos(startOrigin, appendLen);
                    appendBarcodes(startOrigin, appendLen);
                }

                // function: 	 append
                // description: append another contig's specific region to current contig
                // input:		 contigAppend -- another contig in TightString format
                //			 appendLen -- append length of another contig
                // output:	 none
                void append(TightString const& tStrAppend, Len_t appendLen) {

                    if (tightString) {

                        (*tightString).append(tStrAppend, appendLen);
                        Len_t newLength = (*tightString).getLength();

                        //			(*readPositions).increase(newLength);
                        (*readPositions).setLength(newLength);
                        /* delete depth code */
                        /*
                           (*depths).increase(newLength);
                           */
                        /* delete quality code */
                        /*
                           (*quality).increase(newLength);
                           */
                    }
                    else {

                        allocateMemorys(appendLen);
                        tStrAppend.readTightStringFragment(0, appendLen, *tightString);
                    }
                }

                // function: 	 combine
                // description: combine two contigs in Contig format into one and replace current contig
                // input:		 leftContig -- contig at the beginning of combined contig
                //			 rightContig --  contig at the end of combined contig
                // output:	 none
                void combine(Contig const& leftContig, Contig const& rightContig) {

                    Len_t leftLength = leftContig.getLength();
                    Len_t rightLength = rightContig.getLength();
                    Len_t length = leftLength + rightLength;

                    TightString& tStrContig = *new TightString(leftContig.getTightString(), rightContig.getTightString(), blockLenTStr);

                    Len_t combineBlockLen;
                    if (length%contigBlockLen) {
                        combineBlockLen = ((length/contigBlockLen)+1)*contigBlockLen;
                    }
                    else {
                        combineBlockLen = length;
                    }

                    ArrayBlock< LinkedList<Number_t> >& readPositions = *new ArrayBlock< LinkedList<Number_t> >(contigBlockLen, combineBlockLen);
                    readPositions.combine(leftContig.getReadPositions(), leftLength, rightContig.getReadPositions(), rightLength);

                    /* delete depth code */
                    /*
                       ArrayBlock<Len_t>& depths = *new ArrayBlock<Len_t>(contigBlockLen, combineBlockLen);
                       depths.clear(0, combineBlockLen);
                       depths.combine(leftContig.getDepths(), leftLength, rightContig.getDepths(), rightLength);
                       */

                    /* delete quality code */
                    /*
                       ArrayBlock<Len_t>& quality = *new ArrayBlock<Len_t>(contigBlockLen, combineBlockLen);
                       quality.clear(0, combineBlockLen);
                       quality.combine(leftContig.getQuality(), leftLength, rightContig.getQuality(), rightLength);
                       */
                    purge();
                    setTightString(&tStrContig);
                    setReadPositions(&readPositions);
                    /* delete depth code */
                    /*
                       setDepths(&depths);
                       */
                    /* delete quality code */
                    /*
                       setQuality(&quality);
                       */
                    getReadPositionsHash().initialize();
                    appendContigPos();
                    barcode_info.purge();
                    appendBarcodes() ;
                }

                // function: 	 reverse
                // description: reverse specific region of contig
                // input:		 start -- start position
                //			 end -- end position
                // output:	 contigReverse -- reversed contig region
                void reverse(Contig& contigReverse, Len_t start, Len_t end) const {

                    Len_t reverseLen = end - start;
                    TightString& tStrContigReverse = *new TightString(reverseLen, blockLenTStr);
                    (*tightString).readTightStringFragmentAtReverse(start, end, tStrContigReverse);

                    Len_t reverseBlockLen;
                    if (reverseLen%contigBlockLen) {
                        reverseBlockLen = ((reverseLen/contigBlockLen)+1)*contigBlockLen;
                    }
                    else {
                        reverseBlockLen = reverseLen;
                    }

                    ArrayBlock< LinkedList<Number_t> >& readPositionsReverse = *new ArrayBlock< LinkedList<Number_t> >(contigBlockLen, reverseBlockLen);
                    reverseReadPosition(*readPositions, readPositionsReverse, start, end);

                    ArrayBlock<Len_t>& depthsReverse = *new ArrayBlock<Len_t>(contigBlockLen, reverseBlockLen);
                    depthsReverse.clear(0, reverseBlockLen);
                    /* delete depth code */
                    /*
                       (*depths).reverse(depthsReverse, start, end);
                       */
                    ArrayBlock<Len_t>& qualityReverse = *new ArrayBlock<Len_t>(contigBlockLen, reverseBlockLen);
                    qualityReverse.clear(0, reverseBlockLen);
                    /* delete quality code */
                    /*
                       (*quality).reverse(qualityReverse, start, end);
                       */
                    contigReverse.purge();
                    contigReverse.setTightString(&tStrContigReverse);
                    contigReverse.setReadPositions(&readPositionsReverse);
                    /* delete depth code */
                    /*
                       contigReverse.setDepths(&depthsReverse);
                       */
                    /* delete quality code */
                    /*
                       contigReverse.setQuality(&qualityReverse);
                       */
                    contigReverse.getReadPositionsHash().initialize();
                    contigReverse.appendContigPos();
                    contigReverse.barcode_info.purge();
                    contigReverse.appendBarcodes() ;

                }

                // function: 	 reverseReadPosition
                // description: reverse record of aligned reads in specific contig region
                // input:		 readPositions -- record to reverse
                //			 start -- start position of region
                //			 end -- end position of region
                // output:	 readPositionsReverse -- reversed record
                void reverseReadPosition(
                        ArrayBlock< LinkedList<Number_t> >& readPositions, 
                        ArrayBlock< LinkedList<Number_t> >& readPositionsReverse, 
                        Len_t start, 
                        Len_t end
                        ) const {

                    for (Len_t i=start; i<end; i++) {

                        ListElement<Number_t> const* ptr;
                        for (ptr = readPositions[i].getHead(); ptr != 0; ptr = ptr->getNext()) {

                            ReadElement const& readElement =
                                GlobalAccesser::the_read_accessor->getRead(ptr->getDatum());
                            int pos = end-i-readElement.getLen();
                            if (pos >= 0) {
                                readPositionsReverse[pos].append(readElement.getID());
                            }
                            else {
                                continue;
                            }
                        }
                    }
                }

            };


#endif /*CONTIG_HPP_*/
