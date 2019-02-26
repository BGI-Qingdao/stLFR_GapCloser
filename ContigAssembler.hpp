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

#include <fstream>

#include "Common.hpp"
#include "Contig.hpp"

class ContigAssembler
{
    public:
        static const Short_Len_t fixedOverlapMode = 1;
        static const Len_t maxReadsCount = 10000;

    protected:
        char* outfile;
        std::ofstream& fout;
        //	ofstream foutDG;
        std::ofstream foutContig;
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

        ContigAssembler(char* _outfile , std::ofstream& _fout, ReadAccessor& _readAccessor, PairInfo const& _pairInfo, Len_t _threadSum, Len_t _maxReadLength=35, Short_Len_t _overlapMode=fixedOverlapMode, Short_Len_t _overlapParam=25) : 
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

    protected:

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
};

#endif /*CONTIGASSEMBLER_HPP_*/
