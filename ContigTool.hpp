#ifndef CONTIGTOOL_HPP__
#define CONTIGTOOL_HPP__

#include <tuple>
#include "Contig.hpp"
#include "Utils.hpp"
#include "ContigForFill.hpp"
#include "GlobalAccesser.hpp"
#include "readtool/PairInfo.hpp"

struct  ContigTool
{
    public:

        static std::vector<ReadElement> find_all_reads_start_with(
                const Number_t & kmer
                , const Contig & contig 
                , int position /* in 1 base */
                )
        {
            std::vector<ReadElement> ret ;
            LinkedList<ReadElement> reads ;
            reads.purge();
            GlobalAccesser::the_read_accessor->getReadsBeginWith(kmer,reads,true,ReadAccessor::inAll);
            ListElement<ReadElement> const* ptrRead;

            for (ptrRead = reads.getHead(); ptrRead != 0; ptrRead = ptrRead->getNext()) 
            {
                ReadElement const& readElement = ptrRead->getDatum();
                if ((int)readElement.getDepth()
                        > Threshold::max_reads_depth ) 
                    continue ;
                if( ! IsReadMatchContig( readElement , contig , position) )
                    continue ;
                ret.push_back( readElement ) ;
            }
            return ret ;
        }

        static bool IsReadMatchContig( const ReadElement & readElement 
                , const Contig & contig 
                , int position /* in 1 base */
                )
        {
            TightString  read_str(readElement.getLen());
            readElement.getSequence(read_str) ;
            int match_check_num = 0 ;
            int unmatch = CheckUnmatchNum( 
                    contig.getTightString(),
                    contig.getLength() ,
                    read_str ,
                    readElement.getLen() ,
                    position - 1 , /*1base->0base*/
                    match_check_num);

            return  unmatch <= Threshold::max_error_count ;
        };

        static std::tuple<bool,int,int>
            IsGapFinish(const Contig & contig
                , const Contig & nextContig
                , const GapInfo & gap
                , int contigLenPrevious
                )
        {
            int contigLen = contig.getLength() ;
            int nextCtgLen = nextContig.getLength() ;
            int endNumLen = Threshold::the_k ;
            const TightString & tStrContig = contig.getTightString();

            for (int pos=contigLenPrevious-endNumLen;
                    pos<(contigLen-endNumLen); pos++)
            {
                Number_t endNum;
                tStrContig.readTightStringFragment
                    (pos, pos+endNumLen, endNum);
                int endNumPosInNextCtg = 0;	
                if (gap.endNumHash.find(endNum)
                    && ((gap.endNumHash[endNum])!=(Len_t)-1) )
                {   //seed was found in next contig	
                    endNumPosInNextCtg = gap.endNumHash[endNum] ;
                    Len_t ctg1InFrontofCtg2 = 
                        (pos >= endNumPosInNextCtg)?1:0;
                    Len_t ctg1CompareLen = 
                        (ctg1InFrontofCtg2==1)?
                        (contigLen - pos + endNumPosInNextCtg):contigLen;
                    Len_t ctg2CompareLen = 
                        (ctg1InFrontofCtg2==1)?
                        nextCtgLen:(nextCtgLen-(endNumPosInNextCtg-pos));
                    Len_t compareLen = 
                        ctg1CompareLen < ctg2CompareLen ? 
                        ctg1CompareLen : ctg2CompareLen;
                    TightString const& nexttStrContig = 
                        nextContig.getTightString();
                    Len_t errorSum = 0;
                    Len_t posInCtg1 = (pos>=endNumPosInNextCtg)
                        ?(pos-endNumPosInNextCtg):0;
                    Len_t posInCtg2 = (pos>=endNumPosInNextCtg)?
                        0:(endNumPosInNextCtg-pos);

                    // calculate mismatch number
                    for (int i=0; i<(int)compareLen; ++i) {
                        if (tStrContig[posInCtg1] 
                                != nexttStrContig[posInCtg2]){
                            errorSum++;
                            if ((int)errorSum > Threshold::max_error_count)
                            {
                                break;
                            }
                        }
                        posInCtg1++;
                        posInCtg2++;
                    }

                    if ((int)errorSum > Threshold::max_error_count) 
                        continue;
                    return std::make_tuple( true, pos, endNumPosInNextCtg);
                }
            }
            return std::make_tuple(false,-1,-1);
        }

        static bool
            IsEligibleBarcodeCheckRead_NOTPE(
                    ReadElement const& readElement
                    , const Contig& prev_contig
                    , const Contig& next_contig
                    , int pe_pos
                    )
            {
                auto rbarcode = readElement.getBarodes() ;
                if( rbarcode.empty() )
                    return false ;
                const auto & prev_barcode = prev_contig.getBarodes();
                if( prev_barcode.empty() )
                    return false ;
                auto with_prev = SetIntersection(rbarcode,prev_barcode);
                if( with_prev.empty()) 
                    return false ;
                const auto &  next_barcode = next_contig.getBarodes() ;
                if( next_barcode.empty()) 
                    return false ;
                auto  with_next = SetIntersection( rbarcode, next_barcode);
                if( with_next.empty() )
                    return false ;
                auto both = SetIntersection( with_prev, with_next );
                if( both.empty() )
                    return false ;
                // to check not pe barcode
                for( const auto a_barcode : both )
                {
                    const auto & bi = prev_contig.getBarodeInfo() ;
                    const auto & bps = bi.m_barcode_pos.at(a_barcode);
                    if( bps.size() <= 1 )
                        continue ;
                    for( const auto pos : bps )
                    {
                        if( pos != pe_pos )
                            return true ;
                    }
                }
                return false ;
            }

        static bool
            IsEligibleBarcodeCheckRead( ReadElement const& readElement
                    , const Contig& prev_contig 
                    , const Contig& next_contig 
                    )
            {
                auto rbarcode = readElement.getBarodes() ;
                if( rbarcode.empty() )
                    return false ;
                const auto & prev_barcode = prev_contig.getBarodes();
                if( prev_barcode.empty() )
                    return false ;
                auto with_prev = SetIntersection(rbarcode,prev_barcode);
                if( with_prev.empty()) 
                    return false ;
                const auto &  next_barcode = next_contig.getBarodes() ;
                if( next_barcode.empty()) 
                    return false ;
                auto  with_next = SetIntersection( rbarcode, next_barcode);
                if( with_next.empty() )
                    return false ;
                auto both = SetIntersection( with_prev, with_next );
                if( both.empty() )
                    return false ;
                return true  ;
            }

        static bool
            IsEligibleONTCheckRead( ReadElement const& /*readElement*/
                    , const Contig& /*prev_contig */
                    , const Contig& /*next_contig */
                    )
            {
                return false;
            }
        static bool
            IsEligiblePECheckRead( ReadElement const& readElement
                    , const Contig& contig
                    , int read_pos 
                    , int & pe_pos ) 
            {
                std::vector<std::vector<int>> contigPos;
                Len_t arrayLen = GlobalAccesser::the_pair_info->getArrayLen();
                contigPos.resize(arrayLen);
                LinkedList<Number_t>* ids = new LinkedList<Number_t>[arrayLen];
                GlobalAccesser::the_read_accessor->getReadIdsByPair(readElement, ids);
                for (Len_t i=0; i<arrayLen; i++) {
                    getContigPos(ids[i], contigPos[i], contig);
                }
                delete [] ids;
                Len_t PECount = 0;		
                for (Len_t i=0; i<arrayLen; i++) {	//process rank by rank
                    PairInfoElement const& pairInfoElement = GlobalAccesser::the_pair_info->getArray()[i];
                    for (int pos : contigPos[i] ) {

                        //check the actual insert size
                        Len_t actualInsertSize = read_pos - pos + readElement.getLen();
                        if ( 
                                ( actualInsertSize
                                  >= 
                                  (pairInfoElement.insertSize-pairInfoElement.variance) 
                                ) 
                                && 
                                ( actualInsertSize 
                                  <= 
                                  (pairInfoElement.insertSize+pairInfoElement.variance) 
                                ) 
                           ) 
                        {
                            pe_pos = pos ;
                            PECount++;		
                        }
                    }
                }
                return PECount > 0 ;
            }

        static void AddReadIntoContig(
                Contig & contig  ,
                const ReadElement & readElement ,
                int offset /* 0base*/
                )
        {

            auto & readPositions = contig.getReadPositions() ;
            readPositions[offset].append(readElement.getID());
            contig.appendContigPos(readElement, offset);
            contig.appendBarcodes(readElement,offset);
        }

    private :

        static void getContigPos(
                LinkedList<Number_t>& ids
                , std::vector<int>& contigPos
                , const Contig& contig)  
        {
            ListElement<Number_t> const* ptrId;
            for (ptrId = ids.getHead(); ptrId != 0; ptrId = ptrId->getNext()) {
                Number_t id = ptrId->getDatum();
                LinkedList<Len_t> pos;
                contig.getContigPos(id, pos);
                ListElement<Len_t> const* ptr;
                for (ptr = pos.getHead(); ptr != 0; ptr = ptr->getNext()) {
                    contigPos.push_back(ptr->getDatum());
                }
            }
        }

        static int CheckUnmatchNum(
                const TightString & ref 
                , int ref_len 
                , const TightString & read
                , int read_len 
                , int read2ref_1th_pos /* 0 base */
                , int & checked_num
                )
        {
            int ret = 0 ;

            int contig_start = read2ref_1th_pos ;
            int contig_end = contig_start + read_len -1 ;
            if( contig_end > ref_len  )
                contig_end = ref_len ;
            checked_num = contig_end - contig_start +1 ;
            for (int m= 0; m< checked_num ; m++)
                if ( ref[contig_start + m ] != read[m] )
                    ret ++ ;
            return ret ;
        }

};

#endif //CONTIGTOOL_HPP__ 
