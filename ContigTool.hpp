#ifndef CONTIGTOOL_HPP__
#define CONTIGTOOL_HPP__

#include "Contig.hpp"
#include "Utils.hpp"
#include "GlobalAccesser.hpp"
#include "readtool/PairInfo.hpp"

struct  ContigTool
{
    public:

        static std::vector<ReadElement> find_all_reads_start_with(
                const Number_t & kmer
                , const Contig & contig 
                , int position
                , int max_reads_depth 
                , int max_error_count
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
                if ((int)readElement.getDepth() >= max_reads_depth ) 
                {
                    continue ;
                }
                TightString  read_str ;
                readElement.getSequence(read_str) ;
                int match_check_num = 0 ;
                int unmatch =  CheckUnmatchNum( 
                        contig.getTightString(),
                        contig.getLength() ,
                        read_str ,
                        readElement.getLen() ,
                        position - 1 , /*1base->0base*/
                        match_check_num);

                if( unmatch > max_error_count )
                {
                    continue ;
                }
                ret.push_back( readElement ) ;
            }
            return ret ;
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
                    , int read_pos ) 
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
                            PECount++;		
                        }
                    }
                }
                return PECount > 0 ;
            }
    private :

        static void getContigPos(LinkedList<Number_t>& ids, std::vector<int>& contigPos, const Contig& contig) {
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
                , int read2ref_1th_pos
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
