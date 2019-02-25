#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>

#include "Common.hpp"
#include "TightString.hpp"
#include "readtool/ReadElement.hpp"
#include "readtool/ReadAccessor.hpp"
#include "Contig.hpp"
#include "ConsensusConfig.hpp"
#include "ContigForFill.hpp"

struct Kmer2Reads
{
    Number_t the_kmer ;
    std::vector<ReadElement> reads;
    int relative_num ;
};

struct PositionInfoKeeper
{
    private:

        std::map<Number_t , Kmer2Reads> all_kmer2reads;

        int reads_num ;

    public :
        void Init()
        {
            reads_num = 0 ;
        }

        int ReadsNum() const {
            return reads_num ;
        }

        bool KmerExist( const Number_t & kmer)
        {
            auto itr = all_kmer2reads.find(kmer) ;
            if( itr == all_kmer2reads.end() )
                return false ;
            return true ;
        }
        void AddKmer(const Number_t & kmer,
                const std::vector<ReadElement> & reads ,
                int relative_num )
        {
            auto & new_one = all_kmer2reads[kmer] ;
            new_one.the_kmer = kmer ;
            new_one.relative_num = relative_num ;
            new_one.reads = reads ;
            reads_num += new_one.reads.size() ;
        }
};

struct ConsensusResult
{
    static int max_allowed_conflict ;

    std::string       consensus_str;
    std::vector<bool> conflict_flags;
    int conflict_num() const {
        int i = 0 ;
        for( bool x : conflict_flags )
            if( !x )
                i++ ;

        return i ;
    }

    bool is_consensus_done() const 
    {
        return ! (conflict_num() > max_allowed_conflict );
    }
};

int ConsensusResult::max_allowed_conflict = 2 ;

struct ConsensusMatrix
{
    std::vector<int[4]> depth_matrix;

    ConsensusResult GenConsensusResult() const 
    {
        ConsensusResult ret ;

        return ret ;
    }
};

struct ReadMatrix
{
    static int  max_reads_count ;
    static int  min_reads_count ;

    static int min_sub_reads_count;

    private:
    //  position --> reads
    std::map< int  , PositionInfoKeeper> m_raw_reads;

    public: 

    static ReadMatrix Merge( const ReadMatrix & /*left*/ , const ReadMatrix & /*right*/ )
    {
        ReadMatrix ret ;
        //TODO 
        return ret ;
    }

    void AddKmer(const Number_t & kmer,
            const std::vector<ReadElement> & reads , 
            int relative_num , int pos )
    {
        m_raw_reads[pos].AddKmer(kmer,reads,relative_num);
    }
    void tryInitPos(int pos)
    {
        if( m_raw_reads.find(pos) == m_raw_reads.end())
            m_raw_reads[pos].Init() ;
    }

    bool checkKmerInPos(const Number_t & kmer , int pos )
    {
        if( m_raw_reads.find(pos) == m_raw_reads.end())
            return false ;
        return m_raw_reads[pos].KmerExist(kmer);
    }

    ConsensusArea m_area;

    ConsensusMatrix GenConsensusMatrix(  )
    {
        ConsensusMatrix ret ;
        //TODO
        return ret ;
    }
    ReadMatrix GetSubMatrixByPECheck(const Contig & prev_contig)
    {
        ReadMatrix ret ;
        //TODO
        return ret ;
    }

    ReadMatrix GetSubMatrixByBarcodeCheck(
            const Contig & prev_contig,
            const Contig & next_contig
            )
    {
        ReadMatrix ret ;
        //TODO
        return ret ;
    }

    ReadMatrix GenSubMatrixByGap(const Contig & prev_contig ,
            const Contig & next_contig ,
            const GapInfo & gap )
    {

        if( ReadsNum() < min_sub_reads_count )
            return *this ;
        if( gap.is_gap_big() ) {;} else {;}
        //{
        auto sub1 = GetSubMatrixByPECheck( prev_contig) ;
        auto sub2 = GetSubMatrixByBarcodeCheck(prev_contig , next_contig );
        auto sub3 = Merge( sub1 , sub2 );
        if( sub3.ReadsNum() < min_sub_reads_count )
            return *this ;
        else 
            return sub3 ;
        //}
        //else
        //{

        //}
        return *this ;
    }

    ReadMatrix & operator = ( const ReadMatrix & other )
    {
        if( this != &other )
        {
            m_area = other.m_area ;
            m_raw_reads = other.m_raw_reads ;
        }
        return *this ;
    }

    int ReadsNum() const 
    {
        int ret = 0 ;
        for( const auto & pair : m_raw_reads)
        {
            ret += pair.second.ReadsNum() ;
        }
        return ret ;
    }

    bool is_reads_too_little() const 
    {
        if ( ReadsNum() < min_reads_count )
            return true ;
        return false ;
    }

    bool is_reads_too_much() const 
    {
        if ( ReadsNum() > max_reads_count )
            return true ;
        return false ;
    }
};

struct ReadMatrixFactory
{
    static int max_error_count ;
    static int min_match_check ;
    static int max_reads_depth ;
    static int the_k ;
    static ReadAccessor * readAccessor ;

    static ReadMatrix GenReadMatrix(const Contig & contig , const ConsensusArea & area )
    {
        ReadMatrix ret ;
        ret.m_area = area ;
        int relative_num = 0 ;
        std::map<int,std::vector<ReadElement>> new_reads;
        while(true)
        {
            relative_num ++ ;
            if( relative_num == 1 )
            {
                UpdateKmerReadsFromContig(contig
                        ,area
                        ,ret
                        ,new_reads
                        );
            }
            else
            {
                std::map<int,std::vector<ReadElement>> prev_reads;
                std::swap(new_reads,prev_reads);
                UpdateKmerReadsFromReads(prev_reads,
                        area,
                        contig ,
                        relative_num ,
                        ret ,
                        new_reads);
            }

            if( ret.is_reads_too_much() )
                break ;
            if( new_reads.empty() )
                break ;
        };

        return ret ;
    }

    static std::vector<ReadElement> find_all_reads_start_with(
            const Number_t & kmer
            , const Contig & contig 
            , int position
            )
    {
        std::vector<ReadElement> ret ;

        LinkedList<ReadElement> reads ;
        reads.purge();
        readAccessor->getReadsBeginWith(kmer,reads,true,ReadAccessor::inAll);
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

    static void UpdateKmerReadsFromReads(
            const std::map<int ,std::vector<ReadElement>> & prev_reads ,
            const ConsensusArea & area ,
            const Contig & contig ,
            int index,
            ReadMatrix & ret,
            std::map<int ,std::vector<ReadElement>> & new_reads
            )
    {
        for( const auto & pair : prev_reads )
        {
            int pos = pair.first ;
            if( ret.is_reads_too_much() )
                break ;
            if( ! area.valid_starter( pos , the_k ) )
                continue ;
            for( const auto & read : pair.second )
            {
                TightString  tStrRead(read.getLen());
                read.getSequence(tStrRead);
                for( int i = 0 ; i < (int)read.getLen() - the_k ; i++ )
                {
                    if( ret.is_reads_too_much() )
                        break ;
                    if( ! area.valid_starter( pos + i , the_k ) )
                        break  ;
                    Number_t kmer;
                    tStrRead.readTightStringFragment
                        (i, i+the_k, kmer);
                    ret.tryInitPos(pos+i);
                    if( ret.checkKmerInPos( pos+i , kmer ) )
                        continue ;
                    else
                    {
                        auto reads = find_all_reads_start_with
                            (kmer ,contig , pos+ i );
                        if( reads.empty() )
                            continue ;
                        ret.AddKmer(kmer,reads,index,pos+i);
                        new_reads[i]= reads;
                    }
                }
            }
        }
    }

    static void UpdateKmerReadsFromContig(
            const Contig & contig,
            const ConsensusArea & area ,
            ReadMatrix & ret,
            std::map<int ,std::vector<ReadElement>> & new_reads
            )
    {

        const TightString & tStrContig = contig.getTightString() ;
        int pos_start = area.left_most_pos_in_contig ;
        int pos_end = (int)contig.getLength() - the_k + 1 ;
        //int pos_end1 = area.consensus_end_pos_in_contig  - the_k + 1 ;
        for( int i = pos_start ; i <= pos_end ;i++ ) /* i in 1 base */
        {
            Number_t kmer;
            tStrContig.readTightStringFragment
                (i-1, i+the_k-1, kmer);
            ret.tryInitPos(i);
            if( ret.checkKmerInPos(kmer,i) )
                continue ;
            else
            {
                auto reads = find_all_reads_start_with
                    (kmer ,contig , i );
                if( reads.empty() )
                    continue ;
                ret.AddKmer(kmer,reads,1,i);
                new_reads[i]= reads;
            }
        }
    }
};

#endif// READMATRIX_HPP__
