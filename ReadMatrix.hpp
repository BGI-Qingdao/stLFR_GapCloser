#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>
#include <set>

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
    public:

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
        bool AddKmer(const Number_t & kmer, std::vector<ReadElement> & reads , int relative_num )
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

    //  position --> reads
    std::map< int  , PositionInfoKeeper> m_raw_reads;

    ConsensusArea m_area;

    ConsensusMatrix GenConsensusMatrix(  )
    {
        ConsensusMatrix ret ;
        //TODO
        return ret ;
    }

    ReadMatrix GenSubMatrixByGap(const Contig & prev_contig ,
            const Contig & next_contig , 
            const GapInfo & gap )
    {
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

    int ReadsNum()
    {
        int ret = 0 ;
        for( const auto & pair : m_raw_reads)
        {
            ret += pair.second.ReadsNum() ;
        }
    }

    bool is_reads_too_little() const 
    {

    }

    bool is_reads_too_much() const 
    {

    }

};

struct ReadMatrixFactory
{
    static int max_error_count ;
    static int min_match_check ;
    static int max_reads_depth ;
    static ReadAccessor * readAccessor ;

    static ReadMatrix GenReadMatrix(const Contig & contig , const ConsensusArea & area )
    {
        ReadMatrix ret ;
        ret.m_area = area ;
        int relative_num = 0 ;
        std::set<Number_t> unsed_kmers ;
        std::set<ReadElement> new_reads;
        while(true)
        {
            relative_num ++ ;
            if( relative_num == 1 )
            {
                // chop kmer from contig
 //               for( int i = area.
            }
            else
            {
                // chop kmer from new reads

            }

            if( ret.is_reads_too_much() )
                break ;
            if( new_reads.empty() )
                break ;
        };


        //TODO
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
            int unmatch =  CheckUnmatchNum( contig.getTightString(),
                        contig.getLength() ,
                        read_str ,
                        readElement.getLen() ,
                        position ,
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
        for (Len_t m= 0; m< checked_num ; m++)
            if ( ref[contig_start + m ] != read[m] )
                ret ++ ;
        return ret ;
    }

};

#endif// READMATRIX_HPP__
