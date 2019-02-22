#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>

#include "Common.hpp"
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
    bool succ_found_reads ;
    //bool used ;
};

struct PositionInfoKeeper
{
    public:
        static std::vector<ReadElement> find_all_reads_start_with( const Number_t & kmer )
        {
            std::vector<ReadElement> ret ;
            // todo 
            return ret ;
        }

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

        void AddKmer(const Number_t & kmer, int relative_num )
        {
            auto itr = all_kmer2reads.find(kmer) ;
            if( itr == all_kmer2reads.end() )
            {
                auto & new_one = all_kmer2reads[kmer];
                new_one.the_kmer = kmer ;
                new_one.relative_num = relative_num ;
                new_one.reads = find_all_reads_start_with(kmer) ;
                if( new_one.reads.empty() )
                    new_one.succ_found_reads = false ;
                else 
                    new_one.succ_found_reads = true ;
                reads_num += new_one.reads.size() ;
            }
            else
            {
                ;
            }
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
};

struct ReadMatrixFactory
{
    static ReadMatrix GenReadMatrix(const Contig & contig , const ConsensusArea & area )
    {
        ReadMatrix ret ;
        ret.m_area = area ;

        //TODO
        return ret ;
    }
};

#endif// READMATRIX_HPP__
