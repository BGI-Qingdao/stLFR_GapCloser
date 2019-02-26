#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>

#include "Common.hpp"
#include "TightString.hpp"
#include "GlobalAccesser.hpp"
#include "readtool/ReadElement.hpp"
#include "readtool/ReadAccessor.hpp"
#include "Contig.hpp"
#include "ContigTool.hpp"
#include "ConsensusConfig.hpp"
#include "ContigForFill.hpp"

struct Kmer2Reads
{
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
                const ReadElement & reads ,
                int relative_num )
        {
            auto & new_one = all_kmer2reads[kmer] ;
            new_one.relative_num = relative_num ;
            new_one.reads.push_back(reads) ;
            reads_num ++ ;
        }
        void AddKmer(const Number_t & kmer,
                const std::vector<ReadElement> & reads ,
                int relative_num )
        {
            auto & new_one = all_kmer2reads[kmer] ;
            new_one.relative_num = relative_num ;
            new_one.reads = reads ;
            reads_num += new_one.reads.size() ;
        }

        PositionInfoKeeper GetSubReadsByPE( int pos , const Contig & prev_contig ) const
        {
            PositionInfoKeeper ret;
            for( const auto & pair : all_kmer2reads )
            {
                const Number_t & kmer = pair.first ;
                const Kmer2Reads & reads = pair.second ;
                for( const auto & a_read : reads.reads )
                {

                    if( ContigTool::IsEligiblePECheckRead(a_read, prev_contig,pos) )
                    {
                        ret.AddKmer(kmer,a_read, reads.relative_num);
                    }
                }
            }
            return ret ;
        }

        PositionInfoKeeper GetSubReadsByBarcode(
                const Contig & prev_contig ,
                const Contig & next_contig 
                ) const
        {
            PositionInfoKeeper ret;
            for( const auto & pair : all_kmer2reads )
            {
                const Number_t & kmer = pair.first ;
                const Kmer2Reads & reads = pair.second ;
                for( const auto & a_read : reads.reads )
                {
                    if( ContigTool::IsEligibleBarcodeCheckRead(a_read, prev_contig,next_contig) )
                    {
                        ret.AddKmer(kmer,a_read, reads.relative_num);
                    }
                }
            }
            return ret ;
        }

        PositionInfoKeeper GetSubReadsByONT(
                const Contig & prev_contig ,
                const Contig & next_contig 
                ) const 
        {
            PositionInfoKeeper ret;
            //TODO
            return ret ;
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
    public :
        static int  max_reads_count ;
        static int  min_reads_count ;

        static int min_sub_reads_count;

    private:
        //  position --> reads
        std::map< int  , PositionInfoKeeper> m_raw_reads;

        ConsensusArea m_area;


        static ReadMatrix Merge( const ReadMatrix & /*left*/ , const ReadMatrix & /*right*/ )
        {
            ReadMatrix ret ;
            //TODO 
            return ret ;
        }

        ReadMatrix GetSubMatrixByPECheck(const Contig & prev_contig)
        {
            ReadMatrix ret ;
            for( const auto & pair : m_raw_reads )
            {
                int pos = pair.first ;
                const PositionInfoKeeper & reads  =  pair.second ;
                PositionInfoKeeper sub = reads.GetSubReadsByPE( pos , prev_contig);
                if( sub.ReadsNum() > 0 )
                    ret.m_raw_reads[pos]=sub ;
            }
            return ret ;
        }

        ReadMatrix GetSubMatrixByBarcodeCheck(
                const Contig & prev_contig,
                const Contig & next_contig
                )
        {
            ReadMatrix ret ;
            for( const auto & pair : m_raw_reads )
            {
                int pos = pair.first ;
                const PositionInfoKeeper & reads  =  pair.second ;
                PositionInfoKeeper sub = reads.GetSubReadsByBarcode(
                        prev_contig , next_contig);
                if( sub.ReadsNum() > 0 )
                    ret.m_raw_reads[pos]=sub ;
            }
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

    public: 

        ConsensusMatrix GenConsensusMatrix(  )
        {
            ConsensusMatrix ret ;
            //TODO
            return ret ;
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

        void AddKmer(const Number_t & kmer,
                const std::vector<ReadElement> & reads , 
                int relative_num , int pos )
        {
            m_raw_reads[pos].AddKmer(kmer,reads,relative_num);
        }


        bool checkKmerInPos(const Number_t & kmer , int pos )
        {
            if( m_raw_reads.find(pos) == m_raw_reads.end())
                return false ;
            return m_raw_reads[pos].KmerExist(kmer);
        }

        void tryInitPos(int pos)
        {
            if( m_raw_reads.find(pos) == m_raw_reads.end())
                m_raw_reads[pos].Init() ;
        }

        void AddArea(const ConsensusArea & area)
        {
            m_area = area ;
        }
};

struct ReadMatrixFactory
{
    static int max_error_count ;
    static int min_match_check ;
    static int max_reads_depth ;
    static int the_k ;

    static ReadMatrix GenReadMatrix(const Contig & contig , const ConsensusArea & area )
    {
        ReadMatrix ret ;
        ret.AddArea( area) ;
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
                        auto reads = ContigTool::find_all_reads_start_with
                            (kmer ,contig , pos+ i 
                             ,max_reads_depth
                             ,max_error_count);
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
