#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <tuple>

#include "Common.hpp"
#include "ConsensusConfig.hpp"
#include "TightString.hpp"
#include "GlobalAccesser.hpp"
#include "ContigForFill.hpp"
#include "readtool/ReadElement.hpp"
#include "readtool/ReadAccessor.hpp"
#include "Contig.hpp"
#include "ContigTool.hpp"

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

        void AddKmer(const Number_t & kmer,
                const ReadElement & reads ,
                int relative_num )
        {
            auto & new_one = all_kmer2reads[kmer] ;
            new_one.relative_num = relative_num ;
            new_one.reads.push_back(reads) ;
            reads_num ++ ;
        }

    public :

        std::map<Number_t , Kmer2Reads>::const_iterator begin() const
        {
            return all_kmer2reads.begin() ;
        }

        std::map<Number_t , Kmer2Reads>::const_iterator end()  const
        {
            return all_kmer2reads.end() ;
        }

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
                    if( ContigTool::IsEligibleBarcodeCheckRead
                            (a_read, prev_contig,next_contig) )
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
            for( const auto & pair : all_kmer2reads )
            {
                const Number_t & kmer = pair.first ;
                const Kmer2Reads & reads = pair.second ;
                for( const auto & a_read : reads.reads )
                {
                    if( ContigTool::IsEligibleONTCheckRead
                            (a_read, prev_contig,next_contig) )
                    {
                        ret.AddKmer(kmer,a_read, reads.relative_num);
                    }
                }
            }
            return ret ;
        }
};

struct ConsensusResult
{
    std::string       consensus_str;
    std::vector<bool> conflict_flags;
    std::vector<bool> depth_flags;

    TightString getTightString() const 
    {
        TightString ret(consensus_str.c_str() , 
                consensus_str.length());
        return ret ;
    }

    TightString getTightString(int pos , int len ) const 
    {
        TightString ret(consensus_str.c_str() + pos  , 
                len);
        return ret ;
    }

    int getLength() const {
        return consensus_str.length() ;
    }

    int conflict_num() const {
        int i = 0 ;
        for( bool x : conflict_flags )
            if( !x )
                i++ ;

        return i ;
    }

    int low_depth() const {
        int i = 0 ;
        for( bool x : depth_flags )
            if( !x )
                i++ ;

        return i ;
    }

    bool is_consensus_done() const 
    {
        return  (conflict_num() <= Threshold::max_allowed_conflict )
            && ( low_depth() <= Threshold::max_low_depth)
            ;
    }
};


struct ConsensusMatrix
{
    private:
        std::vector<std::array<int,4>> depth_matrix;

        std::tuple<char , bool , bool > consensus_pos(const std::array<int,4> &  a) const
        {
            int total_depth = 0 ;
            int highest_depeth = 0 ;
            int highest_nucleotide = 0 ;
            for( int i = 0 ; i < 4; i ++ )
            {
                total_depth += a[i] ;
                if( highest_depeth < a[i] )
                {
                    highest_depeth = a[i] ;
                    highest_nucleotide = i ; 
                }
            }
            bool is_low_depth = (total_depth < Threshold::min_nucleotide_depth );
            if ( total_depth == 0 )
                return std::make_tuple( 'N' , true , ! is_low_depth ) ;
            char ret =  numberToNucleotide(highest_nucleotide);
            if( highest_depeth >= (float)total_depth * Threshold::NoConflictThreshold )
                return std::make_tuple( ret , true , ! is_low_depth ) ;
            else 
                return std::make_tuple( ret , false , ! is_low_depth ) ;
        }

    public:

        void Init(int length)
        {
            depth_matrix.resize(length);
            for( auto & m : depth_matrix )
            {
                m[0] = 0 , m[1] = 0 , m[2] =0 , m[3] = 0 ;
            }
        }

        ConsensusResult GenConsensusResult() const
        {
            ConsensusResult ret ;
            for ( const auto & m : depth_matrix )
            {
                char nucleotide ;
                bool not_confilict ;
                bool not_low_depth ;
                std::tie( nucleotide,not_confilict,not_low_depth) 
                    = consensus_pos( m );
                if( nucleotide == 'N' )
                    break ;
                ret.conflict_flags.push_back(not_confilict);
                ret.depth_flags.push_back(not_low_depth);
                ret.consensus_str.push_back(nucleotide);
            }
            return ret ;
        }

        void UpdateDepth( int pos , int nucleotide , int depth = 1)
        {
            assert(pos >= 0 && pos < (int)depth_matrix.size() );
            assert( nucleotide < 4 );
            depth_matrix[pos][nucleotide] += depth ;
        }
};

struct ReadMatrix
{
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

    public:

        const ConsensusArea & Area() const {
            return m_area ;
        }

        bool HasReads( int pos ) const {
            return m_raw_reads.find(pos) != m_raw_reads.end() ;
        }

        const PositionInfoKeeper & GetInfoKeeper(int pos ) const
        {
            static PositionInfoKeeper tmp ;
            if( ! HasReads( pos ) )
                return tmp ;
            return m_raw_reads.at(pos);
        }

        ReadMatrix GenSubMatrixByGap(const Contig & prev_contig ,
                const Contig & next_contig ,
                const GapInfo & gap )
        {

            if( ReadsNum() < Threshold::min_sub_reads_count )
                return *this ;
            if( gap.is_gap_big() ) {;} else {;}
            //{
            auto sub1 = GetSubMatrixByPECheck( prev_contig) ;
            auto sub2 = GetSubMatrixByBarcodeCheck(prev_contig , next_contig );
            auto sub3 = Merge( sub1 , sub2 );
            if( sub3.ReadsNum() < Threshold::min_sub_reads_count )
                return *this ;
            else 
                return sub3 ;
            //}
            //else
            //{

            //}
            return *this ;
        }

        ConsensusMatrix GenConsensusMatrix( const Contig & contig )
        {
            ConsensusMatrix ret ;
            ret.Init(m_area.consensus_end_pos_in_contig
                    -m_area.consensus_start_pos_in_contig
                    + 1 ) ;

            for( const auto & pair : m_raw_reads)
            {
                int pos = pair.first ;
                // update all reads's depth
                for( const auto & k_r_pair : pair.second )
                {
                    const auto & reads = k_r_pair.second.reads ;
                    for( const auto & read : reads )
                    {
                        TightString read_str ;
                        read.getSequence(read_str);
                        for( int i = 0 ; i < (int)read.getLen() ; i++ )
                        {
                            int matrix_pos = m_area.pos_translate_contig2martix(pos);
                            if( matrix_pos < 0 )
                                continue ;
                            ret.UpdateDepth(matrix_pos , read_str[i]);
                        }
                    }
                }
            }
            int end = m_area.consensus_end_pos_in_contig ;
            if( end >= (int) contig.getLength() -1 )
                end = contig.getLength() -1 ;
            // update the contig's depth
            for( int i = m_area.consensus_start_pos_in_contig 
                    ; i <= end ; i ++ )
            {
                int matrix_pos = m_area.pos_translate_contig2martix(i);
                ret.UpdateDepth( matrix_pos
                        , contig.getTightString()[i]
                        /*, contig.getDepths()[i] */ 
                        /*new logic that contig depth always 1 */
                        );
            }
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
            if ( ReadsNum() < Threshold::min_reads_count )
                return true ;
            return false ;
        }

        bool is_reads_too_much() const 
        {
            if ( ReadsNum() > Threshold::max_reads_count )
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
            if( ! area.valid_starter( pos ,Threshold::the_k ) )
                continue ;
            for( const auto & read : pair.second )
            {
                TightString  tStrRead(read.getLen());
                read.getSequence(tStrRead);
                for( int i = 0 ; i < (int)read.getLen() - Threshold::the_k ; i++ )
                {
                    if( ret.is_reads_too_much() )
                        break ;
                    if( ! area.valid_starter( pos + i , Threshold::the_k ) )
                        break  ;
                    Number_t kmer;
                    tStrRead.readTightStringFragment
                        (i, i+Threshold::the_k, kmer);
                    ret.tryInitPos(pos+i);
                    if( ret.checkKmerInPos( pos+i , kmer ) )
                        continue ;
                    else
                    {
                        auto reads = ContigTool::find_all_reads_start_with
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
        int pos_end = (int)contig.getLength() - Threshold::the_k + 1 ;
        //int pos_end1 = area.consensus_end_pos_in_contig  - the_k + 1 ;
        for( int i = pos_start ; i <= pos_end ;i++ ) /* i in 1 base */
        {
            Number_t kmer;
            tStrContig.readTightStringFragment
                (i-1, i+Threshold::the_k-1, kmer);
            ret.tryInitPos(i);
            if( ret.checkKmerInPos(kmer,i) )
                continue ;
            else
            {
                auto reads = ContigTool::find_all_reads_start_with
                    (kmer ,contig , i);
                if( reads.empty() )
                    continue ;
                ret.AddKmer(kmer,reads,1,i);
                new_reads[i]= reads;
            }
        }
    }
};

#endif// READMATRIX_HPP__
