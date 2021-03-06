#ifndef READMATRIX_HPP__
#define READMATRIX_HPP__

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <tuple>

#include "Common.hpp"
#include "Misc.h"
#include "ConsensusConfig.hpp"
#include "TightString.hpp"
#include "GlobalAccesser.hpp"
#include "ContigForFill.hpp"
#include "readtool/ReadElement.hpp"
#include "readtool/ReadAccessor.hpp"
#include "Contig.hpp"
#include "ContigTool.hpp"
#include "ConsensusLog.h"

enum SubReadSetType
{
    Unknow = 0 ,
    PE = 1 ,
    PE_Barcode = 2 ,
    Empty = 3 ,
    Basic = 4
};

struct Kmer2Reads
{
    std::vector<ReadElement> reads;
    int relative_num ;
    int ReadsNum() const { return reads.size() ; }
};

struct PositionInfoKeeper
{
    private:

        std::map<Number_t , Kmer2Reads> all_kmer2reads;

//        int reads_num ;

        void AddKmer(const Number_t & kmer,
                const ReadElement & reads ,
                int relative_num )
        {
            auto & new_one = all_kmer2reads[kmer] ;
            new_one.relative_num = relative_num ;
            new_one.reads.push_back(reads) ;
//            reads_num ++ ;
        }

    public :

        static PositionInfoKeeper Merge( const PositionInfoKeeper & left 
                , const PositionInfoKeeper & right )
        {
            PositionInfoKeeper ret ;
            ret.Init()  ;
            // add left only
            for( const auto & pair1 : left.all_kmer2reads )
            {
                if( right.all_kmer2reads.find(pair1.first) == right.all_kmer2reads.end() )
                {
                    if( ! pair1.second.reads.empty() )
                    {
                        ret.AddKmer( pair1.first, pair1.second.reads, pair1.second.relative_num );
                    }
                }
            }
            // add right only
            for( const auto & pair2 : right.all_kmer2reads )
            {
                if( left.all_kmer2reads.find(pair2.first) == left.all_kmer2reads.end() )
                {
                    if( ! pair2.second.reads.empty() )
                    {
                        ret.AddKmer( pair2.first, pair2.second.reads, pair2.second.relative_num );
                    }
                }
            }
            // merge right && right part 
            for( const auto & pair1 : left.all_kmer2reads )
            {
                if( right.all_kmer2reads.find(pair1.first) != right.all_kmer2reads.end() )
                {
                    std::set<ReadElement> all ;
                    for( const auto & i : pair1.second.reads )
                        all.insert(i);
                    for( const auto & j : right.all_kmer2reads.at(pair1.first).reads )
                        all.insert(j);
                    std::vector<ReadElement> items ;
                    for( const auto & m : all )
                        items.push_back(m);
                    if( ! items.empty() )
                        ret.AddKmer( pair1.first, items , pair1.second.relative_num );
                }
            }
            return ret ;
        }

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
//            reads_num = 0 ;
        }

        int ReadsNum() const {
//            return reads_num ;
            int ret = 0 ;
            for( const auto & i : all_kmer2reads )
            {
                ret += i.second.ReadsNum() ;
            }
            return ret ;
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
//            reads_num += new_one.reads.size() ;
        }

        PositionInfoKeeper GetSubReadsByPE( int pos , const Contig & prev_contig ) const
        {
            PositionInfoKeeper ret;
            ret.Init();
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
            ret.Init();
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
            ret.Init();
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

    TightString getSubTightString(
            int skip_len ,
            int len 
            ) const 
    {
        TightString ret(consensus_str.c_str() +skip_len ,
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

    int consensus_length(const SubReadSetType type ) const  /*1base ret*/
    {
        if( type == SubReadSetType::Empty || consensus_str.size() <1 )
            return -1;
        int ret = 0 ;
        int conflict_nums = 0 ;
        for ( ; ret <(int) consensus_str.size() ; ret ++ )
        {
            if( ! depth_flags.at(ret) )
                break ;
            if( ! conflict_flags.at(ret) )
                conflict_nums ++ ;
            if( type == SubReadSetType::Basic )
                if( conflict_nums > Threshold::basic_set_max_conflict )
                    break ;
                else ;
            else
                if( conflict_nums > Threshold::max_allowed_conflict )
                    break ;
                else ;
        }
        assert( ret <= (int) consensus_str.size() && ret >= 0 );
        return ret ;
    }
};


struct ConsensusMatrix
{
    private:
        std::vector<std::array<int,4>> depth_matrix;

        std::tuple<char , bool , bool , int> consensus_pos(const std::array<int,4> &  a , const SubReadSetType type ) const
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
            int lt = 0 ;
            float th = 1.0f ;
            if ( type != SubReadSetType::Basic )
            {
                th = Threshold::NoConflictThreshold ;
                lt = Threshold::min_nucleotide_depth ;
            }
            else 
            {
                th = Threshold::basic_NoConflictThreshold ;
                lt = Threshold::basic_min_nucleotide_depth ;
            }
            bool is_low_depth = (total_depth < lt );
            if( total_depth == highest_depeth )
                is_low_depth = false ;
            if ( total_depth == 0 )
                return std::make_tuple( 'N' , true , ! is_low_depth,total_depth ) ;
            char ret =  numberToNucleotide(highest_nucleotide);
            if( highest_depeth >= (float)total_depth * th )
                return std::make_tuple( ret , true , ! is_low_depth , total_depth) ;
            else 
                return std::make_tuple( ret , false , ! is_low_depth, total_depth ) ;
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

        ConsensusResult GenConsensusResult( const SubReadSetType type
                , std::vector<LOG::PosInfo> & details ) const
        {
            ConsensusResult ret ;
            for ( const auto & m : depth_matrix )
            {
                char nucleotide ;
                bool not_confilict ;
                bool not_low_depth ;
                int total_depth ;
                std::tie( nucleotide,not_confilict,not_low_depth , total_depth)
                    = consensus_pos( m ,type );
                if( nucleotide == 'N' )
                {
                    details.push_back( { total_depth , 'N' } );
                    break ;
                }
                ret.conflict_flags.push_back(not_confilict);
                ret.depth_flags.push_back(not_low_depth);
                ret.consensus_str.push_back(nucleotide);
                if( not_confilict && not_low_depth )
                    details.push_back( { total_depth , 'M' } );
                else if ( ! not_confilict )
                    details.push_back( { total_depth , 'C' } );
                else 
                    details.push_back( { total_depth , 'L' } );
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
        std::map< int  , PositionInfoKeeper> m_raw_reads; /* 1base*/

        ConsensusArea m_area;

        static ReadMatrix Merge( const ReadMatrix & left , const ReadMatrix & right )
        {
            ReadMatrix ret ;
            // add left only
            for( const auto & pair1 : left.m_raw_reads )
            {
                if( right.m_raw_reads.find(pair1.first) == right.m_raw_reads.end() )
                {
                    if( pair1.second.ReadsNum() > 0 )
                        ret.m_raw_reads[pair1.first] = pair1.second ;
                }
            }
            // add right only
            for( const auto & pair2 : right.m_raw_reads )
            {
                if( left.m_raw_reads.find(pair2.first) == left.m_raw_reads.end() )
                {
                    if( pair2.second.ReadsNum() > 0 )
                        ret.m_raw_reads[pair2.first] = pair2.second ;
                }
            }
            // merge left & right both part
            for( const auto & pair2 : left.m_raw_reads )
            {
                if( right.m_raw_reads.find(pair2.first) != right.m_raw_reads.end() )
                {
                    ret.m_raw_reads[pair2.first] =
                        PositionInfoKeeper::Merge(
                                pair2.second
                                ,right.m_raw_reads.at(pair2.first)
                                );
                }
            }
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

        ReadMatrix GenSubMatrix(const Contig & prev_contig ,
                const Contig & next_contig ,
                /*const GapInfo & gap ,*/
                SubReadSetType & ret_type ,
                LOG::ConsensusLog & log
                )
        {
            SubReadsLog tmp ;
            tmp.Init();
            tmp.basic_num = ReadsNum() ;
            log.size_b = tmp.basic_num ;
            // try PE first 
            auto sub1 = GetSubMatrixByPECheck( prev_contig) ;
            int s1 = sub1.ReadsNum() ;
            tmp.pe_num = s1 ;
            log.size_p  = s1 ;
            if( s1 >= Threshold::min_pe_sub_reads_count )
            {
                tmp.used_num = s1 ;
                tmp.type = "PE" ;
                log.used_id = tmp.type;
                GlobalAccesser::sub_read_num.Touch(tmp);
                GlobalAccesser::sub_type.Touch("PE");
                sub1.m_area = m_area ;
                ret_type = SubReadSetType::PE ;
                return sub1 ;
            }
            // add Barcode and try again
            auto sub2 = GetSubMatrixByBarcodeCheck(prev_contig , next_contig );
            int s2 = sub2.ReadsNum() ;
            auto sub3 = Merge( sub1 , sub2 );
            int s13 = sub3.ReadsNum() ;
            tmp.barcode_num = s2 ;
            log.size_bc = s2 ;
            tmp.pe_barcode_num = s13 ;
            log.size_bp = s13 ;
            if( s13>= Threshold::min_pe_barcode_sub_reads_count )
            {
                tmp.used_num = s13 ;
                tmp.type = "PE_BARCODE" ;
                log.used_id = tmp.type;
                GlobalAccesser::sub_read_num.Touch(tmp);
                GlobalAccesser::sub_type.Touch("PE_BARCODE");
                ret_type = SubReadSetType::PE_Barcode ;
                sub3.m_area = m_area ;
                return sub3;
            }
            // return null or basic
            if( Threshold::use_subset_only )
            {
                tmp.used_num = 0 ;
                tmp.type = "NULL" ;
                log.used_id = tmp.type;
                GlobalAccesser::sub_read_num.Touch(tmp);
                GlobalAccesser::sub_type.Touch("NULL");
                ReadMatrix tmp ;
                ret_type = SubReadSetType::Empty ;
                tmp.m_area = m_area ;
                return tmp;
            }
            else
            {
                tmp.used_num = tmp.basic_num  ;
                tmp.type = "BASIC" ;
                log.used_id = tmp.type;
                GlobalAccesser::sub_read_num.Touch(tmp);
                GlobalAccesser::sub_type.Touch("BASIC");
                ret_type = SubReadSetType::Basic ;
                return *this ;
            }
        }

        ConsensusMatrix GenConsensusMatrix( const Contig & contig )
        {
            ConsensusMatrix ret ;
            ret.Init(m_area.consensus_end_pos_in_contig
                    -m_area.consensus_start_pos_in_contig
                    + 1 ) ;

            for( const auto & pair : m_raw_reads )
            {
                int pos = pair.first ; /*1base*/
                // update all reads's depth
                for( const auto & k_r_pair : pair.second )
                {
                    const auto & reads = k_r_pair.second.reads ;
                    for( const auto & read : reads )
                    {
                        TightString read_str(read.getLen()) ;
                        read.getSequence(read_str);
                        for( int i = 0 ; i < (int)read.getLen() ; i++ )
                        {
                            int matrix_pos = m_area.pos_translate_contig2martix(pos+i);
                            if( matrix_pos < 0 )
                                continue ;
                            ret.UpdateDepth(matrix_pos , read_str[i]);
                        }
                    }
                }
            }

            int end = m_area.consensus_end_pos_in_contig ;
            if( end >= (int) contig.getLength()  )
                end = contig.getLength() ;
            // update the contig's depth
            for( int i = m_area.consensus_start_pos_in_contig 
                    ; i <= end ; i ++ ) /*1 in 1 base */
            {
                int matrix_pos = m_area.pos_translate_contig2martix(i);
                ret.UpdateDepth( matrix_pos
                        , contig.getTightString()[i-1]
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
            if ( ReadsNum() >= Threshold::max_reads_count )
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
        std::set<Number_t> used_kmers;
        ReadMatrix ret ;
        ret.AddArea( area) ;
        int relative_num = 0 ;
        // 1 base position --> read set
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
                        ,used_kmers
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
                        new_reads,
                        used_kmers 
                        );
            }

            if( ret.is_reads_too_much() )
                break ;
            if( new_reads.empty() )
                break ;
            if( relative_num > Threshold::max_reads_round )
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
            std::map<int ,std::vector<ReadElement>> & new_reads,
            std::set<Number_t> & used_kmers
            )
    {
        for( const auto & pair : prev_reads )
        {
            if( ret.is_reads_too_much() )
                break ;
            int pos = pair.first ;
            if( ! area.valid_starter( pos ,Threshold::the_k ) )
                continue ;
            for( const auto & read : pair.second )
            {
                if( ret.is_reads_too_much() )
                    break ;
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
                    if( used_kmers.find( kmer ) != used_kmers.end () )
                        continue ;
                    ret.tryInitPos(pos+i);
                    if( ret.checkKmerInPos(  kmer , pos+i ) )
                        continue ;
                    else
                    {
                        auto reads = ContigTool::find_all_reads_start_with
                            (kmer ,contig , pos+ i );
                        if( reads.empty() )
                            continue ;
                        if((int)reads.size() > Threshold::max_kmer_2_read )
                            continue ;
                        GlobalAccesser::kmer_read_count.Touch(reads.size());
                        ret.AddKmer(kmer,reads,index,pos+i);
                        new_reads[pos + i]= reads;
                        used_kmers.insert(kmer);
                    }
                }
            }
        }
    }

    static void UpdateKmerReadsFromContig(
            const Contig & contig,
            const ConsensusArea & area ,
            ReadMatrix & ret,
            std::map<int ,std::vector<ReadElement>> & new_reads ,
            std::set<Number_t> & used_kmers
            )
    {

        const TightString & tStrContig = contig.getTightString() ;
        int pos_start = area.left_most_pos_in_contig ;
        int pos_end = (int)contig.getLength() - Threshold::the_k + 1 ;
        int pos_end1 = area.right_most_pos_in_contig - Threshold::the_k+ 1 ;
        if( pos_end1 < pos_end )
            pos_end = pos_end1;
        for( int i = pos_start ; i <= pos_end ;i++ ) /* i in 1 base */
        {
            if ( ret.is_reads_too_much() )
                break ;
            Number_t kmer;
            tStrContig.readTightStringFragment
                (ConsensusArea::pos2index(i)
                 ,ConsensusArea::pos2index(i+Threshold::the_k)
                 , kmer); /* 0 base */
            if( used_kmers.find( kmer ) != used_kmers.end () )
                continue ;
            ret.tryInitPos(i);

            if( ret.checkKmerInPos(kmer,i) )
                continue ;
            else
            {
                auto reads = ContigTool::find_all_reads_start_with
                    (kmer ,contig , i);
                if( reads.empty() )
                    continue ;
                if((int)reads.size() > Threshold::max_kmer_2_read )
                    continue ;
                GlobalAccesser::kmer_read_count.Touch(reads.size());
                ret.AddKmer(kmer,reads,1,i);
                new_reads[i]= reads;
                used_kmers.insert(kmer);
            }
        }
    }
};

#endif// READMATRIX_HPP__
