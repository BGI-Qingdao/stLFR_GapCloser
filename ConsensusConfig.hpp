#ifndef CONSENSUSCONFIG_HPP__
#define CONSENSUSCONFIG_HPP__

#include <cassert>

struct ConsensusArea
{
    int left_most_pos_in_contig ;  /*1base*/
    int right_most_pos_in_contig ; /*1base*/
    int consensus_start_pos_in_contig;/*1base*/
    int consensus_end_pos_in_contig ;/*1base*/
    int contig_len ;


    ConsensusArea() : 
        left_most_pos_in_contig(-1),
        right_most_pos_in_contig(-1),
        consensus_start_pos_in_contig(-1),
        consensus_end_pos_in_contig(-1),
        contig_len(-1){}

    bool operator == ( const ConsensusArea & other )
    {
        return left_most_pos_in_contig == other.left_most_pos_in_contig 
            && right_most_pos_in_contig == other.right_most_pos_in_contig 
            && consensus_start_pos_in_contig == other.consensus_start_pos_in_contig 
            && consensus_end_pos_in_contig == other.consensus_end_pos_in_contig 
            && contig_len == other.contig_len ;
    }

    static int pos2index( int pos /*1base*/ )
    {
        return pos -1 ;
    }

    int contig2matrix_translate_pos( int index /* 0base */) const
    {
        assert( index >= 0 );
        return index + consensus_start_pos_in_contig ; 
    }

    int pos_translate_contig2martix(int pos /* in 1 base */ ) const
    {
        if( pos < consensus_start_pos_in_contig 
                || pos > consensus_end_pos_in_contig )
        {
            return -1 ;
        }
        else
        {
            return pos - consensus_start_pos_in_contig ;
        }
    }

    bool valid_starter( int pos /* in 1 base */ , int kvalue ) const
    {
        return pos + kvalue < right_most_pos_in_contig;
    }
    int externed_len() const 
    {
        return consensus_end_pos_in_contig - contig_len;
    }
    int remained_len() const
    {
        return contig_len - consensus_start_pos_in_contig + 1 ;
    }

    int total_consensus_len() const 
    {
        return consensus_end_pos_in_contig - consensus_start_pos_in_contig + 1 ;
    }
};

//             | <-----------reads------>|
//             |               0 pos     |
//             |               |         |
// contig      | x1            |       y1|
//<<----------------------------
//                 |x2             y2|
//                 |                 |
//                 |<---consensus--->|
//
struct NewConsensusConfig
{
    static int x1 ;
    static int y1 ;
    static int x2 ;
    static int y2 ;

    static ConsensusArea  GetConsensusArea( int contiglen )
    {
        ConsensusArea ret ;
        ret.contig_len = contiglen ;

        // Calc the reads area.
        ret.left_most_pos_in_contig =
            contiglen - x1 > 1 ? contiglen - x1 : 1  ;
        ret.right_most_pos_in_contig =
            contiglen - y1 > 1 ?  contiglen - y1 : 1 ;
        assert( ret.left_most_pos_in_contig < ret.right_most_pos_in_contig );
        // Calc the consensus area.
        ret.consensus_start_pos_in_contig =
            contiglen - x2 > 1 ? contiglen - x2 : 1  ;
        ret.consensus_end_pos_in_contig =
            contiglen - y2 > 1 ?  contiglen - y2 : 1 ;
        assert( ret.left_most_pos_in_contig < ret.right_most_pos_in_contig );
        return ret ;
    }
};

#endif // CONSENSUSCONFIG_HPP__
