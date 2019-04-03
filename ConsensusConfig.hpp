#ifndef CONSENSUSCONFIG_HPP__
#define CONSENSUSCONFIG_HPP__

#include <cassert>

struct ConsensusArea
{
    /* all position index in 1 base */
    int left_most_pos_in_contig ;
    int right_most_pos_in_contig ;
    int consensus_start_pos_in_contig;
    int consensus_end_pos_in_contig ;
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

    static int pos2index( int pos )
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
        return  pos >= left_most_pos_in_contig &&
            ( pos + kvalue -1 ) <= right_most_pos_in_contig;
    }

    int reads_area_len() const 
    {
        return right_most_pos_in_contig - left_most_pos_in_contig + 1;
    }

    int consensus_area_len() const 
    {
        return consensus_end_pos_in_contig -  consensus_start_pos_in_contig + 1 ;
    }

    bool IsValid( int kvalue ) const 
    {
        return ( reads_area_len() >= kvalue )
            && ( left_most_pos_in_contig
                <= consensus_start_pos_in_contig )
            && ( -left_most_pos_in_contig + contig_len +1
               >= kvalue );
    }
};

/*
struct ConsensusConfig
{
    static int extend_len ;
    static int prev_extra_len;
    static int last_extra_len;
    static int consensus_len ;

    static ConsensusArea  GetConsensusArea( int contiglen )
    {
        assert( contiglen > 0 );
        assert( extend_len > 0 );
        assert( prev_extra_len > 0 );
        assert( last_extra_len > 0 );
        assert( consensus_len > 0 );

        ConsensusArea ret ;
        int in_contig = consensus_len - extend_len ;

        ret.contig_len = contiglen ;
        ret.left_most_pos_in_contig = contiglen - ( prev_extra_len + in_contig ) + 1 ;
        ret.right_most_pos_in_contig = contiglen + extend_len + last_extra_len ;
        ret.consensus_start_pos_in_contig = contiglen - in_contig + 1 ;
        ret.consensus_end_pos_in_contig = contiglen + extend_len ;

        if ( ret.left_most_pos_in_contig <= 0 )
            ret.left_most_pos_in_contig = 1;
        if ( ret.consensus_start_pos_in_contig <= 0 )
            ret.consensus_start_pos_in_contig = 1 ;


        return ret ;
    }
};
*/
//                        0
//                        |
//          x1            |       y1
//<<-----------------------
//             x2             y2
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
            contiglen - x1 > 0 ? contiglen - x1 : 0  ;
        ret.right_most_pos_in_contig =
            contiglen - y1 > 0 ?  contiglen - y1 : 0 ;

        // Calc the consensus area.
        ret.consensus_start_pos_in_contig =
            contiglen - x2 > 0 ? contiglen - x2 : 0  ;
        ret.consensus_end_pos_in_contig =
            contiglen - y2 > 0 ?  contiglen - y2 : 0 ;

        return ret ;
    }
};

#endif // CONSENSUSCONFIG_HPP__
