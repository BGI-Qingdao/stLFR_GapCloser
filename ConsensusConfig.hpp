#ifndef CONSENSUSCONFIG_HPP__
#define CONSENSUSCONFIG_HPP__

#include <cassert>

struct ConsensusArea
{
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

    int pos_translate_contig2martix(int pos) const
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

    bool valid_starter( int pos , int kvalue ) const
    {
        return pos + kvalue < consensus_end_pos_in_contig ;
    }
};

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
        assert( extra_len > 0 );
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

#endif // CONSENSUSCONFIG_HPP__
