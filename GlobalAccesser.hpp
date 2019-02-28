#ifndef GLOBALACCESSER_HPP__
#define GLOBALACCESSER_HPP__

class ReadAccessor ;
class PairInfo ;

struct GlobalAccesser
{
    static ReadAccessor * the_read_accessor;
    static PairInfo * the_pair_info ;
};

struct Threshold
{
    // the minumum proportion that a nucleotide should
    //  hold to be a major nucletide.
    static float    NoConflictThreshold ;

    // the maximum number of confilicts that a consensus
    //  can accept.
    static int      max_allowed_conflict ;

    // the maximum number of reads that a matrix can hold.
    static int      max_reads_count ;

    // the minimum number of reads that can start a consensus .
    //  otherwise , can not consensus.
    static int      min_reads_count ;

    // the minimum number of reads that can use sub matrix to
    //  consensus , otherwhise use the basic matrix .
    static int      min_sub_reads_count;

    // the maximum number of confilicts that a map
    //  can accept.
    static int      max_error_count ;

    // the maxinum number of depth that a read can be accept.
    static int      max_reads_depth ;

    // the kvalue for kmer-read map.
    static int      the_k ;

    // the minimum depth of a nucleotide that 
    //  make it a NOT to be a low depth nucleotide
    static int      min_nucleotide_depth;

    // the maximun number of low depth nucleotide 
    //   that a consensus can accept
    static int      max_low_depth;

};

#endif //GLOBALACCESSER_HPP__
