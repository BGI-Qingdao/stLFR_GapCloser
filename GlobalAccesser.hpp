#ifndef GLOBALACCESSER_HPP__
#define GLOBALACCESSER_HPP__

#include <string>
#include "TagId.hpp"
#include "freq.h"
#include "Misc.h"

class ReadAccessor ;
class PairInfo ;

struct GlobalAccesser
{
    static ReadAccessor * the_read_accessor;
    static PairInfo * the_pair_info ;
    static TagId barcode_ider;
    static BGIQD::FREQ::Freq<int> consensus_result_freq ;
    static BGIQD::FREQ::Freq<std::string> consensus_failed_reason ;
    static BGIQD::FREQ::Freq<int> conflict_freq;
    static BGIQD::FREQ::Freq<int> too_low_freq;
    static BGIQD::FREQ::Freq<int> basic_reads_set_freq;
    static BGIQD::FREQ::Freq<int> used_reads_set_freq;
    static BGIQD::FREQ::Freq<std::string> gap_type;
    static BGIQD::FREQ::Freq<Sub1ReadNum> sub1_read_num ;
    static BGIQD::FREQ::Freq<Sub1_3ReadNum> sub1_3_read_num ;
};

struct Threshold
{
    // the minumum proportion that a nucleotide should
    //  hold to be a major nucletide.
    static float    NoConflictThreshold ;

    // the maximum number of confilicts that a consensus.
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
    //static int      middle_sub_reads_count;

    // the maximum number of confilicts that a map.
    //  can accept.
    static int      max_error_count ;

    // the maxinum number of depth that a read can be accept.
    static int      max_reads_depth ;

    // the kvalue for kmer-read map.
    static int      the_k ;

    // the minimum depth of a nucleotide that 
    //  make it a NOT to be a low depth nucleotide.
    static int      min_nucleotide_depth;

    // the maximun number of low depth nucleotide 
    //   that a consensus can accept.
    static int      max_accept_low_depth;

    // the maximum bp of a gap that can define as
    //   a small gap.
    static int      max_small_gap;
    //static int      max_middle_gap;

    // how many bp of N will be insert into a unfinished gap 
    static int      NNumber;

    // the maximum length of a read that can loaded .
    static int      maxReadLength;

    // 1 means DO NOT FILL GAP which is smaller than extend_len
    static int      filter_too_small_gap ;
};

#endif //GLOBALACCESSER_HPP__
