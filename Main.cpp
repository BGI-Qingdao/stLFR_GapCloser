/**************************************************
 *
 * Main.cpp
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
 *
 * This file is part of GapCloser.
 *
 * GapCloser is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GapCloser is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GapCloser. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************/

//#define NO_REVERSE

#include "Common.hpp"
#include "Utils.hpp"
#include "GlobalAccesser.hpp"
#include "ConsensusConfig.hpp"
#include "readhash/Read.hpp"
#include "GapCloser.hpp"

/*****************************************
 *
 *  All threshold here
 *
 * ***************************************/

// The K Value 
int      Threshold::the_k = 27 ;

// For map read to contig 
int      Threshold::max_error_count = 2 ;
int      Threshold::max_reads_depth = 1000 ;

// For consensus area 
int      ConsensusConfig::extend_len = 10 ;
int      ConsensusConfig::prev_extra_len = 90 ;
int      ConsensusConfig::last_extra_len = 10 ;
int      ConsensusConfig::consensus_len = Threshold::the_k + 10 + ConsensusConfig::extend_len ;

// For reads set
int      Threshold::max_reads_count = 300 ;
int      Threshold::min_reads_count = 1 ;
int      Threshold::max_kmer_2_read = 10 ;

// For sub read set
//int      Threshold::max_small_gap = 10 ;
//int      Threshold::max_middle_gap = 1000 ;
int      Threshold::min_pe_sub_reads_count = 10 ;
int      Threshold::min_pe_barcode_sub_reads_count = 10 ;
//int      Threshold::middle_sub_reads_count = 10 ;

// For consensus
float    Threshold::NoConflictThreshold = 0.8f ;
float    Threshold::basic_NoConflictThreshold = 0.8f ;
int      Threshold::max_allowed_conflict = 2 ;
int      Threshold::min_nucleotide_depth = 2 ;
int      Threshold::max_accept_low_depth = 2 ;

// For gap fill
int      Threshold::NNumber = 1;
int      Threshold::maxReadLength = 150;
//int      Threshold::filter_too_small_gap = 1;
int      Threshold::use_subset_only = 0 ;

int      Threshold::basic_set_max_conflict = 1;
int      Threshold::basic_set_max_low_depth = 1;

// global accessor class pointer here;
ReadAccessor *  GlobalAccesser::the_read_accessor = NULL;
PairInfo *      GlobalAccesser::the_pair_info = NULL;
TagId           GlobalAccesser::barcode_ider ;

BGIQD::FREQ::Freq<int> GlobalAccesser::consensus_result_freq ;
BGIQD::FREQ::Freq<std::string> GlobalAccesser::consensus_failed_reason ;
BGIQD::FREQ::Freq<int> GlobalAccesser::conflict_freq;
BGIQD::FREQ::Freq<int> GlobalAccesser::too_low_freq;
BGIQD::FREQ::Freq<int> GlobalAccesser::basic_reads_set_freq;
BGIQD::FREQ::Freq<int> GlobalAccesser::used_reads_set_freq;
BGIQD::FREQ::Freq<int> GlobalAccesser::kmer_read_count;

BGIQD::FREQ::Freq<std::string> GlobalAccesser::sub_type;
BGIQD::FREQ::Freq<SubReadsLog> GlobalAccesser::sub_read_num ;

//the ReadAccesser will assign this
Read          * ReadElement::read=NULL;

void usage(void)
{
    std::cout << "Version:" << std::endl;
    std::cout << "	1.12" << std::endl;
    std::cout << std::endl;

    std::cout << "Contact:" << std::endl;
    std::cout << "	soap@genomics.org.cn" << std::endl;
    std::cout << std::endl;

    std::cout << "Usage:" << std::endl;
    std::cout << "	GapCloser [options]" << std::endl;
    std::cout << "Basic options :\n";
    std::cout << "	-a	<string>	input scaffold file name, required." << std::endl;
    std::cout << "	-b	<string>	input library info file name, required." << std::endl;
    std::cout << "	-o	<string>	output file name, required." << std::endl;
    std::cout << "	-l	<int>		maximum read length (<=155), default=" 
        << Threshold::maxReadLength << ".\n";
    std::cout << "	-p	<int>		overlap param(<=31) [the kvalue], default=25.\n";
    //std::cout << "	-N	<int>		how many N that insert into a unfinished gap .\n";

    std::cout << "Performance options :\n";
    std::cout << "	-t	<int>		thread number, default=1.\n";
    std::cout << "	-c	<float>		hash load fractor, default=0.75.\n";


    std::cout << " ---------- new parameters below : ---------\n";
    std::cout << "mapping read to contig options:\n";
    std::cout << "	-1	<int>		maximum read depth, default="
        << Threshold::max_reads_depth << ".\n";
    std::cout << "	-2	<int>		maximum mismatch , default=" 
        << Threshold::max_error_count<< ".\n";

    std::cout << "consensus reads set options:\n";
    std::cout << "	-3	<int>		consensus length, default= [ the kvalue ] + 10 + [consensus extend length] ="
        << ConsensusConfig::consensus_len<< ".\n";
    std::cout << "	-4	<int>		consensus prev extra length, default=" 
        << ConsensusConfig::prev_extra_len<< ".\n";
    std::cout << "	-5	<int>		consensus last extra length, default=" 
        << ConsensusConfig::last_extra_len<< ".\n";
    std::cout << "	-6	<int>		consensus extend length, default=" 
        << ConsensusConfig::extend_len<< ".\n";

    std::cout << "consensus reads set options:\n";
    std::cout << "	-7	<int>		min reads count , default=" 
        << Threshold::min_reads_count<< ".\n";
    std::cout << "	-8	<int>		max reads count , default=" 
        << Threshold::max_reads_count<< ".\n";

    std::cout << "extract sub reads set options:\n";
    std::cout << "	-A	<int>		min pe subset reads count threshold , default=" 
        << Threshold::min_pe_sub_reads_count<< ".\n";
    std::cout << "	-B	<int>		min pe & barcode subset reads count threshold , default=" 
        << Threshold::min_pe_barcode_sub_reads_count<< ".\n";

    std::cout << "consensus options:\n";
    std::cout << "	-C	<float>		non-conflict threshold , default=" 
        << Threshold::NoConflictThreshold<< ".\n";
    std::cout << "	-D	<int>		max number of conflict can accepted , default=" 
        << Threshold::max_allowed_conflict<< ".\n";
    std::cout << "	-E	<int>		min depth threshold, default=" 
        << Threshold::min_nucleotide_depth<< ".\n";
    std::cout << "	-F	<int>		max number of low depth nucleotide, default=" 
        << Threshold::max_accept_low_depth<< ".\n";

    std::cout << "other:\n";
    std::cout << "	-G	<int>		the max number of reads that a kmer can find . default=" 
        << Threshold::max_kmer_2_read<< ".\n";
    std::cout << "	-I	<int>		use sub-reads-set only. default=" 
        << Threshold::use_subset_only<< ".\n";
    std::cout << "	-O	<int>		[ basic set ] max number of conflict can accepted . default=" 
        << Threshold::basic_set_max_conflict<< ".\n";
    std::cout << "	-P	<int>	 	[ basic set ] max number of low depth nucleotide . default=" 
        << Threshold::basic_set_max_low_depth<< ".\n";
    std::cout << "	-Q	<float>	 	[ basic set ] the conflict_threshold of baisc set . default=" 
        << Threshold::basic_NoConflictThreshold<< ".\n";
    std::cout << " ---------- new parameters end     ---------\n";

    std::cout << "	-h	-?		output help information." << std::endl;
    std::cout << std::endl;

    exit(1);
}

Len_t Read::DATA_MAXLEN = (Read::DATA_ARRAY_SIZE*bitsizeof(Number_t) + Read::BITLEN_DATA_REMAIN) / 2 ;

int main(int argc, char *argv[])
{
    //input parameters
    char *infile=NULL;
    char *outfile=NULL;
    char *inPairEndInfo=NULL;
    char *inLibInfo=NULL;
    char *infileContig=NULL;
    float deviation=0.5;
    Len_t endNumLen=10;
    Len_t mismatchLen=5;
    float loadFactor = 0.75;
    Len_t threadSum=1;
    bool consensus_len_setted_flag = false ;
    int c;
    while((c=getopt(argc, argv, 
                    "i:o:e:b:a:l:p:c:t:N:1:2:3:4:5:6:7:8:9:A:B:C:D:E:F:G:I:O:P:Q:"))
            !=-1) 
    {
        switch(c) {
            case 'a': infileContig=optarg; break;
            case 'b': inLibInfo=optarg; break;
            case 'o': outfile=optarg; break;
            case 'l': Threshold::maxReadLength=atoi(optarg); break;
            case 'p': Threshold::the_k=atoi(optarg); break;
            case 't': threadSum=atoi(optarg); break;
            case 'c': loadFactor=atof(optarg); break;
            case 'N': Threshold::NNumber = atoi(optarg); break;

            case 'e': inPairEndInfo=optarg; break; /* why */
            case 'i': infile=optarg; break; /* why */

            /*case 'm': overlapMode=atoi(optarg); break;*/

            case '1': Threshold::max_reads_depth= atoi(optarg); break;
            case '2': Threshold::max_error_count= atoi(optarg); break;

            case '3': {
                          ConsensusConfig::consensus_len= atoi(optarg); 
                          consensus_len_setted_flag = true ;
                          break;
                      }
            case '4': ConsensusConfig::prev_extra_len = atoi(optarg); break;
            case '5': ConsensusConfig::last_extra_len = atoi(optarg); break;
            case '6': ConsensusConfig::extend_len = atoi(optarg); break;

            case '7': Threshold::min_reads_count = atoi(optarg); break;
            case '8': Threshold::max_reads_count = atoi(optarg); break;

            case 'A': Threshold::min_pe_sub_reads_count = atoi(optarg); break;
            case 'B': Threshold::min_pe_barcode_sub_reads_count = atoi(optarg); break;

            case 'C': Threshold::NoConflictThreshold = atof(optarg); break;
            case 'D': Threshold::max_allowed_conflict = atoi(optarg); break;
            case 'E': Threshold::min_nucleotide_depth = atoi(optarg); break;
            case 'F': Threshold::max_accept_low_depth = atoi(optarg); break;
            case 'G': Threshold::max_kmer_2_read= atoi(optarg); break;
            case 'I': Threshold::use_subset_only = atoi(optarg); break;
            case 'O': Threshold::basic_set_max_conflict = atoi(optarg); break;
            case 'P': Threshold::basic_set_max_low_depth = atoi(optarg); break;
            case 'Q': Threshold::basic_NoConflictThreshold= atof(optarg); break;

            case 'h': usage(); break;
            case '?': usage(); break;
        }
    }

    if (Threshold::the_k< 13) {
        std::cout << "[WARNING] Overlap length should be >= 13. Program will use 13 instead of " << Threshold::the_k << ".\n";
        Threshold::the_k = 13;
    } else if (Threshold::the_k  > 31) {
        std::cout << "[WARNING] Overlap length should be <= 31. Program will use 31 instead of " << Threshold::the_k << ".\n";
        Threshold::the_k = 31;
    }
    if( ! consensus_len_setted_flag )
    {
        ConsensusConfig::consensus_len = Threshold::the_k + 10 + ConsensusConfig::extend_len;
    }

    endNumLen = Threshold::the_k;

    Number_t hashLen = 3;

    //check input files
    std::ifstream fin;
    std::ifstream finPairEndInfo;
    std::ifstream finLibInfo;
    if(!inLibInfo) {

        fin.open(infile);
        if(!fin) {
            std::cout << "[Error] Can not open input file." << std::endl << std::endl;
            usage();
        }

        finPairEndInfo.open(inPairEndInfo);
        if(!finPairEndInfo) {
            std::cout << "[Error] Can not open input pair-end info file." << std::endl << std::endl;
            usage();
        }
    }
    else {

        finLibInfo.open(inLibInfo);
        if(!finLibInfo) {
            std::cout << "[Error] Can not open input library info file." << std::endl << std::endl;
            usage();
        }
    }

    std::ofstream fout(outfile);
    if(!fout) {
        std::cout << "[Error] Can not creat output file." << std::endl << std::endl;
        usage();
    }

    std::ifstream finContig(infileContig);
    if(!finContig) {
        std::cout << "[Error] Can not open input scaffold file." << std::endl << std::endl;
        usage();
    }

    if ((int)Threshold::maxReadLength > (int)Read::DATA_MAXLEN) {
        std::cout << "[WARNING] Maximum supported read length is 155 bp! Program will use 155 instead of "<< Threshold::maxReadLength << ".\n";
        Threshold::maxReadLength = Read::DATA_MAXLEN;
    }

    Read::DATA_MAXLEN= (int)Read::DATA_MAXLEN>(int)Threshold::maxReadLength?(int)Threshold::maxReadLength:Read::DATA_MAXLEN;	

    std::cout << "Program: GapCloser" << std::endl;
    std::cout << "Version: 1.12" << std::endl << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "    -a (scaffold file): " << infileContig << std::endl;
    std::cout << "    -b (config file):   " << inLibInfo << std::endl;
    std::cout << "    -o (output file):   " << outfile << std::endl;
    std::cout << "    -l (max read len):  " << (int)Read::DATA_MAXLEN << std::endl;
    std::cout << "    -p (overlap para):  " << (int)Threshold::the_k<< std::endl;
    std::cout << "    -t (thread num):    " << (int)threadSum << std::endl << std::endl;
    std::cout << "    -c (map loadFactor):" << (int)loadFactor<< std::endl << std::endl;

    std::cout << " ---------- new parameters below : ---------\n";
    std::cout << "mapping read to contig options:\n";
    std::cout << "	-1	<int>		maximum read depth,  ="
        << Threshold::max_reads_depth << ".\n";
    std::cout << "	-2	<int>		maximum mismatch ,  =" 
        << Threshold::max_error_count<< ".\n";

    std::cout << "consensus reads set options:\n";
    std::cout << "	-3	<int>		consensus length,  =" 
        << ConsensusConfig::consensus_len<< ".\n";
    std::cout << "	-4	<int>		consensus prev extra length,  =" 
        << ConsensusConfig::prev_extra_len<< ".\n";
    std::cout << "	-5	<int>		consensus last extra length,  =" 
        << ConsensusConfig::last_extra_len<< ".\n";
    std::cout << "	-6	<int>		consensus extend length,  =" 
        << ConsensusConfig::extend_len<< ".\n";

    std::cout << "consensus reads set options:\n";
    std::cout << "	-7	<int>		min reads count ,  =" 
        << Threshold::min_reads_count<< ".\n";
    std::cout << "	-8	<int>		max reads count ,  =" 
        << Threshold::max_reads_count<< ".\n";

    std::cout << "abstract sub reads set options:\n";
    std::cout << "	-A	<int>		min pe subset reads count threshold , =" 
        << Threshold::min_pe_sub_reads_count<< ".\n";
    std::cout << "	-B	<int>		min pe & barcode subset reads count threshold , =" 
        << Threshold::min_pe_barcode_sub_reads_count<< ".\n";

    std::cout << "consensus options:\n";
    std::cout << "	-C	<float>		non-conflict threshold ,  =" 
        << Threshold::NoConflictThreshold<< ".\n";
    std::cout << "	-D	<int>		max conflict can accepted ,  =" 
        << Threshold::max_allowed_conflict<< ".\n";
    std::cout << "	-E	<int>		min depth threshold,  =" 
        << Threshold::min_nucleotide_depth<< ".\n";
    std::cout << "	-F	<int>		max low depth nucleotide,  =" 
        << Threshold::max_accept_low_depth<< ".\n";
    std::cout << "other:\n";
    std::cout << "	-G	<int>		the max number of reads that a kmer can find . =" 
        << Threshold::max_kmer_2_read<< ".\n";
    std::cout << "	-I	<int>		use sub-reads-set only. =" 
        << Threshold::use_subset_only<< ".\n";
    std::cout << "	-O	<int>		[ basic set ] max number of conflict can accepted . =" 
        << Threshold::basic_set_max_conflict<< ".\n";
    std::cout << "	-P	<int>	 	[ basic set ] max number of low depth nucleotide . =" 
        << Threshold::basic_set_max_low_depth<< ".\n";
    std::cout << "	-Q	<float>	 	[ basic set ] the conflict_threshold of baisc set . =" 
        << Threshold::basic_NoConflictThreshold<< ".\n";
    std::cout << " ---------- new parameters end     ---------\n";

    PairInfo* pairInfo;
    ReadHash* readHash;
    if(!inLibInfo) {
        pairInfo = new PairInfo(&finPairEndInfo);
        readHash = new ReadHash(&fin, NULL, Threshold::the_k , pairInfo->getReadsSum(), hashLen, loadFactor);
    }
    else {
        LibInfo libInfo(&finLibInfo);
        readHash = new ReadHash(NULL, &libInfo, Threshold::the_k , 0, hashLen, loadFactor);
        pairInfo = new PairInfo(NULL, &libInfo);
    }
    GlobalAccesser::the_pair_info = pairInfo ;
    ReadAccessor readAccessor(*readHash, *pairInfo);
    GlobalAccesser::the_read_accessor = &readAccessor ;
    ContigTable contigTable(finContig, endNumLen, mismatchLen);

    GapCloser gapcloser( outfile, fout, readAccessor,
            *pairInfo, contigTable, threadSum,
            deviation, endNumLen, mismatchLen
           /* , maxReadLength */
           /* , overlapMode */
            , Threshold::the_k );
    gapcloser.assemble();

    delete pairInfo;
    delete readHash;

    if(!inLibInfo) {
        fin.close();
        finPairEndInfo.close();
    }
    else {
        finLibInfo.close();
    }
    fout.close();
    finContig.close();

    std::cerr<<" consensus result freq         :\n"<<GlobalAccesser::consensus_result_freq.ToString() ;
    std::cerr<<" consensus failed reason  freq :\n"<<GlobalAccesser::consensus_failed_reason.ToString() ;
    std::cerr<<" size of base reads set freq        :\n"<<GlobalAccesser::basic_reads_set_freq.ToString() ;
    std::cerr<<" size of used reads set freq        :\n"<<GlobalAccesser::used_reads_set_freq.ToString() ;
    std::cerr<<" conflict freq                 :\n"<<GlobalAccesser::conflict_freq.ToString() ;
    std::cerr<<" too low depth freq            :\n"<<GlobalAccesser::too_low_freq.ToString() ;
    std::cerr<<" sub read set type freq        :\n"<<GlobalAccesser::sub_type.ToString() ;
    std::cerr<<" detail of read set freq        :\n"<<GlobalAccesser::sub_read_num.ToString() ;

    std::cerr<<" detail of kmer read count        :\n"<<GlobalAccesser::kmer_read_count.ToString() ;
    return 0;
}
