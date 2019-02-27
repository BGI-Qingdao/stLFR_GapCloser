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
#include "GapCloser.hpp"
#include "readhash/Read.hpp"
#include "GlobalAccesser.hpp"

// All threshold here
float    Threshold::NoConflictThreshold ;
int      Threshold::max_allowed_conflict ;
int      Threshold::max_reads_count ;
int      Threshold::min_reads_count ;
int      Threshold::min_sub_reads_count;
int      Threshold::max_error_count ;
int      Threshold::max_reads_depth ;
int      Threshold::the_k ;
// global accessor class pointer here;
ReadAccessor *  GlobalAccesser::the_read_accessor;
PairInfo *      GlobalAccesser::the_pair_info ;

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
    std::cout << "	-a	<string>	input scaffold file name, required." << std::endl;
    std::cout << "	-b	<string>	input library info file name, required." << std::endl;
    std::cout << "	-o	<string>	output file name, required." << std::endl;


    std::cout << "	-l	<int>		maximum read length (<=155), default=" << maxReadLength << ".\n";
    /*	std::cout << "	-m	<int>	overlap mode:" << std::endl;
        std::cout << "		1:	fixed overlap mode" << std::endl;
        std::cout << "		2:	max overlap first mode" << std::endl;
        */
    std::cout << "	-p	<int>		overlap param(<=31), default=25.\n";
    std::cout << "	-t	<int>		thread number, default=1.\n";

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
    Short_Len_t overlapMode=ContigAssembler::fixedOverlapMode;
    Short_Len_t overlapParam=25;
    float loadFactor = 0.75;
    Len_t threadSum=1;

    int c;
    while((c=getopt(argc, argv, "i:o:e:b:a:l:m:p:c:t:N:")) !=-1) {
        switch(c) {
            case 'i': infile=optarg; break;
            case 'o': outfile=optarg; break;
            case 'e': inPairEndInfo=optarg; break;
            case 'b': inLibInfo=optarg; break;

            case 'a': infileContig=optarg; break;

            case 'l': maxReadLength=atoi(optarg); break;
            case 'm': overlapMode=atoi(optarg); break;
            case 'p': overlapParam=atoi(optarg); break;
            case 'c': loadFactor=atof(optarg); break;
            case 't': threadSum=atoi(optarg); break;
            case 'N': NNumber = atoi(optarg); break;

            case 'h': usage(); break;
            case '?': usage(); break;
        }
    }

    if (overlapParam < 13) {
        std::cout << "[WARNING] Overlap length should be >= 13. Program will use 13 instead of " << overlapParam << ".\n";
        overlapParam = 13;
    } else if (overlapParam > 31) {
        std::cout << "[WARNING] Overlap length should be <= 31. Program will use 31 instead of " << overlapParam << ".\n";
        overlapParam = 31;
    }


    endNumLen = overlapParam;

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

    if ((int)maxReadLength > (int)Read::DATA_MAXLEN) {
        std::cout << "[WARNING] Maximum supported read length is 155 bp! Program will use 155 instead of "<< maxReadLength << ".\n";
        maxReadLength = Read::DATA_MAXLEN;
    }

    Read::DATA_MAXLEN= (int)Read::DATA_MAXLEN>(int)maxReadLength?(int)maxReadLength:Read::DATA_MAXLEN;	

    std::cout << "Program: GapCloser" << std::endl;
    std::cout << "Version: 1.12" << std::endl << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "    -a (scaffold file): " << infileContig << std::endl;
    std::cout << "    -b (config file):   " << inLibInfo << std::endl;
    std::cout << "    -o (output file):   " << outfile << std::endl;
    std::cout << "    -l (max read len):  " << (int)Read::DATA_MAXLEN << std::endl;
    std::cout << "    -p (overlap para):  " << (int)overlapParam << std::endl;
    std::cout << "    -t (thread num):    " << (int)threadSum << std::endl << std::endl;


    PairInfo* pairInfo;
    ReadHash* readHash;
    if(!inLibInfo) {

        pairInfo = new PairInfo(&finPairEndInfo);
        readHash = new ReadHash(&fin, NULL, overlapParam, pairInfo->getReadsSum(), hashLen, loadFactor);
    }
    else {

        LibInfo libInfo(&finLibInfo);
        readHash = new ReadHash(NULL, &libInfo, overlapParam, 0, hashLen, loadFactor);
        pairInfo = new PairInfo(NULL, &libInfo);
    }

    ReadAccessor readAccessor(*readHash, *pairInfo);
    ContigTable contigTable(finContig, endNumLen, mismatchLen);

    GapCloser gapcloser(outfile, fout, readAccessor, *pairInfo, contigTable, threadSum, deviation, endNumLen, mismatchLen, maxReadLength, overlapMode, overlapParam);
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


    return 0;
}
