/**************************************************
 *
 * ReadHash.hpp
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

#ifndef READHASH_HPP_
#define READHASH_HPP_

#include <string>
#include <fstream>
#include "Utils.hpp"
#include "GlobalAccesser.hpp"
#include "Read.hpp"
#include "HashTable.hpp"
#include "QuickSorterMulti.hpp"
#include "readtool/LibInfo.hpp"
#include "readtool/stLFRReadHeader.hpp"


int GetBarcodeFromReadName( const char* name)
{
    std::string read_name(name);
    stLFRHeader header;
    header.Init(name);
    if( header.type == stLFRHeader::ReadType::Unknow )
        return 0 ;
    else if ( header.type
            == stLFRHeader::ReadType::readName_barcodeStr
            ||
            header.type 
            == stLFRHeader::ReadType::readName_barcodeStr_index
            ||
            header.type 
            == stLFRHeader::ReadType::readName_barcodeStr_index_barcodeNum 
            )
    {
        return GlobalAccesser::barcode_ider.GetId(header.barcode_str);
    }
    else 
        return 0 ;
}


static const Number_t SEED_NUM = 1073741824;


static const unsigned int COUNT_BITLEN = 30;
static const unsigned int DEPTH_BITLEN = 30;
static const unsigned int SERIAL_NUM_BITLEN = 34;


//struct SortElement
//{
//	Number_t serialNum : SERIAL_NUM_BITLEN;
//	Number_t depth : DEPTH_BITLEN;
//};

//struct SerialNumElement
//{
//	Number_t start : SERIAL_NUM_BITLEN;
//	Number_t count : COUNT_BITLEN;
//	
//	SerialNumElement() : 
//		start(0), 
//		count(0)
//	{
//	}
//};

struct SerialNumElement
{
    Len_t start;  //start position

    SerialNumElement() : 
        start(~0)
    {
    }
};

struct ReadsHashElement
{
    SerialNumElement serialNumForward;  //forward start position in read array
    SerialNumElement serialNumReverse;
};


class ReadHash
{
    protected:

        std::ifstream* finp;  //single input file
        LibInfo* pLibInfo;  // libraries information input

        Read* reads;
        Number_t readsSum;

        Len_t* serialNumsForward;
        Len_t* serialNumsReverse;

        HashTable< Number_t, ReadsHashElement > readsHash;

        Short_Len_t overlapParam;

    public:

        ReadHash(std::ifstream* _finp, LibInfo* _pLibInfo=NULL, Short_Len_t _overlapParam=25, Number_t _readsSum=0, Number_t _hashLen=SEED_NUM, float _loadFactor=0.75) : 
            finp(_finp), 
            pLibInfo(_pLibInfo), 
            reads(NULL), 
            readsSum(_readsSum), 
            serialNumsForward(NULL), 
            serialNumsReverse(NULL), 
            readsHash(_hashLen, _loadFactor), 
            overlapParam(_overlapParam)
    {
        initialize();
    }

        virtual ~ReadHash()
        {	
            freeMemorys();
        }

        Read* getReads() const { return reads; }
        Number_t getReadsSum() const { return readsSum; }

        Len_t* getSerialNumsForward() const { return serialNumsForward; }
        Len_t* getSerialNumsReverse() const { return serialNumsReverse; }

        HashTable< Number_t, ReadsHashElement > const& getHash() const { return readsHash; }
        Number_t getKmersSum() const { return readsHash.getCount(); }

        Short_Len_t getOverlapParam() const { return overlapParam; }

        void initialize() {

            std::cout << std::endl;
            std::cout << ">>>>>>>>>>hash initializing<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
            time_t total_start_time=time(NULL);



            time_t start_time;

            std::cout << "counting reads" << std::endl;
            start_time=time(NULL);

            countReads();
            std::cout << std::endl;
            std::cout << "reads sum: " << readsSum <<std::endl;
            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads counted" << std::endl;
            std::cout <<std::endl;

            if (readsSum > 4294967296)
            {
                std::cout << "At most 4294967296 reads are supported.\nPlease use less reads and try again.\n";
                exit (1);
            }

            std::cout << "allocating memorys" << std::endl;
            start_time=time(NULL);

            allocateMemorys();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "memorys allocated" << std::endl;
            std::cout << std::endl;



            std::cout << "initializing reads" << std::endl;
            start_time=time(NULL);

            initializeReads();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads initialized" << std::endl;
            std::cout << std::endl;



            std::cout << "sorting reads in reversed direction" << std::endl;
            start_time=time(NULL);

            sortReadsReverse();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads sorted" << std::endl;
            std::cout << std::endl;



            std::cout << "initializing hash on reversed reads" << std::endl;
            start_time=time(NULL);

            initializeReadsHashReverse();

            std::cout << "kmers sum: " << getKmersSum() <<std::endl;
            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "hash initialized" << std::endl;
            std::cout << std::endl;



            std::cout << "reversing reads" << std::endl;
            start_time=time(NULL);

            reverseReads();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads reversed" << std::endl;
            std::cout << std::endl;



            std::cout << "sorting reads in forward direction" << std::endl;
            start_time=time(NULL);

            sortReadsForward();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads sorted" << std::endl;
            std::cout << std::endl;



            std::cout << "initializing hash on forward reads" << std::endl;
            start_time=time(NULL);

            initializeReadsHashForward();

            std::cout << "kmers sum: " << getKmersSum() <<std::endl;
            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "hash initialized" << std::endl;
            std::cout << std::endl;



            std::cout << "marking unique reads" << std::endl;
            start_time=time(NULL);

            uniqueReads();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "unique reads marked" << std::endl;
            std::cout << std::endl;



            std::cout << "sorting reads by serial number" << std::endl;
            start_time=time(NULL);

            // recover the read order to their input order with all reads in the strand of smaller one
            sortReadsBySerialNum();

            std::cout << "spent time: " << time(NULL)-start_time << "s" <<std::endl;
            std::cout << "reads sorted" << std::endl;
            std::cout << std::endl;



            std::cout << "spent total time: " << time(NULL)-total_start_time << "s" <<std::endl;
            std::cout << ">>>>>>>>>>hash initialization finished<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
        }

    private:

        void freeMemorys() {

            delete [] reads;
            delete [] serialNumsForward;
            delete [] serialNumsReverse;
        }


        void countReads() {

            if (pLibInfo) {

                int i,j;
                std::ifstream _fin;
                for (i=0; i<(int)pLibInfo->getArrayLen(); i++) {

                    LibInfoElement& libInfo = pLibInfo->getArray()[i];

                    std::cout << "in lib: " << libInfo.name << std::endl;

                    if (libInfo.readsSum==0 || libInfo.pairReadsSum==0) {

                        if (libInfo.readsSum==0 && libInfo.pairReadsSum==0) {

                            for(j=0;j<libInfo.num_p_file;j++) {
                                if (libInfo.p_fname[j].substr(libInfo.p_fname[j].size()-3) != ".gz"){
                                    _fin.open(libInfo.p_fname[j].c_str());
                                    if(!_fin) {
                                        std::cout << std::endl << "Error: can not open " << libInfo.p_fname[j] << std::endl;
                                        exit(1);
                                    }
                                    libInfo.pairReadsSum += countReadsFa(_fin);
                                    _fin.close();
                                    _fin.clear();
                                }
                                else{
                                    libInfo.pairReadsSum += countReadsFaGZ(libInfo.p_fname[j]);
                                }
                            }

                            for(j=0;j<libInfo.num_a1_file;j++) {
                                if (libInfo.a1_fname[j].substr(libInfo.a1_fname[j].size()-3) != ".gz"){
                                    _fin.open(libInfo.a1_fname[j].c_str());
                                    if(!_fin) {
                                        std::cout << std::endl << "Error: can not open " << libInfo.a1_fname[j] << std::endl;
                                        exit(1);
                                    }
                                    libInfo.pairReadsSum += countReadsFa(_fin);
                                    _fin.close();
                                    _fin.clear();
                                }
                                else{
                                    libInfo.pairReadsSum += countReadsFaGZ(libInfo.a1_fname[j]);
                                }
                            }

                            for(j=0;j<libInfo.num_a2_file;j++) {
                                if (libInfo.a2_fname[j].substr(libInfo.a2_fname[j].size()-3) != ".gz"){
                                    _fin.open(libInfo.a2_fname[j].c_str());
                                    if(!_fin) {
                                        std::cout << std::endl << "Error: can not open " << libInfo.a2_fname[j] << std::endl;
                                        exit(1);
                                    }
                                    libInfo.pairReadsSum += countReadsFa(_fin);
                                    _fin.close();
                                    _fin.clear();
                                }
                                else{
                                    libInfo.pairReadsSum += countReadsFaGZ(libInfo.a2_fname[j]);
                                }
                            }

                            for(j=0;j<libInfo.num_q1_file;j++) {
                                if (libInfo.q1_fname[j].substr(libInfo.q1_fname[j].size()-3) != ".gz"){
                                    _fin.open(libInfo.q1_fname[j].c_str());
                                    if(!_fin) {
                                        std::cout << std::endl << "Error: can not open " << libInfo.q1_fname[j] << std::endl;
                                        exit(1);
                                    }
                                    libInfo.pairReadsSum += countReadsFq(_fin);
                                    _fin.close();
                                    _fin.clear();
                                }
                                else{
                                    libInfo.pairReadsSum += countReadsFqGZ(libInfo.q1_fname[j]);
                                }
                            }

                            for(j=0;j<libInfo.num_q2_file;j++) {
                                if (libInfo.q2_fname[j].substr(libInfo.q2_fname[j].size()-3) != ".gz"){
                                    _fin.open(libInfo.q2_fname[j].c_str());
                                    if(!_fin) {
                                        std::cout << std::endl << "Error: can not open " << libInfo.q2_fname[j] << std::endl;
                                        exit(1);
                                    }
                                    libInfo.pairReadsSum += countReadsFq(_fin);
                                    _fin.close();
                                    _fin.clear();
                                }
                                else{
                                    libInfo.pairReadsSum += countReadsFqGZ(libInfo.q2_fname[j]);
                                }
                            }
                        }

                        Number_t singleReadsSum=0;
                        for(j=0;j<libInfo.num_s_a_file;j++) {
                            if (libInfo.s_a_fname[j].substr(libInfo.s_a_fname[j].size()-3) != ".gz"){
                                _fin.open(libInfo.s_a_fname[j].c_str());
                                if(!_fin) {
                                    std::cout << std::endl << "Error: can not open " << libInfo.s_a_fname[j] << std::endl;
                                    exit(1);
                                }
                                singleReadsSum += countReadsFa(_fin);
                                _fin.close();
                                _fin.clear();
                            }
                            else{
                                singleReadsSum += countReadsFaGZ(libInfo.s_a_fname[j]);
                            }
                        }

                        for(j=0;j<libInfo.num_s_q_file;j++) {
                            if (libInfo.s_q_fname[j].substr(libInfo.s_q_fname[j].size()-3) != ".gz"){
                                _fin.open(libInfo.s_q_fname[j].c_str());
                                if(!_fin) {
                                    std::cout << std::endl << "Error: can not open " << libInfo.s_q_fname[j] << std::endl;
                                    exit(1);
                                }
                                singleReadsSum += countReadsFq(_fin);
                                _fin.close();
                                _fin.clear();
                            }
                            else{
                                singleReadsSum += countReadsFqGZ(libInfo.s_q_fname[j]);
                            }
                        }

                        if (libInfo.readsSum==0)
                            libInfo.readsSum = libInfo.pairReadsSum + singleReadsSum;
                        if (libInfo.pairReadsSum==0)
                            libInfo.pairReadsSum = libInfo.readsSum - singleReadsSum;
                    }

                    readsSum += libInfo.readsSum;

                    std::cout << "    total reads sum:     " << libInfo.readsSum << std::endl;
                    std::cout <<"    Paired-end reads sum: " << libInfo.pairReadsSum <<std::endl;
                }
            }
            else {

                if (readsSum == 0) {
                    readsSum += countReadsFa(*finp);
                }
            }
        }

        Number_t countReadsFa(std::ifstream& fin) {

            Number_t readsSum=0;
            char str[MAX_STRING_LEN];
            while(!fin.eof()) {

                fin.getline(str, MAX_STRING_LEN);
                if (str[0] == '>') {
                    continue;
                }
                else {
                    chomp(str);
                    size_t size = strlen(str);
                    if (size == 0)
                        continue;
                    readsSum++;
                }
            }
            return readsSum;
        }

        Number_t countReadsFaGZ(std::string& fin) {

            char *cmd = new char[fin.size()+20];
            sprintf(cmd, "gzip -dc %s", fin.c_str());
            FILE *fp = popen(cmd, "r");

            Number_t readsSum=0;
            char str[MAX_STRING_LEN];
            while(fgets(str, MAX_STRING_LEN, fp)) {

                if (str[0] == '>') {
                    continue;
                }
                else {
                    chomp(str);
                    size_t size = strlen(str);
                    if (size == 0)
                        continue;
                    readsSum++;
                }
            }

            delete [] cmd;
            pclose(fp);

            return readsSum;
        }

        Number_t countReadsFq(std::ifstream& fin) {

            Number_t readsSum=0;
            char name[MAX_STRING_LEN];
            char sequence[MAX_STRING_LEN];
            char quality[MAX_STRING_LEN];
            while(!fin.eof()) {

                fin.getline(name, MAX_STRING_LEN);
                fin.getline(sequence, MAX_STRING_LEN);
                fin.getline(name, MAX_STRING_LEN);
                fin.getline(quality, MAX_STRING_LEN);

                chomp(name);
                chomp(sequence);
                chomp(quality);

                size_t size = strlen(sequence);
                if (size == 0)
                    continue;

                readsSum++;
            }
            return readsSum;
        }

        Number_t countReadsFqGZ(std::string& fin) {

            char *cmd = new char[fin.size()+20];
            sprintf(cmd, "gzip -dc %s", fin.c_str());
            FILE *fp = popen(cmd, "r");

            Number_t readsSum=0;
            char name[MAX_STRING_LEN];
            char sequence[MAX_STRING_LEN];
            char quality[MAX_STRING_LEN];
            while(fgets(name, MAX_STRING_LEN, fp)) {

                fgets(sequence, MAX_STRING_LEN, fp);
                fgets(name, MAX_STRING_LEN, fp);
                fgets(quality, MAX_STRING_LEN, fp);


                chomp(name);
                chomp(sequence);
                chomp(quality);

                size_t size = strlen(sequence);
                if (size == 0)
                    continue;

                readsSum++;
            }

            delete [] cmd;
            pclose(fp);

            return readsSum;
        }


        void allocateMemorys() {

            reads = new Read[readsSum];
            serialNumsForward = new Len_t[readsSum];
            serialNumsReverse = new Len_t[readsSum];
        }

        void initializeReads() {

            //for index use
            readsSum = 0;

            if (pLibInfo) {

                int i,j;
                std::ifstream _fin,_fin1,_fin2;
                for (i=0; i<(int)pLibInfo->getArrayLen(); i++) {

                    LibInfoElement& libInfo = pLibInfo->getArray()[i];

                    std::cout << "in lib: " << libInfo.name << std::endl;

                    // type 3, paired-end fasta reads in one file
                    for(j=0;j<libInfo.num_p_file;j++) {
                        if (libInfo.p_fname[j].substr(libInfo.p_fname[j].size()-3) != ".gz"){
                            _fin.open(libInfo.p_fname[j].c_str());
                            if(!_fin) {
                                std::cout << std::endl << "Error: can not open " << libInfo.p_fname[j] << std::endl;
                                exit(1);
                            }
                            initializeReadsFa(_fin);
                            _fin.close();
                            _fin.clear();
                        }
                        else{
                            initializeReadsFaGZ(libInfo.p_fname[j]);
                        }
                    }

                    // type 1, paired-end fasta reads in two files
                    for(j=0;j<libInfo.num_a1_file;j++) {
                        if (libInfo.a1_fname[j].substr(libInfo.a1_fname[j].size()-3) != ".gz"){
                            _fin1.open(libInfo.a1_fname[j].c_str());
                            if(!_fin1) {
                                std::cout << std::endl << "Error: can not open " << libInfo.a1_fname[j] << std::endl;
                                exit(1);
                            }
                            _fin2.open(libInfo.a2_fname[j].c_str());
                            if(!_fin2) {
                                std::cout << std::endl << "Error: can not open " << libInfo.a2_fname[j] << std::endl;
                                exit(1);
                            }
                            initializeReadsFaPair(_fin1, _fin2);
                            _fin1.close();
                            _fin1.clear();
                            _fin2.close();
                            _fin2.clear();
                        }
                        else{
                            initializeReadsFaPairGZ(libInfo.a1_fname[j], libInfo.a2_fname[j]);
                        }
                    }

                    // type 2, paired-end fastq reads in two files
                    for(j=0;j<libInfo.num_q1_file;j++) {
                        if (libInfo.q1_fname[j].substr(libInfo.q1_fname[j].size()-3) != ".gz"){
                            _fin1.open(libInfo.q1_fname[j].c_str());
                            if(!_fin1) {
                                std::cout << std::endl << "Error: can not open " << libInfo.q1_fname[j] << std::endl;
                                exit(1);
                            }
                            _fin2.open(libInfo.q2_fname[j].c_str());
                            if(!_fin2) {
                                std::cout << std::endl << "Error: can not open " << libInfo.q2_fname[j] << std::endl;
                                exit(1);
                            }
                            initializeReadsFqPair(_fin1, _fin2
                                    , libInfo.with_barcode == 1 );
                            _fin1.close();
                            _fin1.clear();
                            _fin2.close();
                            _fin2.clear();
                        }
                        else{
                            initializeReadsFqPairGZ(libInfo.q1_fname[j]
                                    , libInfo.q2_fname[j]
                                    , libInfo.with_barcode == 1 );
                        }
                    }
                }

                for (i=0; i<(int)pLibInfo->getArrayLen(); i++) {

                    LibInfoElement& libInfo = pLibInfo->getArray()[i];

                    // type 4, single fasta read
                    for(j=0;j<libInfo.num_s_a_file;j++) {
                        if (libInfo.s_a_fname[j].substr(libInfo.s_a_fname[j].size()-3) != ".gz"){
                            _fin.open(libInfo.s_a_fname[j].c_str());
                            if(!_fin) {
                                std::cout << std::endl << "Error: can not open " << libInfo.s_a_fname[j] << std::endl;
                                exit(1);
                            }
                            initializeReadsFa(_fin);
                            _fin.close();
                            _fin.clear();
                        }
                        else{
                            initializeReadsFaGZ(libInfo.s_a_fname[j]);
                        }
                    }

                    // type 5, single fastq read
                    for(j=0;j<libInfo.num_s_q_file;j++) {
                        if (libInfo.s_q_fname[j].substr(libInfo.s_q_fname[j].size()-3) != ".gz"){
                            _fin.open(libInfo.s_q_fname[j].c_str());
                            if(!_fin) {
                                std::cout << std::endl << "Error: can not open " << libInfo.s_q_fname[j] << std::endl;
                                exit(1);
                            }
                            initializeReadsFq(_fin
                                    , libInfo.with_barcode == 1 );
                            _fin.close();
                            _fin.clear();
                        }
                        else{
                            initializeReadsFqGZ(libInfo.s_q_fname[j]
                                    , libInfo.with_barcode == 1 );
                        }
                    }
                }
            }
            else {

                (*finp).clear();
                (*finp).seekg(0, std::ios_base::beg);
                initializeReadsFa(*finp);
            }
        }

        void initializeReadsFa(std::ifstream& fin) {


            //		fin.clear();
            //		fin.seekg(0, ios_base::beg);

            char str[MAX_STRING_LEN];
            while(!fin.eof()) {

                fin.getline(str, MAX_STRING_LEN);
                if (str[0] == '>') {
                    continue;
                }
                else {
                    chomp(str);
                    size_t size = strlen(str);
                    if (size == 0)
                        continue;

                    _initializeReads(str, size);
                }
            }
        }

        void initializeReadsFaGZ(std::string& fin) {

            char *cmd = new char[fin.size()+20];
            sprintf(cmd, "gzip -dc %s", fin.c_str());
            FILE *fp = popen(cmd, "r");

            char str[MAX_STRING_LEN];
            while(fgets(str, MAX_STRING_LEN, fp)) {


                if (str[0] == '>') {
                    continue;
                }
                else {
                    chomp(str);
                    size_t size = strlen(str);
                    if (size == 0)
                        continue;

                    _initializeReads(str, size);
                }
            }

            delete [] cmd;
            pclose(fp);
        }

        void initializeReadsFq(std::ifstream& fin, bool with_bc = false) {

            char name[MAX_STRING_LEN];
            char tmp3[MAX_STRING_LEN];
            char sequence[MAX_STRING_LEN];
            char quality[MAX_STRING_LEN];
            while(!fin.eof()) {

                fin.getline(name, MAX_STRING_LEN);
                fin.getline(sequence, MAX_STRING_LEN);
                fin.getline(tmp3, MAX_STRING_LEN);
                fin.getline(quality, MAX_STRING_LEN);

                chomp(name);
                chomp(sequence);
                chomp(quality);

                size_t size = strlen(sequence);
                if (size == 0)
                    continue;

                if (with_bc)
                {
                    _initializeReads(sequence, size
                            ,GetBarcodeFromReadName(name));
                }
                else
                    _initializeReads(sequence, size);
            }
        }

        void initializeReadsFqGZ(std::string& fin , bool with_bc = false ) {
            char *cmd = new char[fin.size()+20];
            sprintf(cmd, "gzip -dc %s", fin.c_str());
            FILE *fp = popen(cmd, "r");

            char name[MAX_STRING_LEN];
            char sequence[MAX_STRING_LEN];
            char quality[MAX_STRING_LEN];
            while(fgets(name, MAX_STRING_LEN, fp)) {

                fgets(sequence, MAX_STRING_LEN, fp);
                fgets(name, MAX_STRING_LEN, fp);
                fgets(quality, MAX_STRING_LEN, fp);


                chomp(name);
                chomp(sequence);
                chomp(quality);

                size_t size = strlen(sequence);
                if (size == 0)
                    continue;
                if (with_bc)
                {
                    _initializeReads(sequence, size
                            ,GetBarcodeFromReadName(name));
                }
                else

                    _initializeReads(sequence, size);
            }

            delete [] cmd;
            pclose(fp);
        }

        void initializeReadsFaPair(std::ifstream& fin1, std::ifstream& fin2) {

            char str1[MAX_STRING_LEN];
            char str2[MAX_STRING_LEN];
            while(!fin1.eof() && !fin2.eof()) {

                fin1.getline(str1, MAX_STRING_LEN);
                fin1.getline(str1, MAX_STRING_LEN);

                chomp(str1);
                size_t size = strlen(str1);
                if (size == 0)
                    continue;

                _initializeReads(str1, size);

                fin2.getline(str2, MAX_STRING_LEN);
                fin2.getline(str2, MAX_STRING_LEN);

                chomp(str2);
                size = strlen(str2);
                if (size == 0)
                    continue;

                _initializeReads(str2, size);
            }
        }

        void initializeReadsFaPairGZ(std::string& fin1, std::string& fin2) {

            char *cmd1 = new char[fin1.size()+20];
            sprintf(cmd1, "gzip -dc %s", fin1.c_str());
            FILE *fp1 = popen(cmd1, "r");

            char *cmd2 = new char[fin2.size()+20];
            sprintf(cmd2, "gzip -dc %s", fin2.c_str());
            FILE *fp2 = popen(cmd2, "r");

            char str1[MAX_STRING_LEN];
            char str2[MAX_STRING_LEN];
            while(fgets(str1, MAX_STRING_LEN, fp1)&& fgets(str2, MAX_STRING_LEN, fp2)) {

                fgets(str1, MAX_STRING_LEN, fp1);

                chomp(str1);
                size_t size = strlen(str1);
                if (size == 0)
                    continue;

                _initializeReads(str1, size);

                fgets(str2, MAX_STRING_LEN, fp2);

                chomp(str2);
                size = strlen(str2);
                if (size == 0)
                    continue;

                _initializeReads(str2, size);
            }

            delete [] cmd1;
            delete [] cmd2;

            pclose(fp1);
            pclose(fp2);
        }

        void initializeReadsFqPairGZ(std::string& fin1, std::string& fin2 , bool with_bc = false ) {

            char *cmd1 = new char[fin1.size()+20];
            sprintf(cmd1, "gzip -dc %s", fin1.c_str());
            FILE *fp1 = popen(cmd1, "r");

            char *cmd2 = new char[fin2.size()+20];
            sprintf(cmd2, "gzip -dc %s", fin2.c_str());
            FILE *fp2 = popen(cmd2, "r");

            char name1[MAX_STRING_LEN];
            char tmp3[MAX_STRING_LEN];
            char sequence1[MAX_STRING_LEN];
            char quality1[MAX_STRING_LEN];

            char name2[MAX_STRING_LEN];
            char sequence2[MAX_STRING_LEN];
            char quality2[MAX_STRING_LEN];

            while(fgets(name1, MAX_STRING_LEN, fp1) && fgets(name2, MAX_STRING_LEN, fp2)) {

                fgets(sequence1, MAX_STRING_LEN, fp1);
                fgets(tmp3, MAX_STRING_LEN, fp1);
                fgets(quality1, MAX_STRING_LEN, fp1);


                chomp(name1);
                chomp(sequence1);
                chomp(quality1);

                size_t size = strlen(sequence1);
                if (size == 0)
                    continue;

                if (with_bc)
                {
                    _initializeReads(sequence1, size
                            ,GetBarcodeFromReadName(name1));
                }
                else
                    _initializeReads(sequence1, size);

                fgets(sequence2, MAX_STRING_LEN, fp2);
                fgets(tmp3, MAX_STRING_LEN, fp2);
                fgets(quality2, MAX_STRING_LEN, fp2);


                chomp(name2);
                chomp(sequence2);
                chomp(quality2);

                size = strlen(sequence2);
                if (size == 0)
                    continue;
                if (with_bc)
                {
                    _initializeReads(sequence2, size
                            ,GetBarcodeFromReadName(name2));
                }
                else
                    _initializeReads(sequence2, size);
            }

            delete [] cmd1;
            delete [] cmd2;

            pclose(fp1);
            pclose(fp2);
        }

        void initializeReadsFqPair(std::ifstream& fin1, std::ifstream& fin2, bool with_bc = false) {

            char name1[MAX_STRING_LEN];
            char tmp3[MAX_STRING_LEN];
            char sequence1[MAX_STRING_LEN];
            char quality1[MAX_STRING_LEN];

            char name2[MAX_STRING_LEN];
            char sequence2[MAX_STRING_LEN];
            char quality2[MAX_STRING_LEN];

            while(!fin1.eof() && !fin2.eof()) {

                fin1.getline(name1, MAX_STRING_LEN);
                fin1.getline(sequence1, MAX_STRING_LEN);
                fin1.getline(tmp3, MAX_STRING_LEN);
                fin1.getline(quality1, MAX_STRING_LEN);

                chomp(name1);
                chomp(sequence1);
                chomp(quality1);

                size_t size = strlen(sequence1);
                if (size == 0)
                    continue;
                if (with_bc)
                {
                    _initializeReads(sequence1, size
                            ,GetBarcodeFromReadName(name1));
                }
                else
                    _initializeReads(sequence1, size);

                fin2.getline(name2, MAX_STRING_LEN);
                fin2.getline(sequence2, MAX_STRING_LEN);
                fin2.getline(tmp3, MAX_STRING_LEN);
                fin2.getline(quality2, MAX_STRING_LEN);

                chomp(name2);
                chomp(sequence2);
                chomp(quality2);

                size = strlen(sequence2);
                if (size == 0)
                    continue;

                if (with_bc)
                {
                    _initializeReads(sequence2, size,GetBarcodeFromReadName(name2));
                }
                else
                    _initializeReads(sequence2, size);
            }
        }
        void _initializeReads(char* str, size_t size, int barcode_num = 0) {

            while (size < overlapParam) {
                str[size] = 'N';
                size++;
            }
            if (size > Read::DATA_MAXLEN) {
                size = Read::DATA_MAXLEN;
            }
            str[size] = '\0';

            //consider reverse sequence, and store big one
            Number_t* numbers;
            Len_t numLen;
            sequenceToNumbers(str, size, numbers, numLen);

            Number_t* numbersRev;
            reverseComplement(numbers, numbersRev, numLen, size);

            if (compare(numbers, size, numbersRev, size) > 0) {
                reads[readsSum].initialize(numbers, size ,barcode_num);
            }
            else {
                reads[readsSum].initialize(numbersRev, size,barcode_num);
            }

            delete [] numbers;
            delete [] numbersRev;

            readsSum++;
        }

        void sortReadsReverse() {

            //		memcpy(serialNumsReverse, serialNumsForward, readsSum*sizeof(Len_t));

            for (Number_t i=0; i<readsSum; ++i) {
                serialNumsReverse[i] = i;
            }

            QuickSorterMulti<Read, Len_t> sorter;
            sorter.doSort(reads, serialNumsReverse, compare, readsSum);
        }

        void initializeReadsHashReverse() {

            for (Number_t i=0; i<readsSum; ++i) {

                Number_t* readData;
                Len_t arraySize;
                reads[i].getReadData(readData, arraySize);

                Number_t kmer;
                Len_t readLen = reads[i].getLen();
                if (readLen < DATA_LEN) {
                    kmer = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
                }
                else {
                    kmer = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
                }

                if (readsHash.find(kmer)) {

                    if (readsHash[kmer].serialNumReverse.start == (Len_t)~0) {
                        readsHash[kmer].serialNumReverse.start = i;
                    }
                    //				readsHash[kmer].serialNumReverse.count++;
                }
                else {
                    readsHash[kmer].serialNumReverse.start = i;
                    //				readsHash[kmer].serialNumReverse.count = 1;
                }

                delete [] readData;
            }
        }

        void reverseReads() {

            for (Number_t i=0; i<readsSum; ++i) {
                reads[i].reverse();
            }
        }

        void sortReadsForward() {

            //		for (Number_t i=0; i<readsSum; ++i) {
            //			serialNumsForward[i] = i;
            //		}

            memcpy(serialNumsForward, serialNumsReverse, readsSum*sizeof(Len_t));

            QuickSorterMulti<Read, Len_t> sorter;
            sorter.doSort(reads, serialNumsForward, compare, readsSum);
        }

        void initializeReadsHashForward() {

            for (Number_t i=0; i<readsSum; ++i) {

                Number_t* readData;
                Len_t arraySize;
                reads[i].getReadData(readData, arraySize);

                Number_t kmer;
                Len_t readLen = reads[i].getLen();
                if (readLen < DATA_LEN) {
                    kmer = readData[0] >> (readLen-overlapParam)*LEN_TO_BITLEN;
                }
                else {
                    kmer = readData[0] >> (DATA_LEN-overlapParam)*LEN_TO_BITLEN;
                }

                if (readsHash.find(kmer)) {

                    if (readsHash[kmer].serialNumForward.start == (Len_t)~0) {
                        readsHash[kmer].serialNumForward.start = i;
                    }
                    //				readsHash[kmer].serialNumForward.count++;
                }
                else {
                    readsHash[kmer].serialNumForward.start = i;
                    //				readsHash[kmer].serialNumForward.count = 1;
                }

                delete [] readData;
            }
        }

        //asign a same readID to identical reads and record the number of identical reads as depth of this readID
        void uniqueReads() {

            Len_t depth = 0;
            Number_t i;
            for (i=0; i<readsSum-1; ++i) {

                ++depth;
                if (reads[i] != reads[i+1]) {

                    for (Len_t j=i+1-depth; j<i+1; ++j) {

                        reads[j].setDepth(depth);
                        reads[j].setID(i+1-depth);
                    }
                    depth = 0;
                }
            }

            ++depth;
            for (Len_t j=i+1-depth; j<i+1; ++j) {

                reads[j].setDepth(depth);
                reads[j].setID(i+1-depth);
            }
        }

        void sortReadsBySerialNum() {

            Len_t* serialNums = new Len_t[readsSum];
            memcpy(serialNums, serialNumsForward, readsSum*sizeof(Len_t));

            QuickSorterMulti<Len_t, Read> sorter;
            sorter.doSort(serialNums, reads, compare, readsSum);

            delete [] serialNums;
        }

};



#endif /*READHASH_HPP_*/
