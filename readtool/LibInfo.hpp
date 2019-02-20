/**************************************************
*
* LibInfo.hpp
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

#ifndef LIBINFO_HPP_
#define LIBINFO_HPP_


/* Basic information of a library and some members contributed to go through all read files of a library */
struct LibInfoElement
{
	// library name
	string name;

	// read number infomation
	Number_t readsSum; // total read number
	Number_t pairReadsSum;
	Number_t singleReadsSum;

	// insert size information
	int min_ins;
	int max_ins;
	int avg_ins;
	
	int rd_len_cutoff;
	int reverse;
	int asm_flag;
	
	// indicate which file is next to read
	int curr_type; // 1 - 5
	int curr_index; // file index of current type
	
	// file handlers to opened files
	ifstream* fp1;
	ifstream* fp2;
	// whether last read is read1 in pair
	int paired;  // 0 -- single; 1 -- read1; 2 -- read2;

	int libUsed;

	// type1, paired-end reads in fasta format with read1 and read2 sperated in two files
	//The order of read1 of a pair of reads in read1 file should equals to read2's order in read2 file.
	string* a1_fname; // read1 file name
	string* a2_fname;
	int num_a1_file; // read1 file number
	int num_a2_file;
	Number_t *a1_seqNum; // read1 number
	Number_t *a2_seqNum;
	int *a1_used;
	int *a2_used;
	

	// type2, paired-end reads in fastq format
	string* q1_fname;
	string* q2_fname;
	int num_q1_file;
	int num_q2_file;
	Number_t *q1_seqNum;
	Number_t *q2_seqNum;
	int *q1_used;
	int *q2_used;

	// type3, paired-end reads in fasta format with read1 and read2 in one file.
	// Read1 should be always followed by read2.
	string* p_fname;
	int num_p_file; 
	Number_t *p_seqNum;
	int *p_used;

	// type4 &5, single read in fasta or fastq format
	string* s_a_fname;
	int num_s_a_file;
	string* s_q_fname;
	int num_s_q_file;
	Number_t *s_a_seqNum;
	Number_t *s_q_seqNum;
	int *s_a_used;
	int *s_q_used;

	Number_t startNum;
	
};


/* Information of all libraries and some operations on library including libraries comparison and libraries initialization. */
class LibInfo
{
	
protected:

	// config file
	ifstream* finp;

	// total read number
	Number_t readsSum;
	
	Len_t maxReadLen;
	
	LibInfoElement* libInfoArray;
	Len_t libsNum;
	char tabs[2][1024];

	// start positon of every library in read array
	Number_t *startNum;
	
public:
	
	LibInfo() : 
		finp(NULL), 
		libInfoArray(NULL), 
		libsNum(0)
	{
	}
	
	LibInfo(ifstream* _finp) : 
		finp(_finp), 
		libInfoArray(NULL), 
		libsNum(0)
	{
		scanLibInfo(*finp);
		
	}
	
	virtual ~LibInfo() {
		
		freeLibs();
	}
	
	Number_t getReadsSum() const { return readsSum; }
	Len_t getMaxReadLen() const { return maxReadLen; }
	
	LibInfoElement* getArray() const { return libInfoArray; }
	Len_t getArrayLen() const { return libsNum; }
	Number_t* getStartNum() { return startNum; }

	//split a string with sperator '='
	bool splitColumn(char const* line) {
		
		int len = strlen(line);
		int i=0,j;
		int tabs_n = 0;
		
		while(i<len){
			if(line[i]>=32&&line[i]<=126&&line[i]!='='){
				j=0;
				while(i<len&&line[i]>=32&&line[i]<=126&&line[i]!='='){
					tabs[tabs_n][j++] = line[i];	
					i++;
				}
				tabs[tabs_n][j] = '\0';	
				tabs_n++;
				if(tabs_n==2)
					return true;
			}
			i++;
		}
		if(tabs_n==2)
			return true;
		else
			return false;
	}
	
	static int cmpLib(const void* a, const void* b)
	{
		LibInfoElement* A,* B;
		A = (LibInfoElement*)a;
		B = (LibInfoElement*)b;

		if(A->avg_ins>B->avg_ins)
			return 1;
		else if(A->avg_ins==B->avg_ins)
			return 0;
		else
			return -1;
	}
	
	/* get libraries information in config file specified by 'fin' */
	void scanLibInfo(ifstream& fin) {
		
		const char* delim;
		size_t found;
		string str, str1;
		int i,index;
		bool flag;
		
		//count lib
		while(!fin.eof()) {
			getline(fin, str);
			if (str.find("[LIB]") != string::npos) {
				libsNum++;
			}
			if(!libsNum){
				
				delim = "max_rd_len=";
				found = str.find(delim);
				if (found!=string::npos) {
					
					str1 = str.substr(found+strlen(delim));
					maxReadLen = atoi(str1.c_str());
					continue;
				}
				delim = "reads_sum=";
				found = str.find(delim);
				if (found!=string::npos) {
					
					str1 = str.substr(found+strlen(delim));
					readsSum = atoll(str1.c_str());
					continue;
				}
			}
		}
		
		//count file numbers of each type
		libInfoArray = new LibInfoElement[libsNum];
		startNum = new Number_t[libsNum+1];
		for(i=0;i<(int)libsNum;i++){
			libInfoArray[i].asm_flag = 3;
			libInfoArray[i].num_s_a_file=0;	
			libInfoArray[i].num_s_q_file=0;	
			libInfoArray[i].num_p_file=0;	
			
			libInfoArray[i].num_a1_file=0;	
			libInfoArray[i].num_a2_file=0;	
			libInfoArray[i].num_q1_file=0;	
			libInfoArray[i].num_q2_file=0;	

			libInfoArray[i].libUsed = 0;

			startNum[i] = 0;
		}

		startNum[i] = 0;	//for single reads
		
		fin.clear();
		fin.seekg(0, ios_base::beg);
		
		i = -1;
		while(!fin.eof()) {
			getline(fin, str);
			if (str.find("[LIB]") != string::npos) {
				i++;
				continue;
			}
			
			if(i<0)
				continue;
			
			flag = splitColumn(str.c_str());
			if(!flag)
				continue;

			if(strcmp(tabs[0],"f1")==0){
				libInfoArray[i].num_a1_file++;
			}else if(strcmp(tabs[0],"q1")==0){
				libInfoArray[i].num_q1_file++;
			}else if(strcmp(tabs[0],"f2")==0){
				libInfoArray[i].num_a2_file++;
			}else if(strcmp(tabs[0],"q2")==0){
				libInfoArray[i].num_q2_file++;
			}else if(strcmp(tabs[0],"f")==0){
				libInfoArray[i].num_s_a_file++;
			}else if(strcmp(tabs[0],"q")==0){
				libInfoArray[i].num_s_q_file++;
			}else if(strcmp(tabs[0],"p")==0){
				libInfoArray[i].num_p_file++;
			}
		}
		
		//allocate memory for filenames
		for(i=0;i<(int)libsNum;i++){
			
			if(libInfoArray[i].num_s_a_file){
				libInfoArray[i].s_a_fname = new string[libInfoArray[i].num_s_a_file];
				libInfoArray[i].s_a_seqNum = new Number_t[libInfoArray[i].num_s_a_file];
				libInfoArray[i].s_a_used = new int[libInfoArray[i].num_s_a_file];
			}
			if(libInfoArray[i].num_s_q_file){
				libInfoArray[i].s_q_fname = new string[libInfoArray[i].num_s_q_file];
				libInfoArray[i].s_q_seqNum = new Number_t[libInfoArray[i].num_s_q_file];
				libInfoArray[i].s_q_used = new int[libInfoArray[i].num_s_q_file];
			}
			if(libInfoArray[i].num_p_file){
				libInfoArray[i].p_fname = new string[libInfoArray[i].num_p_file];
				libInfoArray[i].p_seqNum = new Number_t[libInfoArray[i].num_p_file];
				libInfoArray[i].p_used = new int[libInfoArray[i].num_p_file];
			}
			if(libInfoArray[i].num_a1_file){
				libInfoArray[i].a1_fname = new string[libInfoArray[i].num_a1_file];
				libInfoArray[i].a1_seqNum = new Number_t[libInfoArray[i].num_a1_file];
				libInfoArray[i].a1_used = new int[libInfoArray[i].num_a1_file];
			}
			if(libInfoArray[i].num_a2_file){
				libInfoArray[i].a2_fname = new string[libInfoArray[i].num_a2_file];
				libInfoArray[i].a2_seqNum = new Number_t[libInfoArray[i].num_a2_file];
				libInfoArray[i].a2_used = new int[libInfoArray[i].num_a2_file];
			}
			if(libInfoArray[i].num_q1_file){
				libInfoArray[i].q1_fname = new string[libInfoArray[i].num_q1_file];
				libInfoArray[i].q1_seqNum = new Number_t[libInfoArray[i].num_q1_file];
				libInfoArray[i].q1_used = new int[libInfoArray[i].num_q1_file];
			}
			if(libInfoArray[i].num_q2_file){
				libInfoArray[i].q2_fname = new string[libInfoArray[i].num_q2_file];
				libInfoArray[i].q2_seqNum = new Number_t[libInfoArray[i].num_q2_file];
				libInfoArray[i].q2_used = new int[libInfoArray[i].num_q2_file];
			}
		}
		
		// get file names
		for(i=0;i<(int)libsNum;i++){
			libInfoArray[i].readsSum=0;
			libInfoArray[i].pairReadsSum=0;
			libInfoArray[i].curr_type=1;
			libInfoArray[i].curr_index=0;
			libInfoArray[i].fp1=NULL;
			libInfoArray[i].fp2=NULL;

			libInfoArray[i].num_s_a_file=0;	
			libInfoArray[i].num_s_q_file=0;	
			libInfoArray[i].num_p_file=0;	
			
			libInfoArray[i].num_a1_file=0;	
			libInfoArray[i].num_a2_file=0;	
			libInfoArray[i].num_q1_file=0;	
			libInfoArray[i].num_q2_file=0;	
		}

		fin.clear();
		fin.seekg(0, ios_base::beg);
		
		i = -1;
		while(!fin.eof()) {
			getline(fin, str);
			if (str.find("[LIB]") != string::npos) {
				i++;
				continue;
			}
			
			if(i<0)
				continue;
			
			flag = splitColumn(str.c_str());
			if(!flag)
				continue;

			if(strcmp(tabs[0],"f1")==0){
				index = libInfoArray[i].num_a1_file++;
				libInfoArray[i].a1_fname[index] = tabs[1];
				libInfoArray[i].a1_seqNum[index] = 0;
				libInfoArray[i].a1_used[index] = 0;
			}else if(strcmp(tabs[0],"q1")==0){
				index = libInfoArray[i].num_q1_file++;
				libInfoArray[i].q1_fname[index] = tabs[1];
				libInfoArray[i].q1_seqNum[index] = 0;
				libInfoArray[i].q1_used[index] = 0;
			}else if(strcmp(tabs[0],"f2")==0){
				index = libInfoArray[i].num_a2_file++;
				libInfoArray[i].a2_fname[index] = tabs[1];
				libInfoArray[i].a2_seqNum[index] = 0;
				libInfoArray[i].a2_used[index] = 0;
			}else if(strcmp(tabs[0],"q2")==0){
				index = libInfoArray[i].num_q2_file++;
				libInfoArray[i].q2_fname[index] = tabs[1];
				libInfoArray[i].q2_seqNum[index] = 0;
				libInfoArray[i].q2_used[index] = 0;
			}else if(strcmp(tabs[0],"f")==0){
				index = libInfoArray[i].num_s_a_file++;
				libInfoArray[i].s_a_fname[index] = tabs[1];
				libInfoArray[i].s_a_seqNum[index] = 0;
				libInfoArray[i].s_a_used[index] = 0;
			}else if(strcmp(tabs[0],"q")==0){
				index = libInfoArray[i].num_s_q_file++;
				libInfoArray[i].s_q_fname[index] = tabs[1];
				libInfoArray[i].s_q_seqNum[index] = 0;
				libInfoArray[i].s_q_used[index] = 0;
			}else if(strcmp(tabs[0],"p")==0){
				index = libInfoArray[i].num_p_file++;
				libInfoArray[i].p_fname[index] = tabs[1];
				libInfoArray[i].p_seqNum[index] = 0;
				libInfoArray[i].p_used[index] = 0;
			}else if(strcmp(tabs[0],"name")==0)
				libInfoArray[i].name = tabs[1];
			else if(strcmp(tabs[0],"reads_sum")==0)
				libInfoArray[i].readsSum = atoi(tabs[1]);
			else if(strcmp(tabs[0],"pair_reads_sum")==0)
				libInfoArray[i].pairReadsSum = atoi(tabs[1]);
			else if(strcmp(tabs[0],"min_ins")==0)
				libInfoArray[i].min_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"max_ins")==0)
				libInfoArray[i].max_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"avg_ins")==0)
				libInfoArray[i].avg_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"rd_len_cutoff")==0)
				libInfoArray[i].rd_len_cutoff = atoi(tabs[1]);
			else if(strcmp(tabs[0],"reverse_seq")==0)
				libInfoArray[i].reverse = atoi(tabs[1]);
			else if(strcmp(tabs[0],"asm_flags")==0)
				libInfoArray[i].asm_flag = atoi(tabs[1]);
		}
		
		qsort(&libInfoArray[0],libsNum,sizeof(LibInfoElement),cmpLib);
	}
	
	void freeLibs() {

		if(!libInfoArray)
			return;

		int i;
		for(i=0;i<(int)libsNum;i++){
			if(libInfoArray[i].num_s_a_file){
				delete [] libInfoArray[i].s_a_fname;
				delete [] libInfoArray[i].s_a_seqNum;
				delete [] libInfoArray[i].s_a_used;
			}
			if(libInfoArray[i].num_s_q_file){
				delete [] libInfoArray[i].s_q_fname;
				delete [] libInfoArray[i].s_q_seqNum;
				delete [] libInfoArray[i].s_q_used;
			}
			if(libInfoArray[i].num_p_file){
				delete [] libInfoArray[i].p_fname;
				delete [] libInfoArray[i].p_seqNum;
				delete [] libInfoArray[i].p_used;
			}
			if(libInfoArray[i].num_a1_file){
				delete [] libInfoArray[i].a1_fname;
				delete [] libInfoArray[i].a1_seqNum;
				delete [] libInfoArray[i].a1_used;
			}
			if(libInfoArray[i].num_a2_file){
				delete [] libInfoArray[i].a2_fname;
				delete [] libInfoArray[i].a2_seqNum;
				delete [] libInfoArray[i].a2_used;
			}
			if(libInfoArray[i].num_q1_file){
				delete [] libInfoArray[i].q1_fname;
				delete [] libInfoArray[i].q1_seqNum;
				delete [] libInfoArray[i].q1_used;
			}
			if(libInfoArray[i].num_q2_file){
				delete [] libInfoArray[i].q2_fname;
				delete [] libInfoArray[i].q2_seqNum;
				delete [] libInfoArray[i].q2_used;
			}
		}
		
		libsNum = 0;
		delete [] libInfoArray;

		if (startNum != NULL)
			delete [] startNum;
	}
	
	
	
	/*
	void scanLibInfo(char *libfile)
	{
		FILE *fp;
		char line[1024],ch;
		int i,j,index;
		int libCounter;
		bool flag;
		
		fp = ckopen(libfile,"r");
		libsNum = 0;
		while(fgets(line,1024,fp)){
			ch = line[5];
			line[5] = '\0';
			if(strcmp(line,"[LIB]")==0)
				libsNum++;
			if(!libsNum){
				line[5] = ch;
				flag = splitColumn(line);
				if(!flag)
					continue;
				if(strcmp(tabs[0],"max_rd_len")==0)
					maxReadLen = atoi(tabs[1]);
			}
		}
	//count file numbers of each type
		libInfoArray = (LibInfoElement *)ckalloc(libsNum*sizeof(LibInfoElement));
		for(i=0;i<libsNum;i++){
			libInfoArray[i].asm_flag = 3;
			libInfoArray[i].num_s_a_file=0;	
			libInfoArray[i].num_s_q_file=0;	
			libInfoArray[i].num_p_file=0;	
			
			libInfoArray[i].num_a1_file=0;	
			libInfoArray[i].num_a2_file=0;	
			libInfoArray[i].num_q1_file=0;	
			libInfoArray[i].num_q2_file=0;	
		}
		libCounter = -1;
		rewind(fp);
		
		i = -1;
		while(fgets(line,1024,fp)){
			ch = line[5];
			line[5] = '\0';
			if(strcmp(line,"[LIB]")==0){
				i++;
				continue;
			}
			line[5] = ch;
			flag = splitColumn(line);
			if(!flag)
				continue;

			if(strcmp(tabs[0],"f1")==0){
				libInfoArray[i].num_a1_file++;
			}else if(strcmp(tabs[0],"q1")==0){
				libInfoArray[i].num_q1_file++;
			}else if(strcmp(tabs[0],"f2")==0){
				libInfoArray[i].num_a2_file++;
			}else if(strcmp(tabs[0],"q2")==0){
				libInfoArray[i].num_q2_file++;
			}else if(strcmp(tabs[0],"f")==0){
				libInfoArray[i].num_s_a_file++;
			}else if(strcmp(tabs[0],"q")==0){
				libInfoArray[i].num_s_q_file++;
			}else if(strcmp(tabs[0],"p")==0){
				libInfoArray[i].num_p_file++;
			}
		}
	//allocate memory for filenames
		for(i=0;i<libsNum;i++){
			
			if(libInfoArray[i].num_s_a_file){
				libInfoArray[i].s_a_fname = (char **)ckalloc(libInfoArray[i].num_s_a_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_s_a_file;j++)
					libInfoArray[i].s_a_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_s_q_file){
				libInfoArray[i].s_q_fname = (char **)ckalloc(libInfoArray[i].num_s_q_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_s_q_file;j++)
					libInfoArray[i].s_q_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_p_file){
				libInfoArray[i].p_fname = (char **)ckalloc(libInfoArray[i].num_p_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_p_file;j++)
					libInfoArray[i].p_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_a1_file){
				libInfoArray[i].a1_fname = (char **)ckalloc(libInfoArray[i].num_a1_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_a1_file;j++)
					libInfoArray[i].a1_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_a2_file){
				libInfoArray[i].a2_fname = (char **)ckalloc(libInfoArray[i].num_a2_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_a2_file;j++)
					libInfoArray[i].a2_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_q1_file){
				libInfoArray[i].q1_fname = (char **)ckalloc(libInfoArray[i].num_q1_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_q1_file;j++)
					libInfoArray[i].q1_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
			if(libInfoArray[i].num_q2_file){
				libInfoArray[i].q2_fname = (char **)ckalloc(libInfoArray[i].num_q2_file*sizeof(char *));
				for(j=0;j<libInfoArray[i].num_q2_file;j++)
					libInfoArray[i].q2_fname[j] = (char *)ckalloc(1024*sizeof(char));
			}
		}

	// get file names
		for(i=0;i<libsNum;i++){
			libInfoArray[i].curr_type=1;
			libInfoArray[i].curr_index=0;
			libInfoArray[i].fp1=NULL;
			libInfoArray[i].fp2=NULL;

			libInfoArray[i].num_s_a_file=0;	
			libInfoArray[i].num_s_q_file=0;	
			libInfoArray[i].num_p_file=0;	
			
			libInfoArray[i].num_a1_file=0;	
			libInfoArray[i].num_a2_file=0;	
			libInfoArray[i].num_q1_file=0;	
			libInfoArray[i].num_q2_file=0;	
		}
		libCounter = -1;
		rewind(fp);
		
		i = -1;
		while(fgets(line,1024,fp)){
			ch = line[5];
			line[5] = '\0';
			if(strcmp(line,"[LIB]")==0){
				i++;
				continue;
			}
			line[5] = ch;
			flag = splitColumn(line);
			if(!flag)
				continue;

			if(strcmp(tabs[0],"f1")==0){
				index = libInfoArray[i].num_a1_file++;
				strcpy(libInfoArray[i].a1_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"q1")==0){
				index = libInfoArray[i].num_q1_file++;
				strcpy(libInfoArray[i].q1_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"f2")==0){
				index = libInfoArray[i].num_a2_file++;
				strcpy(libInfoArray[i].a2_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"q2")==0){
				index = libInfoArray[i].num_q2_file++;
				strcpy(libInfoArray[i].q2_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"f")==0){
				index = libInfoArray[i].num_s_a_file++;
				strcpy(libInfoArray[i].s_a_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"q")==0){
				index = libInfoArray[i].num_s_q_file++;
				strcpy(libInfoArray[i].s_q_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"p")==0){
				index = libInfoArray[i].num_p_file++;
				strcpy(libInfoArray[i].p_fname[index],tabs[1]);
			}else if(strcmp(tabs[0],"min_ins")==0)
				libInfoArray[i].min_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"max_ins")==0)
				libInfoArray[i].max_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"avg_ins")==0)
				libInfoArray[i].avg_ins = atoi(tabs[1]);
			else if(strcmp(tabs[0],"rd_len_cutoff")==0)
				libInfoArray[i].rd_len_cutoff = atoi(tabs[1]);
			else if(strcmp(tabs[0],"reverse_seq")==0)
				libInfoArray[i].reverse = atoi(tabs[1]);
			else if(strcmp(tabs[0],"asm_flags")==0)
				libInfoArray[i].asm_flag = atoi(tabs[1]);
			
		}

		fclose(fp);

		qsort(&libInfoArray[0],libsNum,sizeof(LibInfoElement),cmpLib);
	}


	void freeLibs()
	{

		if(!libInfoArray)
			return;

		int i,j;
		for(i=0;i<libsNum;i++){
			printf("[LIB] %d, avg_ins %d, reverse %d \n",i,libInfoArray[i].avg_ins,libInfoArray[i].reverse);
			if(libInfoArray[i].num_s_a_file){
				//printf("%d single fasta files\n",libInfoArray[i].num_s_a_file);
				for(j=0;j<libInfoArray[i].num_s_a_file;j++)
					free((void *)libInfoArray[i].s_a_fname[j]);
				free((void *)libInfoArray[i].s_a_fname);
			}
			if(libInfoArray[i].num_s_q_file){
				//printf("%d single fastq files\n",libInfoArray[i].num_s_q_file);
				for(j=0;j<libInfoArray[i].num_s_q_file;j++)
					free((void *)libInfoArray[i].s_q_fname[j]);
				free((void *)libInfoArray[i].s_q_fname);
			}
			if(libInfoArray[i].num_p_file){
				//printf("%d paired fasta files\n",libInfoArray[i].num_p_file);
				for(j=0;j<libInfoArray[i].num_p_file;j++)
					free((void *)libInfoArray[i].p_fname[j]);
				free((void *)libInfoArray[i].p_fname);
			}
			if(libInfoArray[i].num_a1_file){
				//printf("%d read1 fasta files\n",libInfoArray[i].num_a1_file);
				for(j=0;j<libInfoArray[i].num_a1_file;j++)
					free((void *)libInfoArray[i].a1_fname[j]);
				free((void *)libInfoArray[i].a1_fname);
			}
			if(libInfoArray[i].num_a2_file){
				//printf("%d read2 fasta files\n",libInfoArray[i].num_a2_file);
				for(j=0;j<libInfoArray[i].num_a2_file;j++)
					free((void *)libInfoArray[i].a2_fname[j]);
				free((void *)libInfoArray[i].a2_fname);
			}
			if(libInfoArray[i].num_q1_file){
				//printf("%d read1 fastq files\n",libInfoArray[i].num_q1_file);
				for(j=0;j<libInfoArray[i].num_q1_file;j++)
					free((void *)libInfoArray[i].q1_fname[j]);
				free((void *)libInfoArray[i].q1_fname);
			}
			if(libInfoArray[i].num_q2_file){
				//printf("%d read2 fastq files\n",libInfoArray[i].num_q2_file);
				for(j=0;j<libInfoArray[i].num_q2_file;j++)
					free((void *)libInfoArray[i].q2_fname[j]);
				free((void *)libInfoArray[i].q2_fname);
			}
		}
		
		libsNum = 0;
		free((void *)libInfoArray);
	}*/
};




#endif /*LIBINFO_HPP_*/
