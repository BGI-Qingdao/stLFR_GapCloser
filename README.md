# stLFR GapCloser

## <a name=intro>Introduction</a>

This project is a advanced version of the original GapCloser from SOAPdenovo2[1].
Instead of using PE information from NGS reads, this version use barcode information from stLFR reads to fill gaps.
We alse change the original extendsion strategy to get more accurate results.

## <a name=table>Table of Contents</a>

- [Introduction](#intro)
- [Table of Contents](#table)
- [User's Guide](#user-guide)
    - [Installation](#install)
    - [Preliminary](#pre)
    - [Quick start](#quick-start)
    - [General usage](#usage)
- [Reference](#ref)
- [Contact](#contact)

## <a name=user-guide>User's Guide</a>

### <a name=install>Installation</a>

- How to download the source codes.
```
git clone https://github.com/BGI-QingDao/stLFR_Scaffold_Assembler.git YOUR-DOWNLOAD-DIR
```
- How to compile source codes .
```
> cd YOUR-DOWNLOAD-DIR
> make
> YOUR-DOWNLOAD-DIR/Release/GapCloser # will run the program and show usage.
```

### <a name=pre>Preliminary</a> 

- the stLFR reads are required as 2 files :
    - your-prefix.read1.your-suffix
    - your-prefix.read2.your-suffix

  *We assume your stLFR reads and barcode information have already been splitted.*
  
  *If your data have not been splitted yet, then use the split barcode script below:*
```
# if your raw stLFR reads contain more than 1 lane, you need to cat all lines into a single file first!
YOUR-INSTALL-DIR/split_barcode/split_barcode.sh raw_read1.fq.gz raw_read2.fq.gz
```
*Also, you can try "1.fq_BarcodeSplit" step from stLFR_v1(https://github.com/MGI-tech-bioinformatics/stLFR_v1.git)*


### <a name=quick-start>Quick start</a>

- 1st. prepare the lib.cfg
*The format comes from SOAPdenovo, see details from https://github.com/aquaskyline/SOAPdenovo2 . Here is a basic example:*

```
#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=300
#if sequence needs to be reversed
reverse_seq=0
#a pair of fastq files, read1 file should be followed by read2 file
q1=split_reads.1.fq.gz
q2=split_reads.2.fq.gz
```
- 2nd. run the program

```
YOUR-DOWNLOAD-DIR/Release/GapCloser -a your-scaffold-file -b lib.cfg -o result.scaff
```

### <a name=usage>General usage</a>

```
Usage:
        GapCloser [options]
Basic options :
        -a      <string>        input scaffold file name, required.
        -b      <string>        input library info file name, required.
        -o      <string>        output file name, required.
        -l      <int>           maximum read length (<=155), default=150.
        -p      <int>           overlap param(<=31) [the kvalue], default=27
Performance options :
        -t      <int>           thread number, default=1.
        -c      <float>         hash load fractor, default=0.75.
mapping read to contig options:
        -1      <int>           maximum read depth, default=100.
        -2      <int>           maximum mismatch , default=1.
consensus reads set options:
        -3      <int>           consensus length, default= 41.
        -4      <int>           consensus prev extra length, default=57.
        -5      <int>           consensus last extra length, default=1.
        -6      <int>           consensus extend length, default=40.
consensus reads set options:
        -7      <int>           min reads count , default=1.
        -8      <int>           max reads count , default=5000.
extract sub reads set options:
        -A      <int>           min pe subset reads count threshold , default=1.
        -B      <int>           min pe & barcode subset reads count threshold , default=1.
consensus options:
        -C      <float>         non-conflict threshold , default=0.6.
        -D      <int>           max number of acceptable conflicts, default=2.
        -E      <int>           min depth threshold, default=1.
        -F      <int>           max number of low depth nucleotide, default=1.
other:
        -G      <int>           max number of reads that a kmer can find. default=50.
        -I      <int>           use sub-reads-set only. default=0.
        -O      <int>           [ basic set ] max number of acceptable conflicts. default=2.
        -P      <int>           [ basic set ] max number of low depth nucleotide. default=1.
        -Q      <float>         [ basic set ] the conflict_threshold of baisc set. default=0.8.
        -R      <int>           max read-find-read round . default=1.
 ---------- new parameters end     ---------
        -h      -?              output help information.
```


## <a name=ref>Reference</a>

[1] [SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler][11]

[11]:https://www.ncbi.nlm.nih.gov/pubmed/23587118

## <a name=contact>Contact</a>

- for algrothim detail & discussion
    - please contact dengli1@genomics.cn
- for code detail & bug report 
    - please contact guolidong@genomics.cn 
    - or please creat an issue in github.
