# stLFR GapCloser

## <a name=intro>Introduction</a>

This project is a modified version of original GapCloser from SOAPdenovo2[1].
Instread of use NGS reads to fill gap , this GapCloser can use stLFR reads .
We alse change the original extern strategy to get a more accuracy result .

## <a name=table>Table of Contents</a>

- [Introduction](#intro)
- [Table of Contents](#table)
- [User's Guide](#user-guide)
    - [Installation](#install)
    - [Preliminary](#pre)
    - [Quick start](#quick-start)
    - [General usage](#usage)
- [Reference](#ref)

## <a name=user-guide>User's Guide</a>

### <a name=install>Installation</a>

- How to download the source codes.
```
git clone https://github.com/BGI-QingDao/stLFR_Scaffold_Assembler.git YOUR-DOWNLOAD-DIR
```
- How to compiler source codes .
```
> cd YOUR-DOWNLOAD-DIR
> make
> YOUR-DOWNLOAD-DIR/Release/GapCloser # will run the program and show usage.
```

### <a name=pre>Preliminary</a> 

- the stLFR reads must have and only have 2 file :
    - your-prefix.read1.your-suffix
    - your-prefix.read2.your-suffix

  *We assume your stLFR reads is the reads that after barcode splitted.*
*The official barcode split step is the "1.fq_BarcodeSplit" step of stLFR_v1(https://github.com/MGI-tech-bioinformatics/stLFR_v1.git)*

  *If you don't want to download the stLFR_v1 package , you can try our split barcode script:*

```
//TODO
```

### <a name=quick-start>Quick start</a>

- 1st. prepare the lib.cfg
*The format comes from SOAPdenovo, see details from https://github.com/aquaskyline/SOAPdenovo2 . Here is a basic example :*

```
#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=300
#if sequence needs to be reversed
reverse_seq=0
#a pair of fastq file, read 1 file should always be followed by read 2 file
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
        -D      <int>           max number of conflict can accepted , default=2.
        -E      <int>           min depth threshold, default=1.
        -F      <int>           max number of low depth nucleotide, default=1.
other:
        -G      <int>           the max number of reads that a kmer can find . default=50.
        -I      <int>           use sub-reads-set only. default=0.
        -O      <int>           [ basic set ] max number of conflict can accepted . default=2.
        -P      <int>           [ basic set ] max number of low depth nucleotide . default=1.
        -Q      <float>         [ basic set ] the conflict_threshold of baisc set . default=0.8.
        -R      <int>           max read-find-read round . default=1.
 ---------- new parameters end     ---------
        -h      -?              output help information.
```


## <a name=ref>Reference</a>

[1] [SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler][11]

[11]:https://www.ncbi.nlm.nih.gov/pubmed/23587118

