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

## <a name=ref>Reference</a>

[1] [SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler][11]

[11]:https://www.ncbi.nlm.nih.gov/pubmed/23587118

