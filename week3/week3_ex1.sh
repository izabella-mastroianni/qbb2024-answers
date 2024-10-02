#!/usr/bin/env bash

# conda activate qb24 to run



## Question 1.1 ##

# use FASTQC - produces html, only argument is name of fastq file
# open html file from finder 
FASTQC A01_09.fastq

# The sequencing reads are 76 base pairs long


## Question 1.2 ##

# each read in a fastq file takes up 4 lines
wc -l A01_09.fastq # 2678192
expr 2678192 / 4  # 669548

# There are 669548 reads in the A01_09.fastq file


## Question 1.3 ##

# find lines that start with > which are sequences and then count to get size of reference genome 
grep -v '^>' sacCer3.fa | tr -d '\n' | wc -c

# the S. cerevisiae reference genome is 12157105

# expected average depth of coverage genome length divided by reads
expr 12157105 / 2678192 #4
# the depth of coverage is 4

## Question 1.4 ##

# check file size of fastq files, h is human readable
du -h *.fastq

# The largest file is A01_62.fastq, 149Mb. The smallest file is A01_27.fastq, 110Mb. 

## Question 1.5 ##



## Question 2.1 ##

# commented out lines below so didn't download and run every time

# download sacCer3 reference genome 
# wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
# gunzip sacCer3.fa.gz

# create index for sacCer3 reference 
# bwa index sacCer3.fa

# find number of chromosome using count of number of times chr appears
# grep -c 'chr' sacCer3.fa
# gives 17, but google says 16, if grep for 'chrXVII' get no result, but get one for 'chrXVI'
# grep 'chrXVII' sacCer3.fa

# There are 17 chromosomes in the yeast genome but according to google it should be 16 so i'm not sure


# for loop to align reads to reference genome, also includes step 2.4
for my_sample in *.fastq
    do 
        my_sample_basename=`basename ${my_sample} .fastq`
        read_group="@RG\tID:${my_sample_basename}\tSM:${my_sample_basename}\tPL:illumina"
        bwa mem -R "${read_group}" sacCer3.fa ${my_sample_basename}.fastq > ${my_sample_basename}.sam
        samtools sort -@ 4 -O bam -o "${my_sample_basename}.bam" "${my_sample_basename}.sam"
        samtools index ${my_sample_basename}.bam
    done

## Question 2.2 ##

# to find the number of reads, count number of lines with an @ which correspond to reads
grep -v '^@' A01_09.sam | wc -l

# In the A01_09 sam file, there are 669548 total read alignments 

## Question 2.3 ##

# count how many lines with an @ and III to specify for chromosome III
grep -v '^@' A01_09.sam | grep -c '^III'

# In the A01_09 sam file, there are 0 alignments to loci on chromosome III

## Step 2.4 ##





