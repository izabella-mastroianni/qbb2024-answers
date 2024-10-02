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

# 

## Question 1.4 ##


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

## Question 2.2 ##

# for loop to align reads to reference genome 
for my_sample in *.fastq
    do 
        my_sample_basename=`basename ${my_sample} .fastq`
        read_group="@RG\tID:${my_sample_basename}\tSM:${my_sample_basename}\tPL:illumina"
        bwa mem -R "${read_group}" sacCer3.fa ${my_sample_basename}.fastq > ${my_sample_basename}.sam
        samtools sort -@ 4 -O bam -o "${my_sample_basename}.bam" "${my_sample_basename}.sam"
        samtools index ${my_sample_basename}.bam
    done





