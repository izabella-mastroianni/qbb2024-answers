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

# download sacCer3 reference genome 
wget https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz

# create index for sacCer3 reference 
bwa index sacCer3.fa

# There are chromosomes in the yeast genome 








