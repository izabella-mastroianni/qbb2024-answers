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



