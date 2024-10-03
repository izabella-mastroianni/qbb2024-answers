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
FASTQC A01_09.fastq

# The median base quality along the read is around 34. 
# This means that the probabilty a given base called is an error is around 0.1%.
# I do observe some minimal variation in quality with respect to the position in the read, with slightly worse at the ends and best in the middle. 

