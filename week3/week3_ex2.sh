#!/usr/bin/env bash

# conda activate qb24 to run



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


# for loop to align reads to reference genome, steps 2.2-2.4
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

## Question 2.4 ##
# loaded in A01_09.bam file, zoomed in all the way
# The depth of coverage appears to match my estimate of 4x from 1.3 because some areas have just a few reads and some have none but the average seems to be around 4

## Question 2.5 ##
# I see 3 SNPs that I am reasonably confident in because they are present in all reads that align to the position
# There are potentially a few other SNPs that show up in one or two of the aligned reads, but more sequencing depth is needed to be certain if it's a SNP or sequencing error

## Question 2.6 ##
# The SNP position is chrIV:825,834
# It falls outside a gene because there is no overlapping Refseq Gene below it, and the closest has a stop codon upstream of the SNP. 







