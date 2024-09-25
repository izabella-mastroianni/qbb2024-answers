#!/usr/bin/bash


# week 2 exercise 2

# 2.1
# 4 feature files are merged_exons, introns, merged_cCREs, and other
# note which part of pseudocode each part of following code refers to

# create output results file called snp_counts.txt
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

# loop through each possible MAF file
snp_files=("chr1_snps_0.*.bed") # maybe change wild card
feature_files=("chr1_cCREs.bed", "chr1_exons.bed", "genome_chr1.bed", "other_chr1.bed") 
# per snp file (aka MAF value) per feature file

for file in snp_files
    do
        # SNP coverage of the whole chromosome 
        bedtools coverage -a file -b genome_chr1.bed > temp_coverage.txt
        # sum SNPs from coverage
        coverage_SNPs_sum=$(awk '{s+=$1}END{print s}' temp_coverage.txt)
        # sum total bases from coverage


    done 









# notes: scale = SNPs per base, for each  base we have a probability for that SNP covering that base
# want coverage over whole genome (use whole genome file)
# for each MAF file, calculating same value across genome file
# inner loop subdivide by exons introns and regulator elements
# coverage and calculate sums

