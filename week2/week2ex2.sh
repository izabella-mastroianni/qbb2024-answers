#!/usr/bin/bash


# week 2 exercise 2

# 2.1
# 4 feature files are merged_exons, introns, merged_cCREs, and other
# note which part of pseudocode each part of following code refers to

# create output results file called snp_counts.txt
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

# loop through each possible MAF file
snp_files=("chr1_snps_0.*.bed") # maybe change wild card if it doesn't work
feature_files=("chr1_cCREs.bed" "chr1_exons.bed" "chr1_genes.bed" "other_chr1.bed") 
genome_file=("genome_chr1.bed")
# per snp file (aka MAF value) per feature file

for snpfile in "${snp_files[@]}"
    do
        # SNP coverage of the whole chromosome 
        bedtools coverage -a $genome_file -b $snpfile > temp_coverage.txt
        # sum SNPs from coverage - don't know how to figure out column
        coverage_SNPs_sum=$(awk '{s+=$5}END{print s}' temp_coverage.txt)
        # sum total bases from coverage
        total_bases_sum=$(awk '{s+=$6}END{print s}' temp_coverage.txt)
        # calculate background,  total # SNPs for given MAF / chromosome length
        background=$(echo "$coverage_SNPs_sum / $total_bases_sum" | bc -l)
        for featurefile in "${feature_files[@]}"
            do
                # SNP coverage of feature file
                bedtools coverage -a $featurefile -b $snpfile > temp_featuresnp.txt
                # sum SNPs from coverage

                # sum total bases from coverage

                # calculate enrichment

            done



    done 









# notes: scale = SNPs per base, for each  base we have a probability for that SNP covering that base
# want coverage over whole genome (use whole genome file)
# for each MAF file, calculating same value across genome file
# inner loop subdivide by exons introns and regulator elements
# coverage and calculate sums

