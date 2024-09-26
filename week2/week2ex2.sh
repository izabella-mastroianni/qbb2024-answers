#!/usr/bin/bash


# week 2 exercise 2

# 2.1
# 4 feature files are merged_exons, introns, merged_cCREs, and other
# note which part of pseudocode each part of following code refers to

# create output results file called snp_counts.txt
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

# loop through each possible MAF file
snp_files=("chr1_snps_0.1.bed" "chr1_snps_0.2.bed" "chr1_snps_0.3.bed" "chr1_snps_0.4.bed" "chr1_snps_0.5.bed") 
feature_files=("chr1_cCREs.bed" "chr1_exons2.bed" "chr1_introns.bed" "chr1_other.bed") 
genome_file=("genome_chr1.bed")


 # make result file
echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt


# per snp file (aka MAF value) per feature file
for snpfile in "${snp_files[@]}"
    do
        # SNP coverage of the whole chromosome 
        bedtools coverage -a $genome_file -b $snpfile > temp_coverage.txt
        # sum SNPs from coverage - find column from ppt labels and looked in exon file
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
                feature_coverage_SNPs_sum=$(awk '{s+=$5}END{print s}' temp_featuresnp.txt)
                # sum total bases from coverage
                feature_total_bases_sum=$(awk '{s+=$6}END{print s}' temp_featuresnp.txt)
                # calculate enrichment
                ratio=$(echo "$feature_coverage_SNPs_sum / $feature_total_bases_sum" | bc -l)
                enrichment=$(echo "$ratio / $background" | bc -l)
                # save to result file
                echo -e "$snpfile\t$featurefile\t$enrichment" >> snp_counts.txt 


            done

    done 



# echo "Coverage SNPs: $feature_coverage_SNPs_sum"
# echo "Total Bases: $feature_total_bases_sum"






# notes: scale = SNPs per base, for each  base we have a probability for that SNP covering that base
# want coverage over whole genome (use whole genome file)
# for each MAF file, calculating same value across genome file
# inner loop subdivide by exons introns and regulator elements
# coverage and calculate sums

