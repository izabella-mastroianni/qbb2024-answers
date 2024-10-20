#!/usr/bin/bash

# week2 
## need to have script to load in file, to run in terminal do bash week2.sh
# conda activate qb24, bedtools in it, bedtools sort -h | less -S

# 1.4 
# collapsing anything overlapping within same feature, like multiple isoforms

chrfiles=("chr1_exons.bed" "chr1_genes.bed" "chr1_cCREs.bed")
for I in "${chrfiles[@]}";
    do
        bedtools sort -i ${I} > sorted.bed #want to save here as new, sorted file
        # echo ${I}        
        bedtools merge -i sorted.bed > merged_${I} #save as merged files

    done

# 6th file is whole chromosome partitioned into 10kb chunks

# 1.5
# use bedtools to create intron feature file, A - B

bedtools subtract -a merged_chr1_genes.bed -b merged_chr1_exons.bed > introns_chr1.bed

# 1.6 
genome_file="genome_chr1.bed"
subtractingfiles=("merged_chr1_exons.bed" "introns_chr1.bed" "merged_chr1_cCREs.bed")
subresult="$genome_file"

for subtractingfile in "${subtractingfiles[@]}"; do 
    temp_file="${subresult%.bed}_temp.bed"
    bedtools subtract -a "$subresult" -b "$subtractingfile" > "$temp_file"
    subresult="$temp_file"
done     

mv "$subresult" "other_chr1.bed"



