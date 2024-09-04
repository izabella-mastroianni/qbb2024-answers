# day 2 lunch answers

## Answer 1

 cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c

 There are 19618 protein coding genes

 I would like to learn more about translated processed pseudogenes because since they are processed, I am curious if they have any activity, even though they are classified as pseudogenes

## Answer 2

cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort -n

ENSG00000168036 is the ensembl_gene_id with the most go_ids

grep -w "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k3 | > q2.txt

Based on the GO terms, I think the gene is a transcriptional factor and regulator in many different cell types, including T cells and adherens.

## Answer 1 (section 2)

grep -w "IG_._gene" genes.gtf | cut -f 1 | uniq -c

Number of IG genes on each of the following chromosomes: chr2 = 52, chr14 = 91, chr15 = 16, chr16 = 6, chr21 = 1, chr22 = 48

grep -w "IG_._pseudogene" genes.gtf | cut -f 1 | uniq -c

Number of IG pseudogenes on each of the following chromosomes: chr1 = 1, chr 2 = 45, chr8 = 1, chr9 = 5, chr10 = 1, chr14 = 83, chr15 = 6, chr16 = 8, chr18 = 1, chr22 = 48

Pseudogenes are spread across more chromosomes than IG genes. 

