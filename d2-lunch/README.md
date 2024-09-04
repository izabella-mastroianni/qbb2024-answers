# day 2 lunch answers

## Answer 1

 cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c

 There are 19618 protein coding genes

 I would like to learn more about translated processed pseudogenes because since they are processed, I am curious if they have any activity, even though they are classified as pseudogenes

## Answer 2

cut -f 1 hg38-gene-metadata-go.tsv | uniq -c | sort -n

ENSG00000168036 is the ensembl_gene_id with the most go_ids

grep -w "ENSG00000168036" hg38-gene-metadata-go.tsv | sort -k3 | > q2.txt

Based on the GO terms, I think the gene is a transcriptional regulator and part of the MAPK cascade in T cells.





