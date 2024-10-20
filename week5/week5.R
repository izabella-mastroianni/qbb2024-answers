# week 5 exercise 3

# set working directory 
setwd("/Users/cmdb/qbb2024-answers/week5")

## 3.1 ##

# load needed libraries
library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(hexbin)
library(ggfortify)

# load salmon gene counts file
salmon_merged_gene_counts = readr::read_tsv('salmon.merged.gene_counts.tsv')

# make gene names the row names
salmon_merged_gene_counts = column_to_rownames(salmon_merged_gene_counts, var="gene_name")

# remove the gene id column 
salmon_merged_gene_counts = salmon_merged_gene_counts %>% 
  select(-gene_id)

# convert data that are numeric to integers so DESeq2 can run properly
salmon_merged_gene_counts = salmon_merged_gene_counts %>%
  mutate_if(is.numeric, as.integer)
  
# keep rows with at least 100 reads
salmon_merged_gene_counts = salmon_merged_gene_counts[rowSums(salmon_merged_gene_counts) > 100,]

# select narrow region samples covering midgut 
narrow = salmon_merged_gene_counts %>%
  select("A1_Rep1":"P2-4_Rep3")


## 3.2 ## 

# create metadata tibble
narrow_metadata = tibble(tissue=as.factor(c("A1", "A1", "A1",
                                            "A2-3", "A2-3", "A2-3",
                                            "Cu", "Cu", "Cu",
                                            "LFC-Fe", "LFC-Fe", "Fe",
                                            "LFC-Fe", "Fe", "Fe",
                                            "P1", "P1", "P1",
                                            "P2-4", "P2-4", "P2-4")),
                         rep=as.factor(c("Rep1", "Rep2", "Rep3",
                                         "Rep1", "Rep2", "Rep3",
                                         "Rep1", "Rep2", "Rep3",
                                         "Rep1", "Rep2", "Rep1",
                                         "Rep3", "Rep2", "Rep3",
                                         "Rep1", "Rep2", "Rep3",
                                         "Rep1", "Rep2", "Rep3")))
# make DESeq2 object
narrow_data = DESeqDataSetFromMatrix(countData=as.matrix(narrow), colData=narrow_metadata, design=~tissue)

# use vst to correct for batch effects
narrow_vst = vst(narrow_data)

# plot data to ensure the vst correction worked to remove the mean and variance relationship
meanSdPlot(assay(narrow_vst))

## 3.3 ##


# make PCA plot and save
pca_narrow_data = plotPCA(narrow_vst, intgroup=c("rep","tissue"), returnData=TRUE)
ggplot(pca_narrow_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5) +
  labs(title = "VST Corrected PCA of Drosophila Midgut Tissues")

# ggsave("~/qbb2024-answers/week5/PCAplot.png", width = 8, height = 6, units = "in")

## 3.4 ##

# averaged replicates
narrow_matrix = as.matrix(assay(narrow_vst))
reps_combined = narrow_matrix[,seq(1, 21, 3)]
reps_combined = reps_combined + narrow_matrix[,seq(2, 21, 3)]
reps_combined = reps_combined + narrow_matrix[,seq(3, 21, 3)]
reps_combined = reps_combined / 3
sds = rowSds(reps_combined)

# only keep replicates with std dev greater than 1
narrow_matrix = narrow_matrix[rowSds(reps_combined) > 1,]

## 3.5 ##

# set seed so clustering is reproducible
set.seed(42)

# 12 k-means clusters
k=kmeans(narrow_matrix, centers=12)$cluster

# find ordering of samples to put into clusters
ordering = order(k)

# reorder genes 
k = k[ordering]

# plot heatmap and save
png("heatmap.png", width = 800, height = 600)
heatmap(narrow_matrix[ordering,], Rowv=NA, Colv=NA, RowSideColors=RColorBrewer::brewer.pal(12,"Paired")[k])
dev.off()

## 3.6 ##

# GO enrichment analysis 

# select gene names from a specific cluster
c1_genes = rownames(narrow_matrix[k == 1,])

# save gene names to txt file
write.table(c1_genes, "cluster1_genes.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)



