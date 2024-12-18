# week 8 qb lab homework 
# looking at already analyzed data by fly cell atlas

setwd("~/qbb2024-answers/week8")

### Exercise 1 ###

# load in libraries 
library(zellkonverter)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)

# load in and inspect data - gut data downloaded on web browser from flycellatlas.org
gut <- readH5AD("~/qbb2024-answers/week8/v2_fca_biohub_gut_10x_raw.h5ad")

# change assay name from X to counts
assayNames(gut) <- "counts"

# use logNorm to normalize counts
gut <- logNormCounts(gut)

## Question 1 

# number of genes and cells
dimensions <- dim(gut)
num_rows <- dimensions[1]
num_cols <- dimensions[2]

# There are 13407 genes that are being quantitated (the number of rows)
# There are 11788 cells in the dataset

# names of dimension reduction datasets present
reducedDimNames(gut)
# The dimension reduction datasets present are PCA, TSNE, and UMAP


## Question 2

# number of columns in colData(gut)
num_cols_colData <- ncol(colData(gut))
# There are 39 columns in colData(gut)

# to see column names in colData(gut)
colnames(colData(gut))
# The three most interesting column names to me are tissue, n_genes, and age 
# because combining information from those will give you insight into fly aging.

# UMAP plot
UMAP_broad_annotation_plot <- plotReducedDim(gut, "X_umap", colour_by= "broad_annotation")

ggsave("UMAP_broad_annotation_plot.png", plot = UMAP_broad_annotation_plot, width = 8, height = 6, dpi = 300)


### Exercise 2 ###

# sum expression of each gene across all cells
genecounts <- rowSums(assay(gut))


## Question 3 

# get summary stats of genecounts
summary(genecounts)
# the mean genecount is 3185 and the median is 254. 
# from this, you can conclude that there is a wide distribution of gene counts 
# since the middle value is so different from the average

# sort genecounts from highest to lowest
head(sort(genecounts, decreasing = TRUE))
# the three genes with the highest expression are lncRNA:Hsromega, pre-rRNA:CR45845, and lncRNA:roX1.
# they are all non-coding RNAs, either long or ribosomal. 


## Question 4a

# cell counts using assay to make gut a matrix that colSums can process
cellcounts <- colSums(assay(gut))

# make histogram of cell counts
hist(cellcounts)

# summary stats of cellcounts
summary(cellcounts)
# the mean number of counts per cell is 3622

# for cells with very high total counts (above 10k), this could be due to the biological
# function of the cell if its job is to make a lot of a protein,
# or it could be due to technical errors and quality control issue, like if a cell lysed
# multiple cells are counted as one, or if contaminants not in the cell are counted as in the cell. 


## Question 4b

# number of genes detected in cell vector
celldetected <- colSums(assay(gut)>0)

# histogram of cells detected
hist(celldetected)

# summary stats of celldetected
summary(celldetected)
# the mean number of genes detected per cell in 1059

# 1059/13407, this represents about 8% of the total number of genes

# mito genes vector
mito <- grep("^mt:", rownames(gut), value = TRUE)

# dataframe for mito genes
df <- perCellQCMetrics(gut, subsets = list(Mito = mito))

# convert df to data.frame
df <- as.data.frame(df)
summary(df)

# add metrics to cell data
colData(gut) <- cbind(colData(gut), df)

# check above worked 
print(colnames(colData(gut)))


## Question 5

# colData as df for plotting
col_data_df <- as.data.frame(colData(gut))

# plot subsets_Mito_percent vs broad_annotation
mito_vs_broad_plot <- ggplot(col_data_df, aes(x = broad_annotation, y = subsets_Mito_percent)) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=90)) +
  labs(x = "Broad Tissue Annotation", y = "Mitochondrial Percent",
  title = "Mitochondrial Gene Expression in Broad Tissue Categories") +
  theme_minimal()

ggsave("mito_vs_broad_plot.png", plot = mito_vs_broad_plot, width = 8, height = 6, dpi = 300)

# epithelial cells and gut cells have higher percentages of mitochondrial gene expression
# this is possibly because there is a lot of rapid turnover of these cell types


### Exercise 3 ###

## Question 6a

# subset cells annotated as "epithelial cell"
coi <- colData(gut)$broad_annotation == "epithelial cell"

# make new SingleCellExperiment object subsetting gut by coi
epi <- gut[, coi]

# make epi a dataframe for plotting ease 
epi_df <- as.data.frame(reducedDim(epi, "X_umap"))
epi_df$annotation <- colData(epi)$broad_annotation

# UMAP plot of epi 
epi_umap_plot <- ggplot(epi_df, aes(x = V1, y = V2, color = annotation)) +
  geom_point() +
  labs(title = "UMAP of Epithelial Cells", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("epi_umap_plot.png", plot = epi_umap_plot, width = 8, height = 6, dpi = 300)

# marker genes in anterior midgut
# pairwise comparisons of all annotation categories
marker.info <- scoreMarkers(epi, colData(epi)$annotation)

# identify top marker genes
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])


## Question 6b

# the top six marker genes in the fly anterior midgut are Mal-A6, Men-b, vnd, betaTry, Mal-A1, Nhe2. 
# they are involved in metabolism of glucose and other sugars

# plot expression of top marker genes
top_marker_plot <- plotExpression(epi, c("Mal-A6","Men-b","vnd","betaTry","Mal-A1","Nhe2"), x ="annotation") + 
  labs(x = "Annotation", y = "Expression Log Fold Change",
       title = "Anterior Midgut Expression Log Fold Change") +
  theme(axis.text.x=element_text(angle=90))

ggsave("top_marker_genes_anterior_midgut.png", plot = top_marker_plot)


## Question 7

# analyze somatic percursor cells
# subset somatic precursor cells 
spc_all <- colData(gut)$broad_annotation == "somatic precursor cell"
spc <- gut[, spc_all]

# log normalize counts so scoreMarkers runs 
spc <- logNormCounts(spc)  

# marker genes for precursor cells
spc_marker_info <- scoreMarkers(spc, colData(spc)$annotation)

# selected genes
selected_spc <- spc_marker_info[["intestinal stem cell"]]

# order genes 
ordered_spc <- selected_spc[order(selected_spc$mean.AUC, decreasing=TRUE),]

# create genes of interest vector with names of top 6 marker genes
spc_goi = rownames(ordered_spc)[1:6]

# plot top 6 genes
spc_goi_plot <- plotExpression(spc, features = spc_goi, x = "annotation") +
  labs(x= "Annotation",y="Expression Log Fold Change",
       title = "Log Fold Change of Expression of Top Six Somatic Precursor Cell Marker Genes") +
  theme(axis.text.x=element_text(angle=90)) 

ggsave("top_six_genes_spc.png", plot = spc_goi_plot)

# The two cell types with the most similar expression patterns are enteroblasts and intestinal stem cells
# DI is the marker that looks most specific for intestinal stem cells

