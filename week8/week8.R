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
# The three most interesting column names to me are tissue, n_genes, and age because combining information from those will give you insight into fly aging.

# UMAP plot
plotReducedDim(gut, "X_umap", colour_by= "broad_annotation")











