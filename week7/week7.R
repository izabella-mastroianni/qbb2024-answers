library(tidyverse) # read in data
library(broom) # takes output of stat models and format in table
library(DESeq2) # stat tests for differential expression

# set your working directory to where your data and output will be stored
setwd("~/qbb2024-answers/week7")

### Exercise 1 ###

# 1.1.1
# load in counts data
counts_df <- read_delim("gtex_whole_blood_counts_downsample.txt")

# load in metadata
metadata_df <- read_delim("gtex_metadata_downsample.txt")

# 1.1.2
# move gene names to row names for counts
counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")

# 1.1.3
# move subject IDs to row names for metadata
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")

# 1.1.4
# look at first five rows of counts to make sure it loaded correctly 
counts_df[1:5,]
# look at first five rows of metadata to make sure it loaded correctly 
metadata_df[1:5,]


# 1.2.1
# check that the columns of the counts are identical and in the same order as the rows of the metadata
colnames(counts_df) == rownames(metadata_df)

# table to count number of TRUEs and FALSEs
table(colnames(counts_df) == rownames(metadata_df))

# 1.2.2
# load data into DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                               colData = metadata_df,
                              # column headers from metadata file that aren't sample id
                               design = ~ SEX + AGE + DTHHRDY)

# 1.3.1
# vst to normalize for PCA analysis, vsd = variance-stabilized data
vsd <- vst(dds)

# 1.3.2
# apply and plot principal components, intgroup = factor to color points by
# by sex PCA plot
png("sex_PCAplot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "SEX") 
dev.off()

# by age PCA plot
png("age_PCAplot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "AGE") 
dev.off()

# by dthhrdy PCA plot
png("dthhrdy_PCAplot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = "DTHHRDY") 
dev.off()

# 1.3.3
# PC1 explains 48% of the variance, PC2 explains 7% of the variance.
# PC1 seems to be associated with the subject-level variable of DTHHRDY, since most of the fast death of natural causes is are on the left and most of the ventilator cases are on the right of the graph. 
# For PC1, age and sex appear to be mostly even mixes and are instead probably partly contributing to PC2, but since it only has 7% variance it is harder to tell visually


### Exercise 2 ###




