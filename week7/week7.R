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

# extract VST expression matrix and bind to metadata
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()

vsd_df <- bind_cols(metadata_df, vsd_df)

# differential expression of gene WASH7P
m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

# 2.1.1
# WASH7P does not show significant evidence of sex-differential expression because the p value of 0.279 is more than 0.05 threshold needed for significance

# 2.1.2
# differential expression of gene SLC25A47
m2 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

# SLC25A47 does show significant evidence of sex-differential expression because the p value of 0.0257 is less than 0.05 threshold needed for significance. 
# The direction is that it is higher expressed in males because for SEXmale the estimate is positive

# 2.2.1 
# apply DESeq2 to fit regression model 
dds <- DESeq(dds)

# 2.3.1
# extract differential expression results for variable sex
res <- results(dds, name = "SEX_male_vs_female")  %>%
  as_tibble(rownames = "GENE_NAME")

# 2.3.2
# reorder so most diff expressed genes are at the top and remove NAs
res <- res %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# count number of rows with padj value less than 0.1, aka an FDR of 10% or greater 
sex_genes_de <- res %>% 
  filter(padj < 0.1) %>% 
  tally()
# there are 262 genes that are differentially expressed between males and females at a 10% FDR

# 2.3.3
# load in gene locations 
gene_locations_df <- read_delim("gene_locations.txt")

# left join gene locations to results by gene_name and arrange smallest to largest by padj
res_locations <- res %>%
  left_join(gene_locations_df, by = "GENE_NAME") %>%
  arrange(padj)

# most of the genes that are strongly upregulated in males vs females are on the Y chromosome, 
# some are on the X and it takes until the 27th top hit to get an autosomal chromosome
# there are more male upregulated genes at the top of the list because they are on the Y chromosome which only males have


# 2.3.4
# while the p values and padj values are different between the linear regression and the DESeq2 analysis, 
# the two genes WASH7P and SLC25A47 are still called as not significant or significant in both analyses

# 2.4.1
# results of DESeq2 for DTHHRDY differential gene expression
res_DTHHRDY <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes")  %>%
  as_tibble(rownames = "GENE_NAME")

# reorder so most diff expressed genes are at the top and remove NAs
res_DTHHRDY <- res_DTHHRDY %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# number of genes differentially expressed by death classification with 10% FDR
death_genes_de <- res_DTHHRDY %>% 
  filter(padj < 0.1) %>% 
  tally()

# there are 16069 genes that are differentially expressed between death classification at a 10% FDR

# 2.4.2
# it does make sense that more genes are differentially expressed based on death compared to sex 
# since the body upregulates many different genes and processes during death, and the cause of death impacts which are
# while humans are very similar in gene expression regardless of sex, and only a few genes on the sex chromosomes and 
# their target genes make up this difference 


### Exercise 3 ###









