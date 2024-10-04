# week 2 exercise 2.1
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)

df <- readr::read_tsv("~/qbb2024-answers/week2/snp_counts.txt")

df2 <- df %>%
  mutate(Feature=gsub("chr1_", "",Feature)) %>%
  mutate(Feature=gsub(".bed", "", Feature)) %>%
  mutate(Feature=gsub("2", "", Feature)) %>%
  mutate(MAF=gsub(".bed", "", MAF)) %>%
  mutate(MAF=gsub("chr1_snps_", "", MAF))
  

df2_transformed <- df2

df2_transformed$Enrichment <- log2(df2_transformed$Enrichment)

ggplot(df2_transformed, mapping=aes(x=MAF, y=Enrichment, color = Feature, group = Feature)) +
  geom_point() +
  geom_line() +
  ylab("Log2 Enrichment") +
  ggtitle("SNP Enrichment by Feature")

# save plot
ggsave("/Users/cmdb/qbb2024-answers/week2/snp_enrichment_plot.png")


