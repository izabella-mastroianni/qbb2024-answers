# week 2 exercise 2.1
library(tidyverse)
library(ggplot2)

df <- readr::read_tsv("~/qbb2024-answers/week2/snp_counts.txt")

df2 <- df %>%
  mutate(Feature=gsub("chr1_", "",Feature)) %>%
  mutate(Feature=gsub(".bed", "", Feature)) %>%
  mutate(Feature=gsub("2", "", Feature)) %>%
  mutate(MAF=gsub(".bed", "", MAF)) %>%
  mutate(MAF=gsub("chr1_snps_", "", MAF))
  


ggplot(df2, mapping=aes(x=MAF, y=Enrichment, color = Feature, group = Feature)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(trans="log2") +
  ggtitle("SNP Enrichment by Feature")


