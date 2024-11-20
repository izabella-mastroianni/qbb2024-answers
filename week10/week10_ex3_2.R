# week 10 exercise 3.2

library(tidyverse)
library(dplyr)
library(ggplot2)

# load in csv with data
mean_nucleus_signal <- read.csv("~/qbb2024-answers/week10/mean_nucleus_signal.csv")

# plot each set of data in separate violin plots

# nascent RNA signal
nRNA <- ggplot(mean_nucleus_signal, aes(gene, nascentRNA)) +
  geom_violin(fill = "maroon") +
  ggtitle("Nascent RNA Signal by Gene Knockdown")

ggsave("nRNA_violinplot.pdf", plot = nRNA)
  

# PCNA signal
PCNA <- ggplot(mean_nucleus_signal, aes(gene, PCNA)) +
  geom_violin(fill = "lavender") +
  ggtitle("PCNA Signal by Gene Knockdown")

ggsave("PCNA_violinplot.pdf", plot = PCNA)

# log2 ratio
ratio <- ggplot(mean_nucleus_signal, aes(gene, Ratio)) +
  geom_violin(fill = "lightblue") +
  ggtitle("Log2 Ratio of nascent RNA Signal vs PCNA signal by Gene Knockdown")

ggsave("ratio_violinplot.pdf", plot = PCNA)


