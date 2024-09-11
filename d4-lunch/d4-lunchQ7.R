# Q7 
library(tidyverse)
library(ggplot2)

#paste combines strings, 0 equals no space

dicts_expr <- read_tsv("dicts_expr.tsv") %>% 
  mutate(Tissue_Gene=paste0(Tissue, " ", GeneID)) %>% 
  mutate(log2Expr = log2(Expr))



ggplot(data = dicts_expr, mapping = aes(y= Expr, x=Tissue_Gene)) +
  geom_violin(mapping = aes())+ 
  scale_y_continuous(trans="log10") +
  ylab("Gene Expression") +
  xlab("Tissue and Gene ID") +
  ggtitle("Gene Expression by Tissue and Gene ID")+
  coord_flip() 
  
# No, i am not surprised by the results because high tissue specificity = lower expression in other tissues for a given gene
# The pancreas has the most variety in mean expression, likely partly because it has more genes than the other tissues. Most of the genes for a given tissue have similar mean expression values. The stomach, small intestine, and pituitary have longer bulbs on the violin plot.
# Maybe tissues with low variability need the genes on at a constant level whereas the others have the genes expressed at different levels depending on environmental signals. 

