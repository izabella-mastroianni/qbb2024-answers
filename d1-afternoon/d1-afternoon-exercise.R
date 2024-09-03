#Q1
library(tidyverse)
df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

#Q2
df

glimpse(df)

#Q3
RNAseqOnly <- df %>%
  filter(SMGEBTCHT == "TruSeq.v1")

glimpse(RNAseqOnly)

#Q4
ggplot(data=RNAseqOnly, 
       mapping = aes(x = SMTSD)) +
  geom_histogram(stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue Type") +
  ylab("Number of Samples") +
  ggtitle("Number of Samples for each Tissue Type")

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q4.pdf")

  
#Q5
ggplot(data=RNAseqOnly, 
       mapping = aes(x = SMRIN)) +
  geom_histogram(stat = "count")+
  xlab("RNA Integrity Number") +
  ylab("Number of Samples") +
  ggtitle("Distribution of RNA Integrity Numbers")

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q5.pdf")
 
#Q5: I used a histogram plot to best visualize the distribution. The shape is unimodal because there is one big peak and a small spike at the end.

#Q6
ggplot(data=RNAseqOnly, 
       mapping = aes(x = SMRIN, y = SMTSD)) +
  geom_boxplot()+
  xlab("RNA Integrity Number") +
  ylab("Tissue Type") +
  ggtitle("Distribution of RNA Integrity Numbers by Tissue Type")

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q6.pdf")

#Q6: The RIN numbers vary by tissue type, with some of the brain and artery samples on average having worse RINs, possibly because RNA is degraded more often in those cells. The cell culture lines and cancer samples are outliers with very high RINs, likely because they produce so much RNA. 

#Q7 
ggplot(data=RNAseqOnly, 
       mapping = aes(x = SMGNSDTC, y = SMTSD)) +
  geom_boxplot()+
  xlab("Genes per Sample") +
  ylab("Tissue Type") +
  ggtitle("Distribution of Number of Genes by Tissue Type")

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q7.pdf")

#Q7: Most of the tissues mean number of genes per sample is around 22,000, but the testis is an outlier with significantly more genes expressed - averaging over 32,000, possibly because the widespread gene transcription is used to maintain genome integrity. Whole blood expresses the least number of genes, likely because red blood cells express so few genes and eject their nuclei.  

#Q8 
ggplot(data=RNAseqOnly, 
      mapping = aes(x = SMRIN, y = SMTSISCH)) +
  geom_point(size = 0.5, alpha = 0.5)+
  facet_wrap(SMTSD~.)+
  geom_smooth(method = "lm")+
  xlab("RNA Integrity Number") +
  ylab("Ischemic Time") +
  ggtitle("Distribution of RIN by Ischemic Time")

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q8.pdf")

#Q8: Most distributions have fewer samples with high RIN numbers/number of samples decrease as RIN increases. Some tissue types are relatively flat/consistent distributions, especially in brain samples.

#Q9
ggplot(data=RNAseqOnly, 
       mapping = aes(x = SMRIN, y = SMTSISCH)) +
  geom_point(size = 0.5, alpha = 0.5, 
             mapping = aes(color = SMATSSCR))+
  facet_wrap(SMTSD~.)+
  geom_smooth(method = "lm")+
  xlab("RNA Integrity Number") +
  ylab("Ischemic Time") +
  labs(color="Autolysis") +
  ggtitle("Distribution of RIN by Ischemic Time and Autolysis")
  

ggsave(filename = "~/qbb2024-answers/d1-afternoon/Q9.pdf")


