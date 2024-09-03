library("tidyverse")
df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )

df %>%
  group_by(SUBJECT) %>%
  summarize(n()) %>%
  arrange(-`n()`)

#subjects K-562 and GTEX-NPJ8 have the most samples

df %>%
  group_by(SUBJECT) %>%
  summarize(n()) %>%
  arrange(`n()`)

#subjects GTEX-1JMI6 and GTEX-1PAR6 have the least samples

df %>%
  group_by(SMTSD) %>%
  summarize(n()) %>%
  arrange(-`n()`)

#tissue types Whole Blood and Muscle - Skeletal have the most samples.  Blood is very easy to extract and analyze, and muscle is very plentiful and easy to collect, which is likely why they are the most common sampled tissues.

df %>%
  group_by(SMTSD) %>%
  summarize(n()) %>%
  arrange(`n()`)

#tissue types Kidney - Medulla and Cervix - Ectocervix and Fallopian Tubs have the least samples, likely because not all subjects have a cervix and fallopian tubes, and the tissues may be difficult to sample or not often needed

df_npj8 <- filter(df, SUBJECT == "GTEX-NPJ8")

df_npj8

df_npj8 %>%
  group_by(SMTSD)%>%
  summarize(n()) %>%
  arrange(-`n()`)

#Whole Blood has the most samples. The difference between the samples is the isolation protocol and the analysis assay performed.

SMATSSCR <- df %>%
  filter( !is.na(SMATSSCR))%>%
  group_by(SUBJECT)%>%
  summarize(meanSM=mean(SMATSSCR))

sum(SMATSSCR$meanSM == 0)

#15 subjects have mean SMATSSCR score of 0. Other observations you could make are the standard deviation, how the mean differs by tissue type or other metadata. To present this information in a report, you could make a graph to show the distribution or a table with how many subjects have each SMATSSCR


  

