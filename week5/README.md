Week 5 README

### Question 1.1 ####
Since FastQC was originally designed for use on DNA samples, some of its metrics will fail/not be relevant for RNAseq data. For example, RNA is spliced and a single gene can have different introns spliced out in its transcripts, which can result in sequences that map to the same gene but are inconsistent and therefore be flagged as a sequencing error. 

The GC content is the distribution across all samples, and there can be a lot of variation in GC content between RNA transcripts which leads to this quality control metric failing even though there is not actually a problem with the sample quality. 


### Question 1.2 ####
Serine proteases make up most of the overrepresetned sequences in the sample. This makes sense because these enzymes are found in many essential pathways. 


### Question 2 ####
From the Sequence Counts table from FastQC, all samples are below the 45% unique reads cut off and would all be rejected. You would not keep any samples.

Once you increase the min value slider, you can see the blocks of triplicates clearly. This suggests that the replicates are very consistent, which is encouraging for sample data quality and data reproducability. 


### Question 3.3 ####



### Question 3.6 ####



