# plot for week 3 exercise 3.2

library(ggplot2)


# load in allele frequency data from AF.txt file and the dp data from DP.txt
af_data <- read.table("/Users/cmdb/qbb2024-answers/week3/AF.txt", header = TRUE)
dp_data <- read.table("/Users/cmdb/qbb2024-answers/week3/DP.txt", header = TRUE)


# make sure AF values are numeric 
af_data$Allele_Frequency <- as.numeric(as.character(af_data$Allele_Frequency))


# histogram of the af data
ggplot(af_data, mapping = aes(x = Allele_Frequency)) +
  geom_histogram(bins = 11, fill = "orange", color = "black") +
  aes(y = after_stat(count)/sum(after_stat(count))) + 
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Distribution of Allele Frequencies of Variants",
       y = "Percent of Variants",
       x = "Allele Frequency") +
  theme_minimal()

# save the plot
ggsave("/Users/cmdb/qbb2024-answers/week3/af_plot.png")


## Question 3.1 ##
# The plot is approximately a normal distribution. This makes sense due to selection pressure and genetic drift impacting the presence of different alleles overtime depending on mutations and environmental conditions for example. 


# histogram of dp data
ggplot(dp_data, mapping = aes(x = Read_Depth)) +
  # bins of 21 mean some high depth are cut off, 1094 rows from txt file says error
  geom_histogram(bins = 21, fill = "pink", color = "black") +
  xlim(0, 20) +
  aes(y = after_stat(count)/sum(after_stat(count))) + 
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Distribution of Read Depth for Variants",
       y = "Percent of Variants",
       x = "Read Depth") +
  theme_minimal()

# save the plot
ggsave("/Users/cmdb/qbb2024-answers/week3/dp_plot.png")

## Question 3.2 ##
# This distribution looks as expected, almost a normal distribution. The peak around a read depth of 4 which was what we estimated earlier. 





