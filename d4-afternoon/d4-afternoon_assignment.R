library(tidyverse)
library(broom)

# exercise 1 
# read in the data
dnm <- read_csv(file="~/qbb2024-answers/d4-afternoon/aau1043_dnm.csv")

ages <- read_csv(file="~/qbb2024-answers/d4-afternoon/aau1043_parental_age.csv")

# count paternally vs.  maternally inherited DNMs
dnm_summary <- dnm%>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE), 
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))

dnm_by_paternal_age <- left_join(dnm_summary, ages, by ="Proband_id")


# exercise 2
# the count of maternal de novo mutations vs. maternal age
ggplot(data = dnm_by_paternal_age, 
       mapping = aes(x = Mother_age, y = n_maternal_dnm, )) +
  geom_point() +
  xlab("Maternal Age") +
  ylab("Maternal DNMs") +
  ggtitle("Maternal DNMs vs Maternal Age")

# the count of paternal de novo mutations vs. paternal age
ggplot(data = dnm_by_paternal_age, 
       mapping = aes(x = Father_age, y = n_paternal_dnm, )) +
  geom_point() +
  xlab("Paternal Age") +
  ylab("Paternal DNMs") +
  ggtitle("Paternal DNMs vs Paternal Age")


# lm for maternal DNMs vs age (2.2)
lm_maternal <- lm(data = dnm_by_paternal_age, 
                  formula = n_maternal_dnm ~ 1 + Mother_age) %>% 
  summary()
print(lm_maternal)

# this means that increases in maternal age only increase the number of DNMs by a factor of 0.3776, which is small. this matches the plot of the data
# this relationship is significant because the p value is 2e-16, which means that increases in maternal age significantly increase the number of maternal DNMs

# lm for paternal DNMs vs age (2.3)
lm_paternal <- lm(data = dnm_by_paternal_age, 
                  formula = n_paternal_dnm ~ 1 + Father_age) %>% 
  summary()
print(lm_paternal)

# this means that increases in paternal age increase the number of DNMs by a factor of 1.3, which is larger than maternal. this matches the plot of the data
# this relationship is significant because the p value is 2e-16, which means that increases in paternal age significantly increase the number of paternal DNMs

# 2.4
paternal_DNMs_50.5y_father <- 10.32632 + 1.35384*50.5
print(paternal_DNMs_50.5y_father)
# the predicted number of paternal probands is 78.69524, using b0 + b1*x1 where x1 is father's age


# 2.5
ggplot(data = dnm_by_paternal_age) +
  geom_histogram(mapping = aes(x = n_maternal_dnm), alpha = 0.5, fill = "green") +
  geom_histogram(mapping = aes(x = n_paternal_dnm), alpha = 0.5, fill = "blue")


# 2.6

t.test(dnm_by_paternal_age$n_paternal_dnm, dnm_by_paternal_age$n_maternal_dnm, paired = TRUE)
#I used a paired t test because we are just comparing two items. 
#The p-value < 2.2e-16 is statistically significant, meaning that there is a statistically significant difference between paternal and maternal DNMs 





