# Load necessary libraries
library(dplyr)

# Load the datasets
gdc_train <- read.csv("gdc_train.csv")
cl_data2B <- read.csv("cl_data2B.csv")

# Create a subset of cl_data2B where Sample_ID matches those in gdc_train
# This will be your train90 data frame
train90 <- cl_data2B %>%
  semi_join(gdc_train, by = "Sample_ID")

# Create another subset of cl_data2B where Sample_IDs do not match those in gdc_train
# This will be your test10 data frame
test10 <- cl_data2B %>%
  anti_join(gdc_train, by = "Sample_ID")

# Optionally, you can write these subsets to new CSV files
write.csv(train90, "train90.csv", row.names = FALSE)
write.csv(test10, "test10.csv", row.names = FALSE)
