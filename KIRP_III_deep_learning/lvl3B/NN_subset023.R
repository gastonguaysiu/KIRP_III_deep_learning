# Load necessary libraries
library(readr)
library(dplyr)
library(caret)
library(reticulate)
library(keras)
library(tensorflow)

# Environment settings for TensorFlow to use CPU
Sys.setenv(TF_FORCE_GPU_ALLOW_GROWTH = "false", CUDA_VISIBLE_DEVICES = "")

# Use Python executable
use_python("/usr/bin/python3", required = TRUE)

# Set seed for reproducibility
tensorflow::set_random_seed(42)

# Read and prepare data
matched_subset_df <- read_csv('cl_data2B.csv') %>% rename(Sample_ID = 1) %>% as.data.frame()
ref4_df <- read_csv('gdc_join_clinical2.csv') %>% filter(Sample.Type == "Primary Tumor")
live_score_constant <- median(ref4_df$days_to_death, na.rm = TRUE)

# Data transformations
ref4_df <- ref4_df %>%
  mutate(
    Binary_Survival = ifelse(vital_status == 'Dead', 0, 1),
    live_score = ifelse(vital_status == 'Dead', 1, days_to_last_follow_up / live_score_constant)
  )

# Preparing datasets
prepare_dataset <- function(df, threshold) {
  df <- df %>%
    filter(live_score > threshold | vital_status == "Dead") %>%
    merge(matched_subset_df, by = 'Sample_ID') %>%
    select(-Sample_ID) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(!is.na(Binary_Survival) & !is.na(live_score))
  rownames(df) <- df$Sample_ID
  return(df)
}

ref4_df2 <- ref4_df %>% filter(live_score > 0 | vital_status == "Dead")

# Merge methylation data with survival labels using 'Sample_ID' as the key
merged_data <- merge(matched_subset_df, ref4_df2[, c('Sample_ID', 'Binary_Survival', 'live_score')], by = 'Sample_ID')

# Transform to numeric data
numeric_data <- merged_data %>%
  select(-Sample_ID) %>%
  mutate(across(everything(), as.numeric))

# Set row names to 'Sample_ID' from the merged data
rownames(numeric_data) <- merged_data$Sample_ID

# Split data into training and test sets
set.seed(42)
index <- createDataPartition(numeric_data$Binary_Survival, p = .8, list = FALSE)
trainSet <- numeric_data[index,]
testSet <- numeric_data[-index,]

# Model definition and compilation
model <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = 'tanh', input_shape = c(ncol(trainSet) - 2)) %>%
  layer_dense(units = 32, activation = 'leaky_relu') %>%
  layer_dense(units = 1, activation = 'sigmoid')

model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# Fit model
history <- model %>% fit(
  x = as.matrix(trainSet[,1:(ncol(trainSet) - 2)]),
  y = trainSet$Binary_Survival,
  epochs = 16,
  batch_size = 32,
  validation_split = 0.2
)

# Model evaluation and predictions
performance <- model %>% evaluate(
  x = as.matrix(testSet[,1:(ncol(testSet) - 2)]),
  y = testSet$Binary_Survival
)
print(performance)

# Predict on new data (excluding both the target and error value columns)
predictions <- model %>% predict(as.matrix(testSet[,1:(ncol(testSet) - 2)])) # Exclude the last two columns

# Generate predictions and convert to classes
predicted_probabilities <- model %>% predict(as.matrix(testSet[,1:(ncol(testSet) - 2)]))
predicted_classes <- ifelse(predicted_probabilities > 0.5, 1, 0)



# Create a dataframe with actual vs predicted values and the sample names
temp <- data.frame(
  Sample_ID = rownames(testSet),
  Actual = testSet$Binary_Survival,
  Predicted_Class = predicted_classes
)

validation_results_df = merge(temp, ref4_df, by = 'Sample_ID')

dead_true_tot <- sum(validation_results_df$Predicted_Class == 0 & validation_results_df$vital_status == "Dead")
denominator_recall_tot <- sum(validation_results_df$vital_status == 'Dead')
recall_tot = dead_true_tot / denominator_recall_tot

live_Alpha <- sum(
  validation_results_df$Predicted_Class == 0 & 
    validation_results_df$vital_status == "Alive",
  na.rm = TRUE)
live_weight_tot <- sum(validation_results_df$live_score[
  validation_results_df$Predicted_Class == 0 & 
    validation_results_df$vital_status == "Alive"
], na.rm = TRUE)

precision_alpha = dead_true_tot / (dead_true_tot + live_Alpha)
precision_tot = dead_true_tot / (dead_true_tot + live_weight_tot)
f1_Alpha = if_else((precision_alpha + recall_tot) == 0, 0, 2 * (precision_alpha * recall_tot) / (precision_alpha + recall_tot))
f1_Beta = if_else((precision_tot + recall_tot) == 0, 0, 2 * (precision_tot * recall_tot) / (precision_tot + recall_tot))

print(paste("F1_aplha:", f1_Alpha, "F1_Beta:", f1_Beta))
