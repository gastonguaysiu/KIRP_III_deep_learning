# Load necessary libraries
library(readr)
library(dplyr)
library(caret)
library(gbm)
library(ggplot2)

# Read and prepare data
matched_subset_df <- read_csv('cl_data2B.csv') %>% 
  rename(Sample_ID = 1) %>% 
  as.data.frame()

ref4_df <- read_csv('gdc_join_clinical2.csv') %>% 
  filter(Sample.Type == "Primary Tumor")

live_score_constant <- median(ref4_df$days_to_death, na.rm = TRUE)

# Data transformations
ref4_df <- ref4_df %>%
  mutate(
    Binary_Survival = ifelse(vital_status == 'Dead', 0, 1),
    live_score = ifelse(vital_status == 'Dead', 1, days_to_last_follow_up / live_score_constant)
  )

# Merge methylation data with survival labels using 'Sample_ID' as the key
ref4_df2 <- ref4_df %>% filter(live_score > 0 | vital_status == "Dead")
merged_data <- merge(matched_subset_df, ref4_df2[, c('Sample_ID', 'Binary_Survival', 'live_score')], by = 'Sample_ID')

# Transform to numeric data and remove rows with NA
numeric_data <- merged_data %>%
  select(-Sample_ID) %>%
  mutate(across(everything(), as.numeric)) %>%
  filter(!is.na(Binary_Survival) & !is.na(live_score))

# Set row names to 'Sample_ID' from the merged data
rownames(numeric_data) <- merged_data$Sample_ID

# Split data into training and test sets, excluding live_score from training and testing
set.seed(42)
index <- createDataPartition(numeric_data$Binary_Survival, p = .8, list = FALSE)
trainSet <- numeric_data[index, -which(names(numeric_data) == "live_score")]
testSet <- numeric_data[-index, -which(names(numeric_data) == "live_score")]

# Train the gradient boosting model
gbm_model <- gbm(
  Binary_Survival ~ ., 
  data = trainSet, 
  distribution = "bernoulli",
  n.trees = 500, 
  interaction.depth = 3,
  shrinkage = 0.01,
  cv.folds = 5,
  verbose = FALSE
)

# Predict on test set
predicted_probabilities <- predict(gbm_model, testSet, n.trees = gbm_model$n.trees, type = "response")
predicted_classes <- ifelse(predicted_probabilities > 0.5, 1, 0)

# Create a dataframe with actual vs predicted values and the sample names
temp <- data.frame(
  Sample_ID = rownames(testSet),
  Actual = testSet$Binary_Survival,
  Predicted_Class = predicted_classes
)

validation_results_df <- merge(temp, ref4_df, by = 'Sample_ID')

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

# Assess feature importance
feature_importance <- summary(gbm_model)
feature_importance$rel.inf <- feature_importance$rel.inf/100

# Generate the sideways bar plot with the largest number at the top
ggplot(feature_importance, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat = 'identity', fill = 'blue') +
  coord_flip() +  # Make the bar plot horizontal
  labs(x = '', y = 'Relative Importance (%)') +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  geom_text(aes(label = sprintf("%.1f%%", rel.inf * 100), hjust = -0.1)) +  # Adjust labels to show correct percentage
  theme_minimal()

# Display the plot
# ggsave("feature_importance_plot.png", width = 10, height = 6)

print(feature_importance)
