# Load cl_data - replace this with your actual data loading code if needed
cl_data <- read.csv("cl_data.csv")

cl_data2 <- cl_data

# Modify miRNA_cl column
cl_data$miRNA_cl[cl_data$miRNA_cl == 1] <- "I"
cl_data$miRNA_cl[cl_data$miRNA_cl == 2] <- "B"

# Modify RNAseq_cl column
cl_data$RNAseq_cl[cl_data$RNAseq_cl == 3] <- "B"
cl_data$RNAseq_cl[cl_data$RNAseq_cl == 1] <- "W"



# Replace "I" with 2 in both columns
cl_data$miRNA_cl[cl_data$miRNA_cl == "W"] <- 3
cl_data$RNAseq_cl[cl_data$RNAseq_cl == "W"] <- 3

# Replace "I" with 2 in both columns
cl_data$miRNA_cl[cl_data$miRNA_cl == "I"] <- 2
cl_data$RNAseq_cl[cl_data$RNAseq_cl == "I"] <- 2

# Replace "B" with 1 in both columns
cl_data$miRNA_cl[cl_data$miRNA_cl == "B"] <- 1
cl_data$RNAseq_cl[cl_data$RNAseq_cl == "B"] <- 1

# Optional: View the modified dataset
print(head(cl_data))

write.csv(cl_data,"cl_data2.csv")

