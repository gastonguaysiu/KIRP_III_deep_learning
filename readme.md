# Hybrid Predictive Model Implementation

## Overview
This document elaborates on the enhancement of my existing predictive model, which utilizes machine learning (ML) techniques for analyzing differential methylation data. The novel addition incorporates a neural network layer designed to assimilate and process multi-omic data, including GpG, RN-seq, and miRNA, thereby augmenting the model's predictive capabilities.

## Data Handling

### lvl4B
The data segmentation remains consistent with the previous iterations of the machine learning model:
- **Training Data**: Continues to be 90% of the sample set, aligned with the prior phase's configuration.
- **Testing Data**: Retains the 10% allocation, earmarked for subsequent model validation.

### lvl3B
This stage introduces a new data partitioning scheme, maintaining a blind approach towards the categorization outcomes of the last prediction round:
- **Training Data**: Now adjusted to 80% of the samples, selected anew in a random manner.
- **Testing Data**: Expanded to encompass the remaining 20%, designated for validation purposes.

## Predictive ML Workflow

The neural network architecture for this model (NN_subset023.R & NN_subset033.R) comprises three layers: an initial 32-node layer using the tanh activation function, followed by a 32-node layer with a leaky ReLU activation function, and concluding with a single-node layer employing a sigmoid function to categorize patients into two distinct groups.

An alternative modeling approach  explored is the gradient boosting machine (GBM) model (NN_subset025.R & NN_subset035.R). Unlike traditional ML algorithms, GBM constructs a series of decision trees sequentially, with each tree addressing the residuals of its predecessors, thereby progressively refining the model's accuracy. The GBM model stands out by building on an ensemble of decision trees, enhancing the model's overall predictive strength without relying on single-point predictions or specific node-level activation functions.

## Scoring Metric
The utilization of the F1 score remains central to patient categorization, especially in identifying groups with suboptimal survival outcomes. This metric is pivotal for its balanced assessment of precision and recall, addressing the challenges posed by imbalanced data classifications in binary contexts.

### Calculations
- **Precision** " = D / (D + (l/α)) "
- **Recall** " = D / (D + σ) "
- **F1 Score** " = 2 * (precision * recall) / (precision + recall) "

Where:
- D = Count of patients within the cluster who died due to cancer.
- l = Number of days counted in lc_from_IPDD.
- α = Median survival period before cancer deaths across all patient data.
- σ = Number of known patients that succumbed to cancer in another cluster.

## Discussion & Conclusion
The previous model phase yielded an F1 score of 0.611 during training, which increased to 0.800 in the testing phase for methylation data. Despite efforts in lvl4B to surpass these results, the testing data's atypical sampling precluded improvement. An innovative approach involved a blind resampling and splitting the data (lvl3B), comparing the neural network model's testing F1 score with the comprehensive F1 score from the prior phase. This comparison revealed an enhancement, with the neural network model achieving an F1 score of 0.666 on the blindly chosen 20% test data.

Analysis through the gradient-boosted model shed light on the challenges of surpassing the previous phase's F1 score. The feature importance analysis revealed that the CpG data's predictive power for the poor overall survival group was 53.6% in the lvl3B analysis, indicating a majority reliance on this data for phase predictions.
