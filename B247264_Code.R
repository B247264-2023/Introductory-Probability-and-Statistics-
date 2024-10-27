setwd("~/Desktop/Introductory Probability and Statistics")

# Load necessary library
library(dplyr)
library(readxl)

# Load covariates.xlsx and biomarkers.xlsx from the working dir
covariates <- read_excel("covariates.xlsx")
BIO <- read_excel("biomarkers.xlsx")

#### Task 1
## Merge the two dataframes based on the PatientID
# Filter the biomarkers dataframe to only include rows with "0weeks" in the Biomarker column
biomarkers_0weeks <- BIO[grep("0weeks", BIO$Biomarker), ]
# Extract the patient IDs from the Biomarker column (assuming IDs are the part before "-0weeks")
biomarkers_0weeks$patientID <- sub("-0weeks", "", biomarkers_0weeks$Biomarker)
# Merge the biomarkers_0weeks dataframe with the covariates dataframe on the correct patient ID column
biomarkers_merged <- merge(biomarkers_0weeks, covariates[, c("PatientID", "VAS-at-inclusion")], by.x = "patientID", by.y = "PatientID")

## Classify the patients into high and low pain groups by creating a new column
# This conditional statement assigns patients that have a VAS at inclusionf equal or higher than 5 to the high pain group
biomarkers_merged$Group <- ifelse(biomarkers_merged$`VAS-at-inclusion` >= 5, "High", "Low")
# Drop the Biomarker column from the final dataframe
biomarkers_merged <- biomarkers_merged[, !names(biomarkers_merged) %in% c("Biomarker", "VAS-at-inclusion")]

## Assess the normality of our data using histograms and the shapiro.test
# List of biomarker columns 
biomarkers <- c("IL-6", "VEGF-A", "OPG", "TGF-beta-1", "IL-8", 
                "CXCL9", "CXCL1", "IL-18", "CSF-1")

# Loop through each biomarker
for (biomarker in biomarkers) {
  # Perform the Shapiro-Wilk test for Group "Low"
  shapiro_low <- shapiro.test(biomarkers_merged[biomarkers_merged$Group == "Low", biomarker])
  
  # Perform the Shapiro-Wilk test for Group "High"
  shapiro_high <- shapiro.test(biomarkers_merged[biomarkers_merged$Group == "High", biomarker])
  
  # Calculate variance for Group "Low"
  variance_low <- var(biomarkers_merged[biomarkers_merged$Group == "Low", biomarker], na.rm = TRUE)
  
  # Calculate variance for Group "High"
  variance_high <- var(biomarkers_merged[biomarkers_merged$Group == "High", biomarker], na.rm = TRUE)
  
  # Print the results
  cat("Results for", biomarker, ":\n")
  cat("Group Low - W:", shapiro_low$statistic, ", p-value:", round(shapiro_low$p.value, 4), 
      ", Variance:", round(variance_low, 4), "\n")
  cat("Group High - W:", shapiro_high$statistic, ", p-value:", round(shapiro_high$p.value, 4), 
      ", Variance:", round(variance_high, 4), "\n\n")
  
  # Plot histograms for each group
  par(mfrow = c(1, 2))  # Set layout for side-by-side plots
  
  # Histogram for Group "Low"
  hist(biomarkers_merged[biomarkers_merged$Group == "Low", biomarker], 
       main = paste(biomarker, "- Group Low"), 
       xlab = biomarker, 
       col = "lightblue", 
       border = "black")
  
  # Histogram for Group "High"
  hist(biomarkers_merged[biomarkers_merged$Group == "High", biomarker], 
       main = paste(biomarker, "- Group High"), 
       xlab = biomarker, 
       col = "lightgreen", 
       border = "black")
}

## Perform hypothesis testing using t-test
# Loop through each biomarker
for (biomarker in biomarkers) {
  # Perform t-test between the two groups
  t_test_result <- t.test(biomarkers_merged[biomarkers_merged$Group == "Low", biomarker],
                          biomarkers_merged[biomarkers_merged$Group == "High", biomarker],
                          var.equal = FALSE)  # Allow for unequal variances
  # Print t-test results
  cat("Results for", biomarker, ":\n")
  cat("t-test - t:", round(t_test_result$statistic, 4), ", p-value:", round(t_test_result$p.value, 4), 
      ", Mean Group Low:", round(t_test_result$estimate[1], 4), 
      ", Mean Group High:", round(t_test_result$estimate[2], 4), "\n\n")
}

## Calculate the propability of making at least one type I error assuming that your 
## tests are independent and that all null hypotheses are true.
# Given values
alpha <- 0.05
m <- 9
# Probability of at least one Type I error
prob_at_least_one_error <- 1 - (1 - alpha)^m
# Print the result
cat("The probability of making at least one Type I error is:", round(prob_at_least_one_error, 4), "\n")

## Perform hypothesis testing using t-test with Bonferroni correction
# Number of biomarkers
num_biomarkers <- length(biomarkers)

# Loop through each biomarker
for (biomarker in biomarkers) {
  # Perform t-test between the two groups
  t_test_result <- t.test(biomarkers_merged[biomarkers_merged$Group == "Low", biomarker],
                          biomarkers_merged[biomarkers_merged$Group == "High", biomarker],
                          var.equal = FALSE)  # Allow for unequal variances
  
  # Calculate the Bonferroni corrected p-value
  bonferroni_p_value <- p.adjust(t_test_result$p.value, method = "bonferroni", n = num_biomarkers)
  
  # Print t-test results
  cat("Results for", biomarker, ":\n")
  cat("t-test - t:", round(t_test_result$statistic, 4), 
      ", p-value (Bonferroni corrected):", round(bonferroni_p_value, 4), 
      ", Mean Group Low:", round(t_test_result$estimate[1], 4), 
      ", Mean Group High:", round(t_test_result$estimate[2], 4), "\n\n")
}

#### Task 2
# Merge the biomarkers_0weeks dataframe with the covariates dataframe on the correct patient ID column
biomarkers_merged2 <- merge(biomarkers_0weeks, covariates, by.x = "patientID", by.y = "PatientID", all.x = TRUE)
# Remove columns that are not needed
biomarkers_merged2$Biomarker <- NULL

# Rename columns in biomarkers_merged2 to avoid issues with special characters
colnames(biomarkers_merged2) <- gsub("-", "_", colnames(biomarkers_merged2))  # Replace dashes with underscores
names(biomarkers_merged2)[names(biomarkers_merged2) == "Smoker (1=yes, 2=no)"] <- "Smoker"
names(biomarkers_merged2)[names(biomarkers_merged2) == "Sex (1=male, 2=female)"] <- "Sex"

# Split the data into 80% and 20% for testing by randomly dividing the biomarkers_merged2 dataframe 
set.seed(123) # For reproducibility
sample_index <- sample(seq_len(nrow(biomarkers_merged2)), size = 0.8 * nrow(biomarkers_merged2))
train_data <- biomarkers_merged2[sample_index, ]
test_data <- biomarkers_merged2[-sample_index, ]

# Since Sex and Smoker are categorical variables, we need to tell R that so it does not handle them as
# continuous which are the rest.
train_data$Sex <- factor(train_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
train_data$Smoker <- factor(train_data$Smoker, levels = c(1, 2), labels = c("Yes", "No"))

# List of explanatory variables
explanatory_variables <- c("IL_6", "VEGF_A", "OPG", "TGF_beta_1", "IL_8", 
                           "CXCL9", "CXCL1", "IL_18", "CSF_1", "Age", 
                           "VAS_at_inclusion")

# Create scatterplots for each explanatory variable
par(mfrow = c(2, 2)) 
for (var in explanatory_variables) {
  plot(train_data[[var]], train_data$Vas_12months, 
       main = paste("Scatterplot of Vas_12months vs", var),
       xlab = var, ylab = "Vas_12months")
  
  # Fit a linear regression model for the current variable
  model <- lm(Vas_12months ~ train_data[[var]], data = train_data)
  
  # Add the fitted line to the plot
  abline(model, col = "red")  # Fitted line in red
}

# Fit the linear model using all the explanatory variables of the train_data
model <- lm(Vas_12months ~ IL_6 + VEGF_A + OPG + TGF_beta_1 + IL_8 + CXCL9 + 
              CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + VAS_at_inclusion, 
            data = train_data)

summary(model)

par(mfrow = c(1, 2)) 

# Create a residuals plot
plot(predict(model), residuals(model),
     xlab = "Predicted Values",
     ylab = "Residuals",
     main = "Residuals vs Predicted Values",
     pch = 19, col = "blue")

# Create a histogram of the residuals
hist(residuals(model),
     breaks = 20,  # Number of bins
     main = "Histogram of Residuals",
     xlab = "Residuals",
     col = "lightblue",
     border = "black")

## The untrasformed model appear to not have randomly scattered residuals 

# Square root transformation
train_data$Vas_12months_sqrt <- sqrt(train_data$Vas_12months)

# Fit the linear model with the transformed response variable
model_sqrt <- lm(Vas_12months_sqrt ~ IL_6 + VEGF_A + OPG + TGF_beta_1 + IL_8 + 
                   CXCL9 + CXCL1 + IL_18 + CSF_1 + Age + Sex + Smoker + 
                   VAS_at_inclusion, data = train_data)

# Summary of the model
summary(model_sqrt)

# Create a residuals plot
plot(predict(model_sqrt), residuals(model_sqrt),
     xlab = "Predicted Values",
     ylab = "Residuals",
     main = "Residuals vs Predicted Values",
     pch = 19, col = "blue")

# Create a histogram of the residuals
hist(residuals(model_sqrt),
     breaks = 20,  # Number of bins
     main = "Histogram of Residuals",
     xlab = "Residuals",
     col = "lightblue",
     border = "black")


# Create a residuals plot for the untransformed response variable
plot(predict(model), residuals(model),
     xlab = "Predicted Values",
     ylab = "Residuals",
     main = "Residual Plot - Original scale",
     pch = 19, col = "blue")

# Create a residuals plot for the sqrt-transformed response variable
plot(predict(model_sqrt), residuals(model_sqrt),
     xlab = "Predicted Values",
     ylab = "Residuals",
     main = "Residual Plot - Square root scale",
     pch = 19, col = "blue")

## Predict with the test_data
# Since Sex and Smoker are categorical variables, we need to tell R that so it does not handle them as
# continuous which are the rest.
test_data$Sex <- factor(test_data$Sex, levels = c(1, 2), labels = c("Male", "Female"))
test_data$Smoker <- factor(test_data$Smoker, levels = c(1, 2), labels = c("Yes", "No"))

# Make predictions using the sqrt model
test_data$Vas_12months_sqrt <- sqrt(test_data$Vas_12months)
predictions_sqrt <- predict(model_sqrt, newdata = test_data)

# Calculate R-squared using the actual values
actual_values <- test_data$Vas_12months_sqrt
r_squared <- summary(lm(predictions_sqrt ~ actual_values))$r.squared

par(mfrow = c(1, 1))
# Create a scatter plot of the actual vs Predicted values
plot(actual_values, predictions_sqrt,
     main = "Predicted vs Actual 12-month VAS Scores",
     xlab = "Actual 12-month VAS Scores",
     ylab = "Predicted 12-month VAS Scores",
     pch = 19, col = "blue")                 
# Add a 45-degree reference line (y = x)
abline(a = 0, b = 1, col = "red", lwd = 2) # Add line y = x for reference
