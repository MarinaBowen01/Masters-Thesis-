
rm(list = ls())

#### usnea ####
### USNEA TRANSFORMATIONS FOR HURDLE MODEL 
library(randomForest)
library(quantregForest)
library(data.table)
library(dplyr)
library(rsample)
library(isotone)
library(pROC)
library(pdp)
library(clipr)
library(ggplot2)


select <- dplyr::select

# Load data
df <- read.csv("Usnea.csv")
df <- df %>% mutate(across(where(is.character), as.factor))
df$Presence <- as.factor(ifelse(df$Usnea > 0, 1, 0))

# Define transformation functions
transformations <- list(
  none = list(
    transform = function(x) x,
    inverse = function(x) x,
    name = "None"
  ),
  sqrt = list(
    transform = function(x) sqrt(pmax(x, 0)),
    inverse = function(x) pmax(x^2, 0),
    name = "Square Root"
  ),
  log = list(
    transform = function(x) log(pmax(x, 0.001)),  # Add small constant to avoid log(0)
    inverse = function(x) exp(x),
    name = "Log"
  ),
  asinh = list(
    transform = function(x) asinh(x),
    inverse = function(x) sinh(x),
    name = "Asinh"
  )
)

predictors <- c("PDIS_CLASS", "ASP", "SLOPE_CLASS", "MORF_CLASS", "CD_CLASS", "ALTITUDE", "DISTCOAST",
                "WEXP", "HURS_09", "WEFF", "VDEP", "TRI", "GDGFGD0", "MXCURV", "TTCURV", "TNCURV", "Bio18_81",
                "Easting", "Northing")
factor_vars <- names(df[, predictors])[sapply(df[, predictors], is.factor)]

n_iter <- 100

# Initialize storage for all transformations
all_results_clf <- list()
all_results_qrf <- list()
varimp_list_clf <- list()
varimp_list_qrf <- list()
predictions_all <- list()
training_data_list <- list()
test_data_list <- list()

for (trans_name in names(transformations)) {
  trans_func <- transformations[[trans_name]]
  
  cat("Testing transformation:", trans_func$name, "\n")
  
  # Transform the response variable
  df_trans <- df
  df_trans$Usnea_trans <- trans_func$transform(df$Usnea)
  
  all_results_clf[[trans_name]] <- list()
  all_results_qrf[[trans_name]] <- list()
  varimp_list_clf[[trans_name]] <- list()
  varimp_list_qrf[[trans_name]] <- list()
  predictions_all[[trans_name]] <- list()
  training_data_list[[trans_name]] <- list()
  test_data_list[[trans_name]] <- list()
  
  for (i in 1:n_iter) {
    set.seed(i)
    repeat {
      split <- initial_split(df_trans, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      all_ok <- TRUE
      for (fac in factor_vars) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) {
          all_ok <- FALSE
          break
        }
      }
      if (all_ok) break
    }
    
    for (fac in factor_vars) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
    
    ## STEP 1: RF Classifier (same for all transformations)
    rf_clf <- randomForest(
      x = train[, predictors],
      y = train$Presence,
      ntree = 500,
      importance = TRUE,
      keep.forest = TRUE
    )
    
    pred_probs <- predict(rf_clf, newdata = test[, predictors], type = "prob")[, "1"]
    pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
    obs_class <- test$Presence
    
    # Metrics for RF classifier
    acc <- mean(pred_class == obs_class)
    sens <- mean(pred_class == "1" & obs_class == "1")
    spec <- mean(pred_class == "0" & obs_class == "0")
    logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
    brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
    auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
    
    spread_clf <- abs(0.5 - pred_probs)  # proxy for uncertainty
    
    all_results_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
      LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
      MedianSpread = median(spread_clf)
    )
    
    imp_clf <- importance(rf_clf)
    varimp_list_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Variable = rownames(imp_clf), 
      Importance = imp_clf[, 1], 
      Iteration = i, 
      Model = "RF_Presence"
    )
    
    ## STEP 2: QRF for predicted present cases (with transformation)
    test_presence_idx <- which(pred_class == "1")
    test_present_data <- test[test_presence_idx, ]
    train_present_data <- train[train$Presence == "1", ]
    
    if (nrow(test_present_data) > 10) {
      qrf <- quantregForest(
        x = train_present_data[, predictors],
        y = train_present_data$Usnea_trans,  # Use transformed response
        ntree = 500
      )
      
      preds_test_trans <- predict(qrf, newdata = test_present_data[, predictors], what = c(0.025, 0.5, 0.975))
      
      # Apply isotonic regression on transformed scale
      iso_model <- isoreg(preds_test_trans[, 2], test_present_data$Usnea_trans)
      iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
      calibrated_test_trans <- iso_fun(preds_test_trans[, 2])
      
      # Back-transform predictions and bounds
      preds_test_lower <- trans_func$inverse(preds_test_trans[, 1])
      preds_test_median <- trans_func$inverse(calibrated_test_trans)
      preds_test_upper <- trans_func$inverse(preds_test_trans[, 3])
      
      # Ensure predictions are within reasonable bounds
      preds_test_lower <- pmax(pmin(preds_test_lower, 100), 0)
      preds_test_median <- pmax(pmin(preds_test_median, 100), 0)
      preds_test_upper <- pmax(pmin(preds_test_upper, 100), 0)
      
      spread_test <- preds_test_upper - preds_test_lower
      resids_test <- preds_test_median - test_present_data$Usnea
      
      rmse <- sqrt(mean(resids_test^2))
      mae <- mean(abs(resids_test))
      bias <- mean(resids_test)
      r2 <- cor(test_present_data$Usnea, preds_test_median)^2
      nmae <- mae / mean(test_present_data$Usnea)
      coverage <- mean(test_present_data$Usnea >= preds_test_lower & test_present_data$Usnea <= preds_test_upper)
      median_spread <- median(spread_test)
      low_bias <- mean(preds_test_lower - test_present_data$Usnea)
      high_bias <- mean(preds_test_upper - test_present_data$Usnea)
      
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
        R_squared = r2, nMAE = nmae, Coverage_95PI = coverage,
        Median_Spread_95PI = median_spread,
        Lower_Q2.5_Bias = low_bias,
        Upper_Q97.5_Bias = high_bias
      )
      
      # Store predictions (back-transformed for training data too)
      train_preds_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.5)
      train_spread_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.975) -
        predict(qrf, newdata = train_present_data[, predictors], what = 0.025)
      
      predictions_all[[trans_name]][[i]] <- rbind(
        data.table(set = "test", obs = test_present_data$Usnea, pred = preds_test_median,
                   spread = spread_test, iteration = i, transformation = trans_func$name),
        data.table(set = "train", obs = train_present_data$Usnea,
                   pred = trans_func$inverse(train_preds_trans),
                   spread = trans_func$inverse(train_present_data$Usnea_trans + train_spread_trans/2) - 
                     trans_func$inverse(train_present_data$Usnea_trans - train_spread_trans/2),
                   iteration = i, transformation = trans_func$name)
      )
      
      imp_qrf <- importance(qrf)
      varimp_list_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Variable = rownames(imp_qrf), 
        Importance = imp_qrf[, 1], 
        Iteration = i, 
        Model = "QRF_Abundance"
      )
    } else {
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
        R_squared = NA, nMAE = NA, Coverage_95PI = NA,
        Median_Spread_95PI = NA,
        Lower_Q2.5_Bias = NA,
        Upper_Q97.5_Bias = NA
      )
      predictions_all[[trans_name]][[i]] <- NULL
      varimp_list_qrf[[trans_name]][[i]] <- NULL
    }
    
    training_data_list[[trans_name]][[i]] <- train
    test_data_list[[trans_name]][[i]] <- test
    
    if (i %% 10 == 0) cat("Transformation:", trans_func$name, "- Iteration", i, "\n")
  }
}

# Combine results across transformations
results_clf_combined <- rbindlist(lapply(all_results_clf, function(x) rbindlist(x)))
results_qrf_combined <- rbindlist(lapply(all_results_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_clf_combined <- rbindlist(lapply(varimp_list_clf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_qrf_combined <- rbindlist(lapply(varimp_list_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
predictions_combined <- rbindlist(lapply(predictions_all, function(x) rbindlist(x[!sapply(x, is.null)])))

# Summary statistics by transformation
cat("\n=== CLASSIFICATION RESULTS SUMMARY ===\n")
clf_summary <- results_clf_combined[, .(
  Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
  SD_Accuracy = sd(Accuracy, na.rm = TRUE),
  Mean_AUC = mean(AUC, na.rm = TRUE),
  SD_AUC = sd(AUC, na.rm = TRUE),
  Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
  Mean_Specificity = mean(Specificity, na.rm = TRUE)
), by = Transformation]
print(clf_summary)

cat("\n=== QUANTILE REGRESSION RESULTS SUMMARY ===\n")
qrf_summary <- results_qrf_combined[, .(
  Mean_RMSE = mean(RMSE, na.rm = TRUE),
  SD_RMSE = sd(RMSE, na.rm = TRUE),
  Mean_MAE = mean(MAE, na.rm = TRUE),
  SD_MAE = sd(MAE, na.rm = TRUE),
  Mean_R2 = mean(R_squared, na.rm = TRUE),
  SD_R2 = sd(R_squared, na.rm = TRUE),
  Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
  Mean_Bias = mean(Bias, na.rm = TRUE)
), by = Transformation]
print(qrf_summary)

# 1. PREPARE DATA FOR PLOTTING

# Spread plot data (prediction interval width vs observed)
spread_plot_data <- predictions_combined[set == "test", .(
  obs = obs,
  spread = spread,
  trans = transformation,
  iteration = iteration
)]

# Observed vs predicted data
pred_obs_data <- predictions_combined[set == "test", .(
  obs = obs,
  pred = pred,
  trans = transformation,
  iteration = iteration
)]

# Calculate residuals for each transformation
residuals_data <- predictions_combined[set == "test", .(
  pred = pred,
  residuals = obs - pred,
  trans = transformation,
  iteration = iteration
)]

# 2. SPREAD PLOT (Uncertainty vs Observed)
p1 <- ggplot(spread_plot_data, aes(x = obs, y = spread, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Observed Usnea Abundance", y = "Prediction Interval Width",
       title = "Uncertainty vs Observed (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# 3. VARIABLE IMPORTANCE ANALYSIS
cat("\n=== VARIABLE IMPORTANCE ANALYSIS ===\n")

# Calculate variable importance summary
varimp_summary <- varimp_qrf_combined %>%
  group_by(Transformation, Variable) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 4),
            .groups = "drop") %>%
  arrange(Transformation, desc(MedianImportance))

# Calculate percentage importance
varimp_summary <- varimp_summary %>%
  group_by(Transformation) %>%
  mutate(Importance_pct = 100 * MedianImportance / sum(MedianImportance)) %>%
  ungroup()

# For "None" transformation (equivalent to "raw")
none_imp <- varimp_summary %>%
  filter(Transformation == "None") %>%
  mutate(Cumulative_Importance = round(cumsum(Importance_pct), 1)) %>%
  arrange(desc(MedianImportance))

cat("Variable Importance for No Transformation (Raw):\n")
print(none_imp)
clipr::write_clip(none_imp)

# Variable importance comparison plot
p2 <- varimp_summary %>%
  group_by(Transformation) %>%
  slice_head(n = 20) %>%  # Top 10 variables per transformation
  ggplot(aes(x = reorder(Variable, MedianImportance), y = MedianImportance, fill = Transformation)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Variable", y = "Median Importance", 
       title = "Top 10 Variable Importance by Transformation") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  facet_wrap(~Transformation, scales = "free_y")

print(p2)

# 4. OBSERVED VS PREDICTED PLOT
p3 <- ggplot(pred_obs_data, aes(x = obs, y = pred, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed Usnea Abundance", y = "Predicted Usnea Abundance",
       title = "Observed vs Predicted (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# 5. RESIDUALS VS PREDICTED (for each transformation)
p4 <- ggplot(residuals_data, aes(x = pred, y = residuals)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~trans, scales = "free") +
  theme_minimal() +
  labs(title = "Residuals vs Predicted (by Transformation)",
       x = "Predicted Usnea Abundance",
       y = "Residuals (Observed - Predicted)")

print(p4)

# 6. INDIVIDUAL TRANSFORMATION PLOTS

# Function to create plots for a specific transformation
create_transformation_plots <- function(trans_name) {
  
  # Filter data for specific transformation
  trans_spread <- spread_plot_data[trans == trans_name]
  trans_pred_obs <- pred_obs_data[trans == trans_name]
  trans_residuals <- residuals_data[trans == trans_name]
  
  cat(paste("\n=== PLOTS FOR", trans_name, "TRANSFORMATION ===\n"))
  
  # Residuals vs Predicted for specific transformation
  p_resid <- ggplot(trans_residuals, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Residuals vs Predicted -", trans_name, "Transformation"),
         x = "Predicted Usnea Abundance",
         y = "Residuals (Observed - Predicted)")
  
  # Prediction Interval Width vs Observed for specific transformation
  p_spread <- ggplot(trans_spread, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    theme_minimal() +
    labs(title = paste("Prediction Interval Width vs Observed -", trans_name, "Transformation"),
         x = "Observed Usnea Abundance",
         y = "95% Prediction Interval Width")
  
  print(p_resid)
  print(p_spread)
  
  return(list(residuals = p_resid, spread = p_spread))
}

# Create individual plots for each transformation
transformation_plots <- list()
for (trans in unique(spread_plot_data$trans)) {
  transformation_plots[[trans]] <- create_transformation_plots(trans)
}


# 7. PERFORMANCE SUMMARY COMPARISON PLOT
perf_summary_long <- results_qrf_combined %>%
  as.data.frame() %>%
  dplyr::select(Transformation, RMSE, MAE, R_squared, Coverage_95PI, Bias) %>%
  group_by(Transformation) %>%
  summarize(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Mean_MAE = mean(MAE, na.rm = TRUE),
    Mean_R2 = mean(R_squared, na.rm = TRUE),
    Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
    Mean_Bias = abs(mean(Bias, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p5 <- ggplot(perf_summary_long, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_col() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Model Performance Comparison Across Transformations",
       x = "Transformation", y = "Mean Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p5)

# 8. BOXPLOT OF PERFORMANCE METRICS
perf_boxplot_data <- results_qrf_combined %>%
  as.data.frame() %>%
  dplyr::select(Transformation, RMSE, MAE, R_squared, Coverage_95PI) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p6 <- ggplot(perf_boxplot_data, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Distribution of Performance Metrics Across Iterations",
       x = "Transformation", y = "Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p6)

# 9. CORRELATION BETWEEN OBSERVED AND PREDICTED BY TRANSFORMATION
cor_summary <- pred_obs_data %>%
  group_by(trans) %>%
  summarize(
    correlation = cor(obs, pred, use = "complete.obs"),
    rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    mae = mean(abs(obs - pred), na.rm = TRUE),
    .groups = "drop")

print(cor_summary)

#### EXTRACT TOP PREDICTORS FROM TRANSFORMATION ANALYSIS ####

# Extract top predictors from presence model (RF Classifier) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_presence <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

# Extract top predictors from abundance model (QRF) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_abundance <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

cat("\n=== TOP PREDICTORS EXTRACTED FROM RAW TRANSFORMATION (70% IMPORTANCE) ===\n")
cat("Top predictors for PRESENCE model (explaining 70% of importance):\n")
print(top_predictors_presence)
cat("Number of presence predictors:", length(top_predictors_presence), "\n")

cat("\nTop predictors for ABUNDANCE model (explaining 70% of importance):\n")
print(top_predictors_abundance)
cat("Number of abundance predictors:", length(top_predictors_abundance), "\n")

# Optional: Show the importance breakdown
cat("\n=== IMPORTANCE BREAKDOWN FOR SELECTED PREDICTORS ===\n")

# Presence model importance breakdown
presence_importance_breakdown <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_presence
  ) %>%
  filter(Selected)

cat("PRESENCE model - Selected variables and their importance:\n")
print(presence_importance_breakdown)

# Abundance model importance breakdown
abundance_importance_breakdown <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_abundance
  ) %>%
  filter(Selected)

cat("\nABUNDANCE model - Selected variables and their importance:\n")
print(abundance_importance_breakdown)





################### USNEA HURDLE FINAL MODEL USING TOP PREDICTORS #########################################################

# Use the extracted top predictors for each model
predictors_presence <- c("ALTITUDE", "WEXP", "Easting", "Northing", "HURS_09", "VDEP", "DISTCOAST", "WEFF", "ASP", "Bio18_81")
predictors_abundance <- c("Northing", "ALTITUDE", "WEXP", "ASP", "DISTCOAST", "Easting", "WEFF", "MXCURV", "VDEP")

# Get factor variables for each predictor set
factor_vars_presence <- names(df[, predictors_presence])[sapply(df[, predictors_presence], is.factor)]
factor_vars_abundance <- names(df[, predictors_abundance])[sapply(df[, predictors_abundance], is.factor)]

n_iter <- 100
all_results_clf_final <- list()
all_results_qrf_final <- list()
varimp_list_clf_final <- list()
varimp_list_qrf_final <- list()
predictions_all_final <- list()
training_data_list_final <- list()
test_data_list_final <- list()

cat("\n=== RUNNING FINAL MODEL WITH TOP PREDICTORS ===\n")

for (i in 1:n_iter) {
  set.seed(i)
  repeat {
    split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
    train <- training(split)
    test <- testing(split)
    
    # Check factor levels for BOTH predictor sets
    all_ok <- TRUE
    for (fac in factor_vars_presence) {
      train_levels <- levels(train[[fac]])
      test_levels <- levels(factor(test[[fac]]))
      if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
    }
    if (all_ok) {
      for (fac in factor_vars_abundance) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
      }
    }
    if (all_ok) break
  }
  
  # Set factor levels
  for (fac in factor_vars_presence) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  for (fac in factor_vars_abundance) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  
  ## STEP 1: RF Classifier
  rf_clf <- randomForest(
    x = train[, predictors_presence],
    y = train$Presence,
    ntree = 500,
    importance = TRUE,
    keep.forest = TRUE
  )
  
  pred_probs <- predict(rf_clf, newdata = test[, predictors_presence], type = "prob")[, "1"]
  pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
  obs_class <- test$Presence
  
  # Metrics for classifier
  acc <- mean(pred_class == obs_class)
  sens <- mean(pred_class == "1" & obs_class == "1")
  spec <- mean(pred_class == "0" & obs_class == "0")
  logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
  brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
  auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
  spread_clf <- abs(0.5 - pred_probs)
  
  all_results_clf_final[[i]] <- data.table(
    Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
    LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
    MedianSpread = median(spread_clf)
  )
  
  varimp_list_clf_final[[i]] <- data.table(
    Variable = rownames(importance(rf_clf)),
    Importance = importance(rf_clf)[,1],
    Iteration = i, Model = "RF_Presence"
  )
  
  ## STEP 2: QRF for abundance
  test_presence_idx <- which(pred_class == "1")
  test_present_data <- test[test_presence_idx, ]
  train_present_data <- train[train$Presence == "1", ]
  
  if (nrow(test_present_data) > 10) {
    qrf <- quantregForest(
      x = train_present_data[, predictors_abundance],
      y = train_present_data$Usnea,
      ntree = 500
    )
    
    preds_test <- predict(qrf, newdata = test_present_data[, predictors_abundance], what = c(0.025, 0.5, 0.975))
    iso_model <- isoreg(preds_test[, 2], test_present_data$Usnea)
    iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
    calibrated_test <- pmax(pmin(iso_fun(preds_test[, 2]), 100), 0)
    
    spread_test <- preds_test[, 3] - preds_test[, 1]
    resids_test <- calibrated_test - test_present_data$Usnea
    
    rmse <- sqrt(mean(resids_test^2))
    mae <- mean(abs(resids_test))
    bias <- mean(resids_test)
    r2 <- cor(test_present_data$Usnea, calibrated_test)^2
    nrmse <- rmse / mean(test_present_data$Usnea, na.rm = TRUE)
    nrmse_median <- rmse / median(test_present_data$Usnea, na.rm = TRUE)
    nmae <- mae / mean(test_present_data$Usnea, na.rm = TRUE)
    nmae_median <- mae / median(test_present_data$Usnea, na.rm = TRUE)
    coverage <- mean(test_present_data$Usnea >= preds_test[,1] & test_present_data$Usnea <= preds_test[,3])
    median_spread <- median(spread_test)
    low_bias <- mean(preds_test[,1] - test_present_data$Usnea)
    high_bias <- mean(preds_test[,3] - test_present_data$Usnea)
    
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
      R_squared = r2,
      nRMSE = nrmse,
      nRMSE_Median = nrmse_median,
      nMAE = nmae,
      nMAE_Median = nmae_median,
      Coverage_95PI = coverage,
      Median_Spread_95PI = median_spread,
      Lower_Q2.5_Bias = low_bias,
      Upper_Q97.5_Bias = high_bias
    )
    
    predictions_all_final[[i]] <- rbind(
      data.table(set = "test", obs = test_present_data$Usnea, pred = calibrated_test,
                 spread = spread_test, iteration = i),
      data.table(set = "train", obs = train_present_data$Usnea,
                 pred = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.5),
                 spread = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.975) -
                   predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.025),
                 iteration = i)
    )
    
    varimp_list_qrf_final[[i]] <- data.table(
      Variable = rownames(importance(qrf)),
      Importance = importance(qrf)[,1],
      Iteration = i, Model = "QRF_Abundance"
    )
  } else {
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
      R_squared = NA, nRMSE = NA, nRMSE_Median = NA, nMAE = NA, nMAE_Median = NA,
      Coverage_95PI = NA, Median_Spread_95PI = NA, Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA
    )
    predictions_all_final[[i]] <- NULL
    varimp_list_qrf_final[[i]] <- NULL
  }
  
  training_data_list_final[[i]] <- train
  test_data_list_final[[i]] <- test
  
  if (i %% 10 == 0) cat("Final model iteration", i, "\n")
}

# Combine results
results_clf_final_combined <- rbindlist(all_results_clf_final)
results_qrf_final_combined <- rbindlist(all_results_qrf_final[!sapply(all_results_qrf_final, is.null)])
varimp_clf_final_combined <- rbindlist(varimp_list_clf_final[!sapply(varimp_list_clf_final, is.null)])
varimp_qrf_final_combined <- rbindlist(varimp_list_qrf_final[!sapply(varimp_list_qrf_final, is.null)])
predictions_final_combined <- rbindlist(predictions_all_final[!sapply(predictions_all_final, is.null)])

# CLASSIFICATION SUMMARY
clf_final_summary <- results_clf_final_combined[, .(
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 3),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 3),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 3),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 3),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 3),
  SD_Sensitivity = round(sd(Sensitivity, na.rm = TRUE), 3),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 3),
  SD_Specificity = round(sd(Specificity, na.rm = TRUE), 3),
  Mean_LogLoss = round(mean(LogLoss, na.rm = TRUE), 3),
  SD_LogLoss = round(sd(LogLoss, na.rm = TRUE), 3),
  Mean_BrierScore = round(mean(BrierScore, na.rm = TRUE), 3),
  SD_BrierScore = round(sd(BrierScore, na.rm = TRUE), 3),
  Mean_MedianSpread = round(mean(MedianSpread, na.rm = TRUE), 3),
  SD_MedianSpread = round(sd(MedianSpread, na.rm = TRUE), 3)
)]
print(clf_final_summary)

qrf_final_summary <- results_qrf_final_combined[, .(
  Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
  SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
  Median_RMSE = round(median(RMSE, na.rm = TRUE), 3),
  Mean_MAE = round(mean(MAE, na.rm = TRUE), 3),
  SD_MAE = round(sd(MAE, na.rm = TRUE), 3),
  Median_MAE = round(median(MAE, na.rm = TRUE), 3),
  Mean_R2 = round(mean(R_squared, na.rm = TRUE), 3),
  SD_R2 = round(sd(R_squared, na.rm = TRUE), 3),
  Median_R2 = round(median(R_squared, na.rm = TRUE), 3),
  Median_nRMSE = round(median(nRMSE_Median, na.rm = TRUE), 3),
  IQR_nRMSE   = round(IQR(nRMSE_Median, na.rm = TRUE), 3),
  Mean_Bias = round(mean(Bias, na.rm = TRUE), 3),
  SD_Bias = round(sd(Bias, na.rm = TRUE), 3),
  Median_Bias = round(median(Bias, na.rm = TRUE), 3),
  Mean_Coverage_95PI = round(mean(Coverage_95PI, na.rm = TRUE), 3),
  SD_Coverage_95PI = round(sd(Coverage_95PI, na.rm = TRUE), 3),
  Median_Coverage_95PI = round(median(Coverage_95PI, na.rm = TRUE), 3),
  Mean_Median_Spread_95PI = round(mean(Median_Spread_95PI, na.rm = TRUE), 3),
  SD_Median_Spread_95PI = round(sd(Median_Spread_95PI, na.rm = TRUE), 3),
  Median_Spread_95PI = round(median(Median_Spread_95PI, na.rm = TRUE), 3),
  Mean_Lower_Q2.5_Bias = round(mean(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  SD_Lower_Q2.5_Bias = round(sd(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Median_Lower_Q2.5_Bias = round(median(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Mean_Upper_Q97.5_Bias = round(mean(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  SD_Upper_Q97.5_Bias = round(sd(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  Median_Upper_Q97.5_Bias = round(median(Upper_Q97.5_Bias, na.rm = TRUE), 3)
)]
print(qrf_final_summary)





#### VARIABLE IMPORTANCE ANALYSIS FOR FINAL MODELS ####

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== VARIABLE IMPORTANCE FROM FINAL MODELS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create comprehensive importance tables from the final model runs
cat("\nTABLE 1: PRESENCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

presence_final_importance <- varimp_clf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  dplyr::select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(presence_final_importance)

cat("\nPRESENCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(presence_final_importance), "\n"))
cat(paste("- Top variable (", presence_final_importance$Variable[1], ") explains:", 
          presence_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", presence_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 2: ABUNDANCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

abundance_final_importance <- varimp_qrf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  dplyr::select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(abundance_final_importance)

cat("\nABUNDANCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(abundance_final_importance), "\n"))
cat(paste("- Top variable (", abundance_final_importance$Variable[1], ") explains:", 
          abundance_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", abundance_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

# Create a comparison table showing which variables are most important in each model
cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 3: COMPARISON - Top Variables in Each Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

# Get top 5 from each model for comparison
top_presence_vars <- presence_final_importance$Variable[1:min(5, nrow(presence_final_importance))]
top_abundance_vars <- abundance_final_importance$Variable[1:min(5, nrow(abundance_final_importance))]

# Create comparison dataframe
comparison_table <- data.frame(
  Rank = 1:5,
  Presence_Model = c(top_presence_vars, rep("", 5 - length(top_presence_vars)))[1:5],
  Presence_Importance = c(presence_final_importance$Importance_pct[1:min(5, nrow(presence_final_importance))], 
                          rep(NA, 5 - min(5, nrow(presence_final_importance))))[1:5],
  Abundance_Model = c(top_abundance_vars, rep("", 5 - length(top_abundance_vars)))[1:5],
  Abundance_Importance = c(abundance_final_importance$Importance_pct[1:min(5, nrow(abundance_final_importance))], 
                           rep(NA, 5 - min(5, nrow(abundance_final_importance))))[1:5]
)


print(comparison_table)

# Identify shared important variables
shared_vars <- intersect(top_presence_vars, top_abundance_vars)
if(length(shared_vars) > 0) {
  cat("\nShared important variables between models:", paste(shared_vars, collapse = ", "), "\n")
} else {
  cat("\nNo shared variables in top 5 between models\n")
}

#### PLOTTING FUNCTIONS - RUN THIS FIRST ####

# Define color palette for species
cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00")

library(ggplot2)
library(dplyr)

# Function to extract 70% importance variables from FINAL MODEL results
extract_final_70_percent <- function(varimp_final_combined, species_name, model_type) {
  
  # Calculate importance from final model runs and apply 70% threshold
  importance_data <- varimp_final_combined %>%
    group_by(Variable) %>%
    summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MedianImportance)) %>%
    mutate(
      Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
      Cumulative_pct = round(cumsum(Importance_pct), 2),
      Selected_70pct = Cumulative_pct <= 70,
      Rank = row_number(),
      Species = species_name,
      Model = model_type
    ) %>%
    filter(Selected_70pct) %>%  # Only keep variables that explain ≤70% cumulative importance
    select(Species, Model, Variable, Rank, Importance_pct, Cumulative_pct)
  
  return(importance_data)
}

create_individual_plot <- function(species_name, model_type, data_all) {
  
  data_subset <- data_all %>%
    filter(Species == species_name, Model == model_type) %>%
    arrange(desc(Rank))  # Most important at top
  
  if(nrow(data_subset) == 0) {
    cat(paste("No data found for", species_name, model_type, "\n"))
    return(NULL)
  }
  
  p <- ggplot(data_subset, aes(x = reorder(Variable, Rank), y = Importance_pct)) +
    geom_col(fill = cores[species_name], alpha = 0.8, color = "white", size = 1) +
    geom_text(aes(label = paste0(round(Importance_pct, 1), "%")), 
              hjust = -0.1, size = 5, color = "black") +
    coord_flip() +
    labs(
      title = paste(species_name, "-", model_type, "Final Model (70% Variables)"),
      subtitle = paste("Variables explaining", round(max(data_subset$Cumulative_pct), 1), "% importance in FINAL model"),
      x = "Variables (Ranked by Final Model Importance)",
      y = "Individual Importance (%)",
      caption = paste("n =", nrow(data_subset), "variables from final model")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", color = cores[species_name]),
      plot.subtitle = element_text(size = 13),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

create_final_model_plots <- function(species_name, varimp_clf_final, varimp_qrf_final) {
  
  cat(paste("\n=== Creating plots for", species_name, "final model ===\n"))
  
  # Extract 70% variables from final models
  presence_data <- extract_final_70_percent(varimp_clf_final, species_name, "Presence")
  abundance_data <- extract_final_70_percent(varimp_qrf_final, species_name, "Abundance")
  
  # Combine data
  combined_data <- rbind(presence_data, abundance_data)
  
  # Create individual plots
  presence_plot <- create_individual_plot(species_name, "Presence", combined_data)
  abundance_plot <- create_individual_plot(species_name, "Abundance", combined_data)
  
  # Print plots
  if(!is.null(presence_plot)) {
    print(presence_plot)
    cat(paste("Created", species_name, "presence plot with", nrow(presence_data), "variables\n"))
  }
  
  if(!is.null(abundance_plot)) {
    print(abundance_plot) 
    cat(paste("Created", species_name, "abundance plot with", nrow(abundance_data), "variables\n"))
  }
  
  # Create comparison plot
  if(nrow(combined_data) > 0) {
    comparison_plot <- ggplot(combined_data, aes(x = reorder(Variable, Rank), y = Importance_pct, fill = Model)) +
      geom_col(position = "dodge", alpha = 0.8, color = "white", size = 1) +
      coord_flip() +
      labs(
        title = paste(species_name, "- Final Model: Presence vs Abundance (70% Variables)"),
        subtitle = "Variable importance from final models (70% threshold applied)",
        x = "Variables",
        y = "Individual Importance (%)",
        fill = "Model Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
        plot.subtitle = element_text(size = 13),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      ) +
      scale_fill_manual(values = c("Presence" = cores[species_name], "Abundance" = adjustcolor(cores[species_name], alpha.f = 0.6))) +
      facet_wrap(~Model, scales = "free_y", ncol = 2)
    
    print(comparison_plot)
    cat(paste("Created", species_name, "comparison plot\n"))
  }
  
  return(combined_data)
}
#### INSERT THE PLOTTING CODE HERE ####
# Create plots using the final model importance results
Usnea_data <- create_final_model_plots("Usnea", varimp_clf_final_combined, varimp_qrf_final_combined)

# SAVE FINAL MODEL OUTPUTS
saveRDS(all_results_clf_final, "hurdle_rf_classifier_results_USNEA_top_predictors.rds")
saveRDS(all_results_qrf_final, "hurdle_qrf_abundance_results_USNEA_top_predictors.rds")
saveRDS(varimp_list_clf_final, "hurdle_varimp_rf_presence_USNEA_top_predictors.rds")
saveRDS(varimp_list_qrf_final, "hurdle_varimp_qrf_abundance_USNEA_top_predictors.rds")
saveRDS(predictions_all_final, "hurdle_predictions_all_USNEA_top_predictors.rds")
saveRDS(training_data_list_final, "hurdle_training_data_list_USNEA_top_predictors.rds")
saveRDS(test_data_list_final, "hurdle_test_data_list_USNEA_top_predictors.rds")


#### QUASI-POISSON MODEL COMPARISON ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RUNNING QUASI-POISSON MODEL FOR COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

library(MASS)     # For stepwise selection
library(car)      # For VIF
library(broom)    # For tidy model output
library(pROC)     # For AUC calculation

# Initialize storage for Quasi-Poisson results
all_results_qp <- list()
predictions_all_qp <- list()
model_summaries_qp <- list()

# Helper function to calculate Quasi-Poisson performance metrics
calc_qp_metrics <- function(predicted, observed, predicted_probs = NULL) {
  # For count data, we'll calculate both classification metrics (presence/absence) and regression metrics
  
  # Convert to presence/absence for classification metrics
  pred_presence <- ifelse(predicted > 0, 1, 0)
  obs_presence <- ifelse(observed > 0, 1, 0)
  
  # Classification metrics (presence/absence)
  acc <- mean(pred_presence == obs_presence)
  
  # Calculate confusion matrix elements
  tp <- sum(pred_presence == 1 & obs_presence == 1)
  tn <- sum(pred_presence == 0 & obs_presence == 0)
  fp <- sum(pred_presence == 1 & obs_presence == 0)
  fn <- sum(pred_presence == 0 & obs_presence == 1)
  
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  
  # AUC for presence/absence prediction using predicted counts as proxy
  auc <- tryCatch(pROC::auc(obs_presence, predicted, quiet = TRUE), 
                  error = function(e) NA)
  
  # Regression metrics for count data
  resids <- predicted - observed
  rmse <- sqrt(mean(resids^2, na.rm = TRUE))
  mae <- mean(abs(resids), na.rm = TRUE)
  bias <- mean(resids, na.rm = TRUE)
  
  # R-squared for count data
  r2 <- tryCatch({
    cor(observed, predicted, use = "complete.obs")^2
  }, error = function(e) NA)
  
  # Normalized MAE
  nmae <- ifelse(mean(observed, na.rm = TRUE) > 0, 
                 mae / mean(observed, na.rm = TRUE), NA)
  
  # Zero prediction metrics
  obs_zeros <- sum(observed == 0)
  pred_zeros <- sum(predicted == 0)
  zero_accuracy <- mean((observed == 0) == (predicted == 0))
  
  return(list(
    # Classification metrics
    acc = acc, sens = sens, spec = spec, auc = as.numeric(auc),
    # Regression metrics
    rmse = rmse, mae = mae, bias = bias, r2 = r2, nmae = nmae,
    # Zero prediction metrics
    obs_zeros = obs_zeros, pred_zeros = pred_zeros, zero_accuracy = zero_accuracy
  ))
}

cat("Running Quasi-Poisson model with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  set.seed(i)
  
  # Use same train-test splits as RF model for fair comparison
  if (i <= length(training_data_list_final) && i <= length(test_data_list_final)) {
    train <- training_data_list_final[[i]]
    test <- test_data_list_final[[i]]
  } else {
    # Fallback: create new split
    repeat {
      split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      # Check factor levels
      all_ok <- TRUE
      for (fac in factor_vars_presence) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
      if (all_ok && length(unique(train$Presence)) >= 2) break
    }
    
    # Align factor levels
    for (fac in factor_vars_presence) {
      if (fac %in% names(test)) {
        test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
      }
    }
  }
  
  ## Quasi-Poisson Model
  tryCatch({
    # Prepare count data - convert abundance to integer counts
    # You may need to scale/round your abundance data appropriately
    train$Usnea_count <- round(train$Usnea)
    test$Usnea_count <- round(test$Usnea)
    
    # Ensure non-negative integers
    train$Usnea_count <- pmax(0, train$Usnea_count)
    test$Usnea_count <- pmax(0, test$Usnea_count)
    
    # Create formula combining both presence and abundance predictors
    all_predictors <- unique(c(predictors_presence, predictors_abundance))
    qp_formula <- as.formula(paste("Usnea_count ~", paste(all_predictors, collapse = " + ")))
    
    # Fit Quasi-Poisson model
    qp_model <- glm(qp_formula, data = train, family = quasipoisson())
    
    # Check for convergence issues
    if (!qp_model$converged) {
      warning(paste("Quasi-Poisson model did not converge in iteration", i))
    }
    
    # Predict on test set
    pred_counts_qp <- predict(qp_model, newdata = test, type = "response")
    
    # Calculate metrics
    metrics_qp <- calc_qp_metrics(pred_counts_qp, test$Usnea_count)
    
    # Get dispersion parameter
    dispersion <- summary(qp_model)$dispersion
    
    all_results_qp[[i]] <- data.table(
      Iteration = i,
      # Classification metrics (presence/absence)
      Accuracy = metrics_qp$acc,
      Sensitivity = metrics_qp$sens,
      Specificity = metrics_qp$spec,
      AUC = metrics_qp$auc,
      # Regression metrics (count prediction)
      RMSE = metrics_qp$rmse,
      MAE = metrics_qp$mae,
      Bias = metrics_qp$bias,
      R_squared = metrics_qp$r2,
      nMAE = metrics_qp$nmae,
      # Zero prediction metrics
      Obs_Zeros = metrics_qp$obs_zeros,
      Pred_Zeros = metrics_qp$pred_zeros,
      Zero_Accuracy = metrics_qp$zero_accuracy,
      # Model diagnostics
      Dispersion = dispersion,
      Converged = qp_model$converged,
      N_Test = nrow(test),
      N_Train = nrow(train)
    )
    
    # Store predictions
    train_pred_counts <- predict(qp_model, newdata = train, type = "response")
    
    predictions_all_qp[[i]] <- rbind(
      data.table(set = "test", obs = test$Usnea_count, pred = pred_counts_qp, 
                 obs_orig = test$Usnea, iteration = i),
      data.table(set = "train", obs = train$Usnea_count, pred = train_pred_counts, 
                 obs_orig = train$Usnea, iteration = i)
    )
    
    # Store model for first iteration
    if (i == 1) {
      model_summaries_qp[["qp_model"]] <- list(
        model = qp_model,
        summary = summary(qp_model),
        aic = AIC(qp_model),
        formula = qp_formula,
        dispersion = dispersion,
        deviance = qp_model$deviance,
        null_deviance = qp_model$null.deviance
      )
    }
    
  }, error = function(e) {
    cat("Quasi-Poisson Error in iteration", i, ":", e$message, "\n")
    
    # Store NA results for failed iterations
    all_results_qp[[i]] <- data.table(
      Iteration = i, Accuracy = NA, Sensitivity = NA, Specificity = NA, AUC = NA,
      RMSE = NA, MAE = NA, Bias = NA, R_squared = NA, nMAE = NA,
      Obs_Zeros = NA, Pred_Zeros = NA, Zero_Accuracy = NA,
      Dispersion = NA, Converged = FALSE, N_Test = NA, N_Train = NA
    )
    
    predictions_all_qp[[i]] <- NULL
  })
  
  if (i %% 10 == 0) cat("Quasi-Poisson iteration", i, "completed\n")
}

# Combine Quasi-Poisson results
results_qp_combined <- rbindlist(all_results_qp)
predictions_qp_combined <- rbindlist(predictions_all_qp[!sapply(predictions_all_qp, is.null)])

#### MODEL COMPARISON ANALYSIS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RANDOM FOREST vs QUASI-POISSON COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# PRESENCE PREDICTION COMPARISON
cat("\n--- PRESENCE PREDICTION COMPARISON ---\n")
cat("Random Forest vs Quasi-Poisson\n\n")

rf_clf_summary <- results_clf_final_combined[, .(
  Model = "Random Forest",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

qp_clf_summary <- results_qp_combined[, .(
  Model = "Quasi-Poisson",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

presence_comparison <- rbind(rf_clf_summary, qp_clf_summary)
print(presence_comparison)

# ABUNDANCE COMPARISON
cat("\n--- ABUNDANCE PREDICTION COMPARISON ---\n")
cat("Quantile Random Forest vs Quasi-Poisson\n\n")

if (nrow(results_qrf_final_combined) > 0 && nrow(results_qp_combined) > 0) {
  rf_abundance_summary <- results_qrf_final_combined[, .(
    Model = "Quantile RF",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  qp_abundance_summary <- results_qp_combined[, .(
    Model = "Quasi-Poisson",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  abundance_comparison <- rbind(rf_abundance_summary, qp_abundance_summary)
  print(abundance_comparison)
  
  # Statistical tests for abundance metrics
  cat("\n--- STATISTICAL SIGNIFICANCE TESTS ---\n")
  
  # Presence metrics
  accuracy_test <- t.test(results_clf_final_combined$Accuracy, results_qp_combined$Accuracy)
  auc_test <- t.test(results_clf_final_combined$AUC, results_qp_combined$AUC)
  
  cat("Accuracy difference (RF - QP):\n")
  cat("  Mean difference:", round(accuracy_test$estimate[1] - accuracy_test$estimate[2], 4), "\n")
  cat("  p-value:", round(accuracy_test$p.value, 4), ifelse(accuracy_test$p.value < 0.05, "*", ""), "\n")
  
  cat("AUC difference (RF - QP):\n")
  cat("  Mean difference:", round(auc_test$estimate[1] - auc_test$estimate[2], 4), "\n")
  cat("  p-value:", round(auc_test$p.value, 4), ifelse(auc_test$p.value < 0.05, "*", ""), "\n")
  
  # Abundance metrics
  rf_rmse_valid <- results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]
  qp_rmse_valid <- results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)]
  
  if (length(rf_rmse_valid) > 10 && length(qp_rmse_valid) > 10) {
    rmse_test <- t.test(rf_rmse_valid, qp_rmse_valid)
    cat("RMSE difference (QRF - QP):\n")
    cat("  Mean difference:", round(rmse_test$estimate[1] - rmse_test$estimate[2], 4), "\n")
    cat("  p-value:", round(rmse_test$p.value, 4), ifelse(rmse_test$p.value < 0.05, "*", ""), "\n")
    
    rf_r2_valid <- results_qrf_final_combined$R_squared[!is.na(results_qrf_final_combined$R_squared)]
    qp_r2_valid <- results_qp_combined$R_squared[!is.na(results_qp_combined$R_squared)]
    
    if (length(rf_r2_valid) > 10 && length(qp_r2_valid) > 10) {  
      r2_test <- t.test(rf_r2_valid, qp_r2_valid)
      cat("R² difference (QRF - QP):\n")
      cat("  Mean difference:", round(r2_test$estimate[1] - r2_test$estimate[2], 4), "\n")
      cat("  p-value:", round(r2_test$p.value, 4), ifelse(r2_test$p.value < 0.05, "*", ""), "\n")
    }
  }
}

#### OVERDISPERSION ANALYSIS ####
cat("\n--- OVERDISPERSION ANALYSIS ---\n")
qp_dispersion_summary <- results_qp_combined[, .(
  Mean_Dispersion = round(mean(Dispersion, na.rm = TRUE), 3),
  SD_Dispersion = round(sd(Dispersion, na.rm = TRUE), 3),
  Min_Dispersion = round(min(Dispersion, na.rm = TRUE), 3),
  Max_Dispersion = round(max(Dispersion, na.rm = TRUE), 3),
  Convergence_Rate = round(mean(Converged, na.rm = TRUE), 3)
)]

cat("Overdispersion analysis:\n")
print(qp_dispersion_summary)

if (qp_dispersion_summary$Mean_Dispersion > 1.5) {
  cat("\nNote: Mean dispersion > 1.5 indicates overdispersion in the data\n")
} else if (qp_dispersion_summary$Mean_Dispersion < 0.8) {
  cat("\nNote: Mean dispersion < 0.8 indicates potential underdispersion\n")
} else {
  cat("\nNote: Dispersion close to 1, suggesting Poisson assumption may be reasonable\n")
}

#### MODEL INTERPRETATION ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== QUASI-POISSON MODEL INTERPRETATION ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

if (!is.null(model_summaries_qp[["qp_model"]])) {
  cat("\n--- QUASI-POISSON COEFFICIENTS ---\n")
  
  qp_summary <- model_summaries_qp[["qp_model"]]$summary
  
  # Model coefficients
  cat("\nQUASI-POISSON MODEL COEFFICIENTS:\n")
  coefs <- qp_summary$coefficients
  coef_df <- data.frame(
    Variable = rownames(coefs),
    Estimate = round(coefs[, "Estimate"], 4),
    Std_Error = round(coefs[, "Std. Error"], 4),
    t_value = round(coefs[, "t value"], 3),
    p_value = round(coefs[, "Pr(>|t|)"], 4),
    Significant = ifelse(coefs[, "Pr(>|t|)"] < 0.05, "*", "")
  )
  print(coef_df)
  
  cat("\nInterpretation: Coefficients represent log rate ratios\n")
  cat("Positive coefficients increase expected count, negative decrease it\n")
  
  # Model fit statistics
  cat("\nMODEL FIT STATISTICS:\n")
  cat("Dispersion parameter:", round(model_summaries_qp[["qp_model"]]$dispersion, 3), "\n")
  cat("Residual deviance:", round(model_summaries_qp[["qp_model"]]$deviance, 2), "\n")
  cat("Null deviance:", round(model_summaries_qp[["qp_model"]]$null_deviance, 2), "\n")
  
  # Pseudo R-squared
  pseudo_r2 <- 1 - (model_summaries_qp[["qp_model"]]$deviance / model_summaries_qp[["qp_model"]]$null_deviance)
  cat("Pseudo R-squared:", round(pseudo_r2, 3), "\n")
  
  cat("AIC:", round(model_summaries_qp[["qp_model"]]$aic, 2), "\n")
}

# RF Variable Importance comparison
cat("\n--- RANDOM FOREST VARIABLE IMPORTANCE ---\n")
cat("Top 10 variables for presence prediction:\n")
print(head(presence_final_importance[, c("Variable", "Importance_pct")], 10))

if (nrow(abundance_final_importance) > 0) {
  cat("\nTop 10 variables for abundance prediction:\n")
  print(head(abundance_final_importance[, c("Variable", "Importance_pct")], 10))
}

#### PREDICTION COMPARISON PLOTS ####
cat("\n=== CREATING COMPARISON PLOTS ===\n")

if (nrow(predictions_qp_combined) > 0 && nrow(predictions_final_combined) > 0) {
  # Combine predictions for plotting (using original scale)
  rf_pred_plot <- predictions_final_combined[set == "test", .(obs, pred, model = "Random Forest")]
  qp_pred_plot <- predictions_qp_combined[set == "test", .(obs = obs_orig, pred, model = "Quasi-Poisson")]
  
  combined_pred_plot <- rbind(rf_pred_plot, qp_pred_plot)
  
  # Scatter plot comparison
  p_comparison <- ggplot(combined_pred_plot, aes(x = obs, y = pred, color = model)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~model) +
    labs(
      title = "Usnea Abundance: Random Forest vs Quasi-Poisson",
      subtitle = "Test set predictions (perfect predictions would fall on diagonal line)",
      x = "Observed Abundance",
      y = "Predicted Abundance",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Random Forest" = "#00A08A", "Quasi-Poisson" = "#F21A00"))
  
  print(p_comparison)
  
  # Performance metrics comparison plot
  metrics_comparison <- rbind(
    data.table(Model = "Random Forest", 
               RMSE = results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]),
    data.table(Model = "Quasi-Poisson", 
               RMSE = results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)])
  )
  
  if (nrow(metrics_comparison) > 0) {
    p_metrics <- ggplot(metrics_comparison, aes(x = Model, y = RMSE, fill = Model)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = "RMSE Distribution: Random Forest vs Quasi-Poisson",
        subtitle = "Lower RMSE indicates better performance",
        y = "Root Mean Square Error"
      ) +
      theme_minimal() +
      scale_fill_manual(values = c("Random Forest" = "#00A08A", "Quasi-Poisson" = "#F21A00"))
    
    print(p_metrics)
  }
  
  # Residuals vs fitted plot for Quasi-Poisson
  if (nrow(predictions_qp_combined[set == "test"]) > 0) {
    qp_test_data <- predictions_qp_combined[set == "test"]
    qp_test_data[, residuals := obs_orig - pred]
    
    p_residuals <- ggplot(qp_test_data, aes(x = pred, y = residuals)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue") +
      labs(
        title = "Quasi-Poisson Model: Residuals vs Fitted Values",
        subtitle = "Test set predictions - should show random scatter around zero",
        x = "Fitted Values",
        y = "Residuals"
      ) +
      theme_minimal()
    
    print(p_residuals)
  }
}

#### SUMMARY AND RECOMMENDATIONS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== SUMMARY AND RECOMMENDATIONS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\nMODEL COMPARISON:\n")
cat("- Random Forest: Non-parametric, handles complex interactions, ensemble method\n")
cat("- Quasi-Poisson: Parametric, accounts for overdispersion, interpretable coefficients\n")

cat("\nQUASI-POISSON ADVANTAGES:\n")
cat("- Handles overdispersion in count data (dispersion parameter ≠ 1)\n")
cat("- More flexible than standard Poisson for overdispersed data\n")
cat("- Provides interpretable coefficients (log rate ratios)\n")
cat("- Well-established statistical foundation for count data\n")
cat("- Uses same mean-variance relationship as Poisson but allows different variance\n")

cat("\nPERFORMANCE SUMMARY:\n")
if (exists("presence_comparison")) {
  rf_better_acc <- presence_comparison$Mean_Accuracy[1] > presence_comparison$Mean_Accuracy[2]
  rf_better_auc <- presence_comparison$Mean_AUC[1] > presence_comparison$Mean_AUC[2]
  
  cat("Presence prediction:\n")
  cat("  - Better accuracy:", ifelse(rf_better_acc, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Better AUC:", ifelse(rf_better_auc, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("abundance_comparison") && nrow(abundance_comparison) > 0) {
  rf_better_rmse <- abundance_comparison$Mean_RMSE[1] < abundance_comparison$Mean_RMSE[2]
  rf_better_r2 <- abundance_comparison$Mean_R2[1] > abundance_comparison$Mean_R2[2]
  
  cat("Count prediction:\n")
  cat("  - Lower RMSE (better):", ifelse(rf_better_rmse, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Higher R² (better):", ifelse(rf_better_r2, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("qp_dispersion_summary")) {
  cat("Overdispersion:\n")
  cat("  - Mean dispersion parameter:", qp_dispersion_summary$Mean_Dispersion, "\n")
  if (qp_dispersion_summary$Mean_Dispersion > 1.2) {
    cat("  - Data shows overdispersion, Quasi-Poisson is appropriate\n")
  }
}

cat("\nRECOMMENDATION:\n")
cat("Choose Random Forest if: Complex non-linear relationships, prediction accuracy priority\n")
cat("Choose Quasi-Poisson if: Need interpretable coefficients, understand linear relationships\n")
cat("Consider Quasi-Poisson if: Count data with overdispersion, parametric approach preferred\n")





########################################## SANIONIA #######################################################################################

rm(list = ls())

#### Sanionia ####
### SANIONIA TRANSFORMATIONS FOR HURDLE MODEL 
library(randomForest)
library(quantregForest)
library(data.table)
library(dplyr)
library(rsample)
library(isotone)
library(pROC)
library(pdp)
library(clipr)
library(ggplot2)

select <- dplyr::select

# Load data
df <- read.csv("Sanionia.csv")
df <- df %>% mutate(across(where(is.character), as.factor))
df$Presence <- as.factor(ifelse(df$Sanionia > 0, 1, 0))

# Define transformation functions
transformations <- list(
  none = list(
    transform = function(x) x,
    inverse = function(x) x,
    name = "None"
  ),
  sqrt = list(
    transform = function(x) sqrt(pmax(x, 0)),
    inverse = function(x) pmax(x^2, 0),
    name = "Square Root"
  ),
  log = list(
    transform = function(x) log(pmax(x, 0.001)),  # Add small constant to avoid log(0)
    inverse = function(x) exp(x),
    name = "Log"
  ),
  asinh = list(
    transform = function(x) asinh(x),
    inverse = function(x) sinh(x),
    name = "Asinh"
  )
)

predictors <- c("PDIS_CLASS", "ASP", "SLOPE_CLASS", "MORF_CLASS", "CD_CLASS", "WEXP",
                "ALTITUDE", "DISTCOAST", "MXCURV", "VDEP", "WEFF", "MBI", "Bio15_81", "PNCURV",
                "TRI", "TTCURV", "FLAC", "TPW", "FCF", "VDISCN", "SVF",
                "Easting", "Northing")

factor_vars <- names(df[, predictors])[sapply(df[, predictors], is.factor)]

n_iter <- 100

# Initialize storage for all transformations
all_results_clf <- list()
all_results_qrf <- list()
varimp_list_clf <- list()
varimp_list_qrf <- list()
predictions_all <- list()
training_data_list <- list()
test_data_list <- list()

for (trans_name in names(transformations)) {
  trans_func <- transformations[[trans_name]]
  
  cat("Testing transformation:", trans_func$name, "\n")
  
  # Transform the response variable
  df_trans <- df
  df_trans$Sanionia_trans <- trans_func$transform(df$Sanionia)
  
  all_results_clf[[trans_name]] <- list()
  all_results_qrf[[trans_name]] <- list()
  varimp_list_clf[[trans_name]] <- list()
  varimp_list_qrf[[trans_name]] <- list()
  predictions_all[[trans_name]] <- list()
  training_data_list[[trans_name]] <- list()
  test_data_list[[trans_name]] <- list()
  
  for (i in 1:n_iter) {
    set.seed(i)
    repeat {
      split <- initial_split(df_trans, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      all_ok <- TRUE
      for (fac in factor_vars) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) {
          all_ok <- FALSE
          break
        }
      }
      if (all_ok) break
    }
    
    for (fac in factor_vars) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
    
    ## STEP 1: RF Classifier (same for all transformations)
    rf_clf <- randomForest(
      x = train[, predictors],
      y = train$Presence,
      ntree = 500,
      importance = TRUE,
      keep.forest = TRUE
    )
    
    pred_probs <- predict(rf_clf, newdata = test[, predictors], type = "prob")[, "1"]
    pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
    obs_class <- test$Presence
    
    # Metrics for RF classifier
    acc <- mean(pred_class == obs_class)
    sens <- mean(pred_class == "1" & obs_class == "1")
    spec <- mean(pred_class == "0" & obs_class == "0")
    logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
    brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
    auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
    
    spread_clf <- abs(0.5 - pred_probs)  # proxy for uncertainty
    
    all_results_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
      LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
      MedianSpread = median(spread_clf)
    )
    
    imp_clf <- importance(rf_clf)
    varimp_list_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Variable = rownames(imp_clf), 
      Importance = imp_clf[, 1], 
      Iteration = i, 
      Model = "RF_Presence"
    )
    
    ## STEP 2: QRF for predicted present cases (with transformation)
    test_presence_idx <- which(pred_class == "1")
    test_present_data <- test[test_presence_idx, ]
    train_present_data <- train[train$Presence == "1", ]
    
    if (nrow(test_present_data) > 10) {
      qrf <- quantregForest(
        x = train_present_data[, predictors],
        y = train_present_data$Sanionia_trans,  # Use transformed response
        ntree = 500
      )
      
      preds_test_trans <- predict(qrf, newdata = test_present_data[, predictors], what = c(0.025, 0.5, 0.975))
      
      # Apply isotonic regression on transformed scale
      iso_model <- isoreg(preds_test_trans[, 2], test_present_data$Sanionia_trans)
      iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
      calibrated_test_trans <- iso_fun(preds_test_trans[, 2])
      
      # Back-transform predictions and bounds
      preds_test_lower <- trans_func$inverse(preds_test_trans[, 1])
      preds_test_median <- trans_func$inverse(calibrated_test_trans)
      preds_test_upper <- trans_func$inverse(preds_test_trans[, 3])
      
      # Ensure predictions are within reasonable bounds
      preds_test_lower <- pmax(pmin(preds_test_lower, 100), 0)
      preds_test_median <- pmax(pmin(preds_test_median, 100), 0)
      preds_test_upper <- pmax(pmin(preds_test_upper, 100), 0)
      
      spread_test <- preds_test_upper - preds_test_lower
      resids_test <- preds_test_median - test_present_data$Sanionia
      
      rmse <- sqrt(mean(resids_test^2))
      mae <- mean(abs(resids_test))
      bias <- mean(resids_test)
      r2 <- cor(test_present_data$Sanionia, preds_test_median)^2
      nmae <- mae / mean(test_present_data$Sanionia)
      coverage <- mean(test_present_data$Sanionia >= preds_test_lower & test_present_data$Sanionia <= preds_test_upper)
      median_spread <- median(spread_test)
      low_bias <- mean(preds_test_lower - test_present_data$Sanionia)
      high_bias <- mean(preds_test_upper - test_present_data$Sanionia)
      
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
        R_squared = r2, nMAE = nmae, Coverage_95PI = coverage,
        Median_Spread_95PI = median_spread,
        Lower_Q2.5_Bias = low_bias,
        Upper_Q97.5_Bias = high_bias
      )
      
      # Store predictions (back-transformed for training data too)
      train_preds_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.5)
      train_spread_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.975) -
        predict(qrf, newdata = train_present_data[, predictors], what = 0.025)
      
      predictions_all[[trans_name]][[i]] <- rbind(
        data.table(set = "test", obs = test_present_data$Sanionia, pred = preds_test_median,
                   spread = spread_test, iteration = i, transformation = trans_func$name),
        data.table(set = "train", obs = train_present_data$Sanionia,
                   pred = trans_func$inverse(train_preds_trans),
                   spread = trans_func$inverse(train_present_data$Sanionia_trans + train_spread_trans/2) - 
                     trans_func$inverse(train_present_data$Sanionia_trans - train_spread_trans/2),
                   iteration = i, transformation = trans_func$name)
      )
      
      imp_qrf <- importance(qrf)
      varimp_list_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Variable = rownames(imp_qrf), 
        Importance = imp_qrf[, 1], 
        Iteration = i, 
        Model = "QRF_Abundance"
      )
    } else {
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
        R_squared = NA, nMAE = NA, Coverage_95PI = NA,
        Median_Spread_95PI = NA,
        Lower_Q2.5_Bias = NA,
        Upper_Q97.5_Bias = NA
      )
      predictions_all[[trans_name]][[i]] <- NULL
      varimp_list_qrf[[trans_name]][[i]] <- NULL
    }
    
    training_data_list[[trans_name]][[i]] <- train
    test_data_list[[trans_name]][[i]] <- test
    
    if (i %% 10 == 0) cat("Transformation:", trans_func$name, "- Iteration", i, "\n")
  }
}


# Combine results across transformations
results_clf_combined <- rbindlist(lapply(all_results_clf, function(x) rbindlist(x)))
results_qrf_combined <- rbindlist(lapply(all_results_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_clf_combined <- rbindlist(lapply(varimp_list_clf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_qrf_combined <- rbindlist(lapply(varimp_list_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
predictions_combined <- rbindlist(lapply(predictions_all, function(x) rbindlist(x[!sapply(x, is.null)])))

# Summary statistics by transformation
cat("\n=== CLASSIFICATION RESULTS SUMMARY ===\n")
clf_summary <- results_clf_combined[, .(
  Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
  SD_Accuracy = sd(Accuracy, na.rm = TRUE),
  Mean_AUC = mean(AUC, na.rm = TRUE),
  SD_AUC = sd(AUC, na.rm = TRUE),
  Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
  Mean_Specificity = mean(Specificity, na.rm = TRUE)
), by = Transformation]
print(clf_summary)

cat("\n=== QUANTILE REGRESSION RESULTS SUMMARY ===\n")
qrf_summary <- results_qrf_combined[, .(
  Mean_RMSE = mean(RMSE, na.rm = TRUE),
  SD_RMSE = sd(RMSE, na.rm = TRUE),
  Mean_MAE = mean(MAE, na.rm = TRUE),
  SD_MAE = sd(MAE, na.rm = TRUE),
  Mean_R2 = mean(R_squared, na.rm = TRUE),
  SD_R2 = sd(R_squared, na.rm = TRUE),
  Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
  Mean_Bias = mean(Bias, na.rm = TRUE)
), by = Transformation]
print(qrf_summary)

# 1. PREPARE DATA FOR PLOTTING

# Spread plot data (prediction interval width vs observed)
spread_plot_data <- predictions_combined[set == "test", .(
  obs = obs,
  spread = spread,
  trans = transformation,
  iteration = iteration
)]

# Observed vs predicted data
pred_obs_data <- predictions_combined[set == "test", .(
  obs = obs,
  pred = pred,
  trans = transformation,
  iteration = iteration
)]

# Calculate residuals for each transformation
residuals_data <- predictions_combined[set == "test", .(
  pred = pred,
  residuals = obs - pred,
  trans = transformation,
  iteration = iteration
)]

# 2. SPREAD PLOT (Uncertainty vs Observed)
p1 <- ggplot(spread_plot_data, aes(x = obs, y = spread, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Observed Sanionia Abundance", y = "Prediction Interval Width",
       title = "Uncertainty vs Observed (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# 3. VARIABLE IMPORTANCE ANALYSIS
cat("\n=== VARIABLE IMPORTANCE ANALYSIS ===\n")

# Calculate variable importance summary
varimp_summary <- varimp_qrf_combined %>%
  group_by(Transformation, Variable) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 4),
            .groups = "drop") %>%
  arrange(Transformation, desc(MedianImportance))

# Calculate percentage importance
varimp_summary <- varimp_summary %>%
  group_by(Transformation) %>%
  mutate(Importance_pct = 100 * MedianImportance / sum(MedianImportance)) %>%
  ungroup()

# For "None" transformation (equivalent to "raw")
none_imp <- varimp_summary %>%
  filter(Transformation == "None") %>%
  mutate(Cumulative_Importance = round(cumsum(Importance_pct), 1)) %>%
  arrange(desc(MedianImportance))

cat("Variable Importance for No Transformation (Raw):\n")
print(none_imp)
clipr::write_clip(none_imp)

# Variable importance comparison plot
p2 <- varimp_summary %>%
  group_by(Transformation) %>%
  slice_head(n = 10) %>%  # Top 10 variables per transformation
  ggplot(aes(x = reorder(Variable, MedianImportance), y = MedianImportance, fill = Transformation)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Variable", y = "Median Importance", 
       title = "Top 10 Variable Importance by Transformation") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  facet_wrap(~Transformation, scales = "free_y")

print(p2)

# 4. OBSERVED VS PREDICTED PLOT
p3 <- ggplot(pred_obs_data, aes(x = obs, y = pred, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed Sanionia Abundance", y = "Predicted Sanionia Abundance",
       title = "Observed vs Predicted (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# 5. RESIDUALS VS PREDICTED (for each transformation)
p4 <- ggplot(residuals_data, aes(x = pred, y = residuals)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~trans, scales = "free") +
  theme_minimal() +
  labs(title = "Residuals vs Predicted (by Transformation)",
       x = "Predicted Sanionia Abundance",
       y = "Residuals (Observed - Predicted)")

print(p4)

# 6. INDIVIDUAL TRANSFORMATION PLOTS

# Function to create plots for a specific transformation
create_transformation_plots <- function(trans_name) {
  
  # Filter data for specific transformation
  trans_spread <- spread_plot_data[trans == trans_name]
  trans_pred_obs <- pred_obs_data[trans == trans_name]
  trans_residuals <- residuals_data[trans == trans_name]
  
  cat(paste("\n=== PLOTS FOR", trans_name, "TRANSFORMATION ===\n"))
  
  # Residuals vs Predicted for specific transformation
  p_resid <- ggplot(trans_residuals, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Residuals vs Predicted -", trans_name, "Transformation"),
         x = "Predicted Sanionia Abundance",
         y = "Residuals (Observed - Predicted)")
  
  # Prediction Interval Width vs Observed for specific transformation
  p_spread <- ggplot(trans_spread, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    theme_minimal() +
    labs(title = paste("Prediction Interval Width vs Observed -", trans_name, "Transformation"),
         x = "Observed Sanionia Abundance",
         y = "95% Prediction Interval Width")
  
  print(p_resid)
  print(p_spread)
  
  return(list(residuals = p_resid, spread = p_spread))
}

# Create individual plots for each transformation
transformation_plots <- list()
for (trans in unique(spread_plot_data$trans)) {
  transformation_plots[[trans]] <- create_transformation_plots(trans)
}

# 7. PERFORMANCE SUMMARY COMPARISON PLOT
perf_summary_long <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI, Bias) %>%
  group_by(Transformation) %>%
  summarize(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Mean_MAE = mean(MAE, na.rm = TRUE),
    Mean_R2 = mean(R_squared, na.rm = TRUE),
    Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
    Mean_Bias = abs(mean(Bias, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p5 <- ggplot(perf_summary_long, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_col() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Model Performance Comparison Across Transformations",
       x = "Transformation", y = "Mean Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p5)

# 8. BOXPLOT OF PERFORMANCE METRICS
perf_boxplot_data <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p6 <- ggplot(perf_boxplot_data, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Distribution of Performance Metrics Across Iterations",
       x = "Transformation", y = "Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p6)

# 9. CORRELATION BETWEEN OBSERVED AND PREDICTED BY TRANSFORMATION
cor_summary <- pred_obs_data %>%
  group_by(trans) %>%
  summarize(
    correlation = cor(obs, pred, use = "complete.obs"),
    rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    mae = mean(abs(obs - pred), na.rm = TRUE),
    .groups = "drop")
print(cor_summary)

#### EXTRACT TOP PREDICTORS FROM TRANSFORMATION ANALYSIS ####

# Extract top predictors from presence model (RF Classifier) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_presence <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

# Extract top predictors from abundance model (QRF) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_abundance <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

cat("\n=== TOP PREDICTORS EXTRACTED FROM RAW TRANSFORMATION (70% IMPORTANCE) ===\n")
cat("Top predictors for PRESENCE model (explaining 70% of importance):\n")
print(top_predictors_presence)
cat("Number of presence predictors:", length(top_predictors_presence), "\n")

cat("\nTop predictors for ABUNDANCE model (explaining 70% of importance):\n")
print(top_predictors_abundance)
cat("Number of abundance predictors:", length(top_predictors_abundance), "\n")

# Optional: Show the importance breakdown
cat("\n=== IMPORTANCE BREAKDOWN FOR SELECTED PREDICTORS ===\n")

# Presence model importance breakdown
presence_importance_breakdown <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_presence
  ) %>%
  filter(Selected)

cat("PRESENCE model - Selected variables and their importance:\n")
print(presence_importance_breakdown)

# Abundance model importance breakdown
abundance_importance_breakdown <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_abundance
  ) %>%
  filter(Selected)

cat("\nABUNDANCE model - Selected variables and their importance:\n")
print(abundance_importance_breakdown)


#### HURDLE FINAL MODEL USING TOP PREDICTORS ####

# Use the extracted top predictors for each model
predictors_presence <- c("ALTITUDE", "Easting", "WEFF", "WEXP", "DISTCOAST", "Northing", "ASP", "VDEP", "Bio15_81", "FLAC", "TRI")
predictors_abundance <- c("VDEP", "ASP", "DISTCOAST", "SVF", "WEFF", "Easting", "FLAC", "TRI", "ALTITUDE", "Northing")

# Get factor variables for each predictor set
factor_vars_presence <- names(df[, predictors_presence])[sapply(df[, predictors_presence], is.factor)]
factor_vars_abundance <- names(df[, predictors_abundance])[sapply(df[, predictors_abundance], is.factor)]

n_iter <- 100
all_results_clf_final <- list()
all_results_qrf_final <- list()
varimp_list_clf_final <- list()
varimp_list_qrf_final <- list()
predictions_all_final <- list()
training_data_list_final <- list()
test_data_list_final <- list()

cat("\n=== RUNNING FINAL MODEL WITH TOP PREDICTORS ===\n")

for (i in 1:n_iter) {
  set.seed(i)
  repeat {
    split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
    train <- training(split)
    test <- testing(split)
    
    # Check factor levels for BOTH predictor sets
    all_ok <- TRUE
    for (fac in factor_vars_presence) {
      train_levels <- levels(train[[fac]])
      test_levels <- levels(factor(test[[fac]]))
      if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
    }
    if (all_ok) {
      for (fac in factor_vars_abundance) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
      }
    }
    if (all_ok) break
  }
  
  # Set factor levels
  for (fac in factor_vars_presence) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  for (fac in factor_vars_abundance) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  
  ## STEP 1: RF Classifier
  rf_clf <- randomForest(
    x = train[, predictors_presence],
    y = train$Presence,
    ntree = 500,
    importance = TRUE,
    keep.forest = TRUE
  )
  
  pred_probs <- predict(rf_clf, newdata = test[, predictors_presence], type = "prob")[, "1"]
  pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
  obs_class <- test$Presence
  
  # Metrics for classifier
  acc <- mean(pred_class == obs_class)
  sens <- mean(pred_class == "1" & obs_class == "1")
  spec <- mean(pred_class == "0" & obs_class == "0")
  logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
  brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
  auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
  spread_clf <- abs(0.5 - pred_probs)
  
  all_results_clf_final[[i]] <- data.table(
    Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
    LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
    MedianSpread = median(spread_clf)
  )
  
  varimp_list_clf_final[[i]] <- data.table(
    Variable = rownames(importance(rf_clf)),
    Importance = importance(rf_clf)[,1],
    Iteration = i, Model = "RF_Presence"
  )
  
  ## STEP 2: QRF for abundance
  test_presence_idx <- which(pred_class == "1")
  test_present_data <- test[test_presence_idx, ]
  train_present_data <- train[train$Presence == "1", ]
  
  if (nrow(test_present_data) > 10) {
    qrf <- quantregForest(
      x = train_present_data[, predictors_abundance],
      y = train_present_data$Sanionia,
      ntree = 500
    )
    
    preds_test <- predict(qrf, newdata = test_present_data[, predictors_abundance], what = c(0.025, 0.5, 0.975))
    iso_model <- isoreg(preds_test[, 2], test_present_data$Sanionia)
    iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
    calibrated_test <- pmax(pmin(iso_fun(preds_test[, 2]), 100), 0)
    
    spread_test <- preds_test[, 3] - preds_test[, 1]
    resids_test <- calibrated_test - test_present_data$Sanionia
    
    rmse <- sqrt(mean(resids_test^2))
    mae <- mean(abs(resids_test))
    bias <- mean(resids_test)
    r2 <- cor(test_present_data$Sanionia, calibrated_test)^2
    nrmse <- rmse / mean(test_present_data$Sanionia, na.rm = TRUE)
    nrmse_median <- rmse / median(test_present_data$Sanionia, na.rm = TRUE)
    nmae <- mae / mean(test_present_data$Sanionia, na.rm = TRUE)
    nmae_median <- mae / median(test_present_data$Sanionia, na.rm = TRUE)
    coverage <- mean(test_present_data$Sanionia >= preds_test[,1] & test_present_data$Sanionia <= preds_test[,3])
    median_spread <- median(spread_test)
    low_bias <- mean(preds_test[,1] - test_present_data$Sanionia)
    high_bias <- mean(preds_test[,3] - test_present_data$Sanionia)
    
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
      R_squared = r2,
      nRMSE = nrmse,
      nRMSE_Median = nrmse_median,
      nMAE = nmae,
      nMAE_Median = nmae_median,
      Coverage_95PI = coverage,
      Median_Spread_95PI = median_spread,
      Lower_Q2.5_Bias = low_bias,
      Upper_Q97.5_Bias = high_bias
    )
    
    predictions_all_final[[i]] <- rbind(
      data.table(set = "test", obs = test_present_data$Sanionia, pred = calibrated_test,
                 spread = spread_test, iteration = i),
      data.table(set = "train", obs = train_present_data$Sanionia,
                 pred = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.5),
                 spread = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.975) -
                   predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.025),
                 iteration = i)
    )
    
    varimp_list_qrf_final[[i]] <- data.table(
      Variable = rownames(importance(qrf)),
      Importance = importance(qrf)[,1],
      Iteration = i, Model = "QRF_Abundance"
    )
  } else {
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
      R_squared = NA, nRMSE = NA, nRMSE_Median = NA, nMAE = NA, nMAE_Median = NA,
      Coverage_95PI = NA, Median_Spread_95PI = NA, Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA
    )
    predictions_all_final[[i]] <- NULL
    varimp_list_qrf_final[[i]] <- NULL
  }
  
  training_data_list_final[[i]] <- train
  test_data_list_final[[i]] <- test
  
  if (i %% 10 == 0) cat("Final model iteration", i, "\n")
}

# Combine results
results_clf_final_combined <- rbindlist(all_results_clf_final)
results_qrf_final_combined <- rbindlist(all_results_qrf_final[!sapply(all_results_qrf_final, is.null)])
varimp_clf_final_combined <- rbindlist(varimp_list_clf_final[!sapply(varimp_list_clf_final, is.null)])
varimp_qrf_final_combined <- rbindlist(varimp_list_qrf_final[!sapply(varimp_list_qrf_final, is.null)])
predictions_final_combined <- rbindlist(predictions_all_final[!sapply(predictions_all_final, is.null)])

# CLASSIFICATION SUMMARY
clf_final_summary <- results_clf_final_combined[, .(
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 3),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 3),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 3),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 3),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 3),
  SD_Sensitivity = round(sd(Sensitivity, na.rm = TRUE), 3),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 3),
  SD_Specificity = round(sd(Specificity, na.rm = TRUE), 3),
  Mean_LogLoss = round(mean(LogLoss, na.rm = TRUE), 3),
  SD_LogLoss = round(sd(LogLoss, na.rm = TRUE), 3),
  Mean_BrierScore = round(mean(BrierScore, na.rm = TRUE), 3),
  SD_BrierScore = round(sd(BrierScore, na.rm = TRUE), 3),
  Mean_MedianSpread = round(mean(MedianSpread, na.rm = TRUE), 3),
  SD_MedianSpread = round(sd(MedianSpread, na.rm = TRUE), 3)
)]
print(clf_final_summary)

qrf_final_summary <- results_qrf_final_combined[, .(
  Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
  SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
  Median_RMSE = round(median(RMSE, na.rm = TRUE), 3),
  Mean_MAE = round(mean(MAE, na.rm = TRUE), 3),
  SD_MAE = round(sd(MAE, na.rm = TRUE), 3),
  Median_MAE = round(median(MAE, na.rm = TRUE), 3),
  Mean_R2 = round(mean(R_squared, na.rm = TRUE), 3),
  SD_R2 = round(sd(R_squared, na.rm = TRUE), 3),
  Median_R2 = round(median(R_squared, na.rm = TRUE), 3),
  Median_nRMSE = round(median(nRMSE_Median, na.rm = TRUE), 3),
  IQR_nRMSE   = round(IQR(nRMSE_Median, na.rm = TRUE), 3),
  Mean_Bias = round(mean(Bias, na.rm = TRUE), 3),
  SD_Bias = round(sd(Bias, na.rm = TRUE), 3),
  Median_Bias = round(median(Bias, na.rm = TRUE), 3),
  Mean_Coverage_95PI = round(mean(Coverage_95PI, na.rm = TRUE), 3),
  SD_Coverage_95PI = round(sd(Coverage_95PI, na.rm = TRUE), 3),
  Median_Coverage_95PI = round(median(Coverage_95PI, na.rm = TRUE), 3),
  Mean_Median_Spread_95PI = round(mean(Median_Spread_95PI, na.rm = TRUE), 3),
  SD_Median_Spread_95PI = round(sd(Median_Spread_95PI, na.rm = TRUE), 3),
  Median_Spread_95PI = round(median(Median_Spread_95PI, na.rm = TRUE), 3),
  Mean_Lower_Q2.5_Bias = round(mean(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  SD_Lower_Q2.5_Bias = round(sd(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Median_Lower_Q2.5_Bias = round(median(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Mean_Upper_Q97.5_Bias = round(mean(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  SD_Upper_Q97.5_Bias = round(sd(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  Median_Upper_Q97.5_Bias = round(median(Upper_Q97.5_Bias, na.rm = TRUE), 3)
)]
print(qrf_final_summary)


#### VARIABLE IMPORTANCE ANALYSIS FOR FINAL MODELS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== VARIABLE IMPORTANCE FROM FINAL MODELS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create comprehensive importance tables from the final model runs
cat("\nTABLE 1: PRESENCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

presence_final_importance <- varimp_clf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(presence_final_importance)

cat("\nPRESENCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(presence_final_importance), "\n"))
cat(paste("- Top variable (", presence_final_importance$Variable[1], ") explains:", 
          presence_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", presence_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 2: ABUNDANCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

abundance_final_importance <- varimp_qrf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(abundance_final_importance)

cat("\nABUNDANCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(abundance_final_importance), "\n"))
cat(paste("- Top variable (", abundance_final_importance$Variable[1], ") explains:", 
          abundance_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", abundance_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

# Create a comparison table showing which variables are most important in each model
cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 3: COMPARISON - Top Variables in Each Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

# Get top 5 from each model for comparison
top_presence_vars <- presence_final_importance$Variable[1:min(5, nrow(presence_final_importance))]
top_abundance_vars <- abundance_final_importance$Variable[1:min(5, nrow(abundance_final_importance))]

# Create comparison dataframe
comparison_table <- data.frame(
  Rank = 1:5,
  Presence_Model = c(top_presence_vars, rep("", 5 - length(top_presence_vars)))[1:5],
  Presence_Importance = c(presence_final_importance$Importance_pct[1:min(5, nrow(presence_final_importance))], 
                          rep(NA, 5 - min(5, nrow(presence_final_importance))))[1:5],
  Abundance_Model = c(top_abundance_vars, rep("", 5 - length(top_abundance_vars)))[1:5],
  Abundance_Importance = c(abundance_final_importance$Importance_pct[1:min(5, nrow(abundance_final_importance))], 
                           rep(NA, 5 - min(5, nrow(abundance_final_importance))))[1:5]
)

print(comparison_table)

# Identify shared important variables
shared_vars <- intersect(top_presence_vars, top_abundance_vars)
if(length(shared_vars) > 0) {
  cat("\nShared important variables between models:", paste(shared_vars, collapse = ", "), "\n")
} else {
  cat("\nNo shared variables in top 5 between models\n")
}

#### PLOTTING FUNCTIONS - RUN THIS FIRST ####

# Define color palette for species
cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00")

library(ggplot2)
library(dplyr)

# Function to extract 70% importance variables from FINAL MODEL results
extract_final_70_percent <- function(varimp_final_combined, species_name, model_type) {
  
  # Calculate importance from final model runs and apply 70% threshold
  importance_data <- varimp_final_combined %>%
    group_by(Variable) %>%
    summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MedianImportance)) %>%
    mutate(
      Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
      Cumulative_pct = round(cumsum(Importance_pct), 2),
      Selected_70pct = Cumulative_pct <= 70,
      Rank = row_number(),
      Species = species_name,
      Model = model_type
    ) %>%
    filter(Selected_70pct) %>%  # Only keep variables that explain ≤70% cumulative importance
    select(Species, Model, Variable, Rank, Importance_pct, Cumulative_pct)
  
  return(importance_data)
}

create_individual_plot <- function(species_name, model_type, data_all) {
  
  data_subset <- data_all %>%
    filter(Species == species_name, Model == model_type) %>%
    arrange(desc(Rank))  # Most important at top
  
  if(nrow(data_subset) == 0) {
    cat(paste("No data found for", species_name, model_type, "\n"))
    return(NULL)
  }
  
  p <- ggplot(data_subset, aes(x = reorder(Variable, Rank), y = Importance_pct)) +
    geom_col(fill = cores[species_name], alpha = 0.8, color = "white", size = 0.5) +
    geom_text(aes(label = paste0(round(Importance_pct, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    labs(
      title = paste(species_name, "-", model_type, "Final Model (70% Variables)"),
      subtitle = paste("Variables explaining", round(max(data_subset$Cumulative_pct), 1), "% importance in FINAL model"),
      x = "Variables (Ranked by Final Model Importance)",
      y = "Individual Importance (%)",
      caption = paste("n =", nrow(data_subset), "variables from final model")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
      plot.subtitle = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

create_final_model_plots <- function(species_name, varimp_clf_final, varimp_qrf_final) {
  
  cat(paste("\n=== Creating plots for", species_name, "final model ===\n"))
  
  # Extract 70% variables from final models
  presence_data <- extract_final_70_percent(varimp_clf_final, species_name, "Presence")
  abundance_data <- extract_final_70_percent(varimp_qrf_final, species_name, "Abundance")
  
  # Combine data
  combined_data <- rbind(presence_data, abundance_data)
  
  # Create individual plots
  presence_plot <- create_individual_plot(species_name, "Presence", combined_data)
  abundance_plot <- create_individual_plot(species_name, "Abundance", combined_data)
  
  # Print plots
  if(!is.null(presence_plot)) {
    print(presence_plot)
    cat(paste("Created", species_name, "presence plot with", nrow(presence_data), "variables\n"))
  }
  
  if(!is.null(abundance_plot)) {
    print(abundance_plot) 
    cat(paste("Created", species_name, "abundance plot with", nrow(abundance_data), "variables\n"))
  }
  
  # Create comparison plot
  if(nrow(combined_data) > 0) {
    comparison_plot <- ggplot(combined_data, aes(x = reorder(Variable, Rank), y = Importance_pct, fill = Model)) +
      geom_col(position = "dodge", alpha = 0.8, color = "white", size = 0.5) +
      coord_flip() +
      labs(
        title = paste(species_name, "- Final Model: Presence vs Abundance (70% Variables)"),
        subtitle = "Variable importance from final models (70% threshold applied)",
        x = "Variables",
        y = "Individual Importance (%)",
        fill = "Model Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
        plot.subtitle = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      ) +
      scale_fill_manual(values = c("Presence" = cores[species_name], "Abundance" = adjustcolor(cores[species_name], alpha.f = 0.6))) +
      facet_wrap(~Model, scales = "free_y", ncol = 2)
    
    print(comparison_plot)
    cat(paste("Created", species_name, "comparison plot\n"))
  }
  
  return(combined_data)
}
#### INSERT THE PLOTTING CODE HERE ####
# Create plots using the final model importance results
Sanionia_data <- create_final_model_plots("Sanionia", varimp_clf_final_combined, varimp_qrf_final_combined)


# SAVE FINAL MODEL OUTPUTS
saveRDS(all_results_clf_final, "hurdle_rf_classifier_results_SANIONIA_top_predictors.rds")
saveRDS(all_results_qrf_final, "hurdle_qrf_abundance_results_SANIONIA_top_predictors.rds")
saveRDS(varimp_list_clf_final, "hurdle_varimp_rf_presence_SANIONIA_top_predictors.rds")
saveRDS(varimp_list_qrf_final, "hurdle_varimp_qrf_abundance_SANIONIA_top_predictors.rds")
saveRDS(predictions_all_final, "hurdle_predictions_all_SANIONIA_top_predictors.rds")
saveRDS(training_data_list_final, "hurdle_training_data_list_SANIONIA_top_predictors.rds")
saveRDS(test_data_list_final, "hurdle_test_data_list_SANIONIA_top_predictors.rds")


#### QUASI-POISSON MODEL COMPARISON ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RUNNING QUASI-POISSON MODEL FOR COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

library(MASS)     # For stepwise selection
library(car)      # For VIF
library(broom)    # For tidy model output
library(pROC)     # For AUC calculation

# Use Sanionia predictors 
predictors_presence <- c("ALTITUDE", "Easting", "WEFF", "WEXP", "DISTCOAST", "Northing", "ASP", "VDEP", "Bio15_81", "FLAC", "TRI")
predictors_abundance <- c("VDEP", "ASP", "DISTCOAST", "SVF", "WEFF", "Easting", "FLAC", "TRI", "ALTITUDE", "Northing")

# Initialize storage for Quasi-Poisson results
all_results_qp <- list()
predictions_all_qp <- list()
model_summaries_qp <- list()

# Helper function to calculate Quasi-Poisson performance metrics
calc_qp_metrics <- function(predicted, observed, predicted_probs = NULL) {
  # For count data, we'll calculate both classification metrics (presence/absence) and regression metrics
  
  # Convert to presence/absence for classification metrics
  pred_presence <- ifelse(predicted > 0, 1, 0)
  obs_presence <- ifelse(observed > 0, 1, 0)
  
  # Classification metrics (presence/absence)
  acc <- mean(pred_presence == obs_presence)
  
  # Calculate confusion matrix elements
  tp <- sum(pred_presence == 1 & obs_presence == 1)
  tn <- sum(pred_presence == 0 & obs_presence == 0)
  fp <- sum(pred_presence == 1 & obs_presence == 0)
  fn <- sum(pred_presence == 0 & obs_presence == 1)
  
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  
  # AUC for presence/absence prediction using predicted counts as proxy
  auc <- tryCatch(pROC::auc(obs_presence, predicted, quiet = TRUE), 
                  error = function(e) NA)
  
  # Regression metrics for count data
  resids <- predicted - observed
  rmse <- sqrt(mean(resids^2, na.rm = TRUE))
  mae <- mean(abs(resids), na.rm = TRUE)
  bias <- mean(resids, na.rm = TRUE)
  
  # R-squared for count data
  r2 <- tryCatch({
    cor(observed, predicted, use = "complete.obs")^2
  }, error = function(e) NA)
  
  # Normalized MAE
  nmae <- ifelse(mean(observed, na.rm = TRUE) > 0, 
                 mae / mean(observed, na.rm = TRUE), NA)
  
  # Zero prediction metrics
  obs_zeros <- sum(observed == 0)
  pred_zeros <- sum(predicted == 0)
  zero_accuracy <- mean((observed == 0) == (predicted == 0))
  
  return(list(
    # Classification metrics
    acc = acc, sens = sens, spec = spec, auc = as.numeric(auc),
    # Regression metrics
    rmse = rmse, mae = mae, bias = bias, r2 = r2, nmae = nmae,
    # Zero prediction metrics
    obs_zeros = obs_zeros, pred_zeros = pred_zeros, zero_accuracy = zero_accuracy
  ))
}

cat("Running Quasi-Poisson model with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  set.seed(i)
  
  # Use same train-test splits as RF model for fair comparison
  if (i <= length(training_data_list_final) && i <= length(test_data_list_final)) {
    train <- training_data_list_final[[i]]
    test <- test_data_list_final[[i]]
  } else {
    # Fallback: create new split
    repeat {
      split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      # Check factor levels
      all_ok <- TRUE
      for (fac in factor_vars_presence) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
      if (all_ok && length(unique(train$Presence)) >= 2) break
    }
    
    # Align factor levels
    for (fac in factor_vars_presence) {
      if (fac %in% names(test)) {
        test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
      }
    }
  }
  
  ## Quasi-Poisson Model
  tryCatch({
    # Prepare count data - convert abundance to integer counts
    # You may need to scale/round your abundance data appropriately
    train$Sanionia_count <- round(train$Sanionia)
    test$Sanionia_count <- round(test$Sanionia)
    
    # Ensure non-negative integers
    train$Sanionia_count <- pmax(0, train$Sanionia_count)
    test$Sanionia_count <- pmax(0, test$Sanionia_count)
    
    # Create formula combining both presence and abundance predictors
    all_predictors <- unique(c(predictors_presence, predictors_abundance))
    qp_formula <- as.formula(paste("Sanionia_count ~", paste(all_predictors, collapse = " + ")))
    
    # Fit Quasi-Poisson model
    qp_model <- glm(qp_formula, data = train, family = quasipoisson())
    
    # Check for convergence issues
    if (!qp_model$converged) {
      warning(paste("Quasi-Poisson model did not converge in iteration", i))
    }
    
    # Predict on test set
    pred_counts_qp <- predict(qp_model, newdata = test, type = "response")
    
    # Calculate metrics
    metrics_qp <- calc_qp_metrics(pred_counts_qp, test$Sanionia_count)
    
    # Get dispersion parameter
    dispersion <- summary(qp_model)$dispersion
    
    all_results_qp[[i]] <- data.table(
      Iteration = i,
      # Classification metrics (presence/absence)
      Accuracy = metrics_qp$acc,
      Sensitivity = metrics_qp$sens,
      Specificity = metrics_qp$spec,
      AUC = metrics_qp$auc,
      # Regression metrics (count prediction)
      RMSE = metrics_qp$rmse,
      MAE = metrics_qp$mae,
      Bias = metrics_qp$bias,
      R_squared = metrics_qp$r2,
      nMAE = metrics_qp$nmae,
      # Zero prediction metrics
      Obs_Zeros = metrics_qp$obs_zeros,
      Pred_Zeros = metrics_qp$pred_zeros,
      Zero_Accuracy = metrics_qp$zero_accuracy,
      # Model diagnostics
      Dispersion = dispersion,
      Converged = qp_model$converged,
      N_Test = nrow(test),
      N_Train = nrow(train)
    )
    
    # Store predictions
    train_pred_counts <- predict(qp_model, newdata = train, type = "response")
    
    predictions_all_qp[[i]] <- rbind(
      data.table(set = "test", obs = test$Sanionia_count, pred = pred_counts_qp, 
                 obs_orig = test$Sanionia, iteration = i),
      data.table(set = "train", obs = train$Sanionia_count, pred = train_pred_counts, 
                 obs_orig = train$Sanionia, iteration = i)
    )
    
    # Store model for first iteration
    if (i == 1) {
      model_summaries_qp[["qp_model"]] <- list(
        model = qp_model,
        summary = summary(qp_model),
        aic = AIC(qp_model),
        formula = qp_formula,
        dispersion = dispersion,
        deviance = qp_model$deviance,
        null_deviance = qp_model$null.deviance
      )
    }
    
  }, error = function(e) {
    cat("Quasi-Poisson Error in iteration", i, ":", e$message, "\n")
    
    # Store NA results for failed iterations
    all_results_qp[[i]] <- data.table(
      Iteration = i, Accuracy = NA, Sensitivity = NA, Specificity = NA, AUC = NA,
      RMSE = NA, MAE = NA, Bias = NA, R_squared = NA, nMAE = NA,
      Obs_Zeros = NA, Pred_Zeros = NA, Zero_Accuracy = NA,
      Dispersion = NA, Converged = FALSE, N_Test = NA, N_Train = NA
    )
    
    predictions_all_qp[[i]] <- NULL
  })
  
  if (i %% 10 == 0) cat("Quasi-Poisson iteration", i, "completed\n")
}

# Combine Quasi-Poisson results
results_qp_combined <- rbindlist(all_results_qp)
predictions_qp_combined <- rbindlist(predictions_all_qp[!sapply(predictions_all_qp, is.null)])

#### MODEL COMPARISON ANALYSIS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RANDOM FOREST vs QUASI-POISSON COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# PRESENCE PREDICTION COMPARISON
cat("\n--- PRESENCE PREDICTION COMPARISON ---\n")
cat("Random Forest vs Quasi-Poisson\n\n")

rf_clf_summary <- results_clf_final_combined[, .(
  Model = "Random Forest",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

qp_clf_summary <- results_qp_combined[, .(
  Model = "Quasi-Poisson",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

presence_comparison <- rbind(rf_clf_summary, qp_clf_summary)
print(presence_comparison)

# ABUNDANCE COMPARISON
cat("\n--- ABUNDANCE PREDICTION COMPARISON ---\n")
cat("Quantile Random Forest vs Quasi-Poisson\n\n")

if (nrow(results_qrf_final_combined) > 0 && nrow(results_qp_combined) > 0) {
  rf_abundance_summary <- results_qrf_final_combined[, .(
    Model = "Quantile RF",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  qp_abundance_summary <- results_qp_combined[, .(
    Model = "Quasi-Poisson",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  abundance_comparison <- rbind(rf_abundance_summary, qp_abundance_summary)
  print(abundance_comparison)
  
  # Statistical tests for abundance metrics
  cat("\n--- STATISTICAL SIGNIFICANCE TESTS ---\n")
  
  # Presence metrics
  accuracy_test <- t.test(results_clf_final_combined$Accuracy, results_qp_combined$Accuracy)
  auc_test <- t.test(results_clf_final_combined$AUC, results_qp_combined$AUC)
  
  cat("Accuracy difference (RF - QP):\n")
  cat("  Mean difference:", round(accuracy_test$estimate[1] - accuracy_test$estimate[2], 4), "\n")
  cat("  p-value:", round(accuracy_test$p.value, 4), ifelse(accuracy_test$p.value < 0.05, "*", ""), "\n")
  
  cat("AUC difference (RF - QP):\n")
  cat("  Mean difference:", round(auc_test$estimate[1] - auc_test$estimate[2], 4), "\n")
  cat("  p-value:", round(auc_test$p.value, 4), ifelse(auc_test$p.value < 0.05, "*", ""), "\n")
  
  # Abundance metrics
  rf_rmse_valid <- results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]
  qp_rmse_valid <- results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)]
  
  if (length(rf_rmse_valid) > 10 && length(qp_rmse_valid) > 10) {
    rmse_test <- t.test(rf_rmse_valid, qp_rmse_valid)
    cat("RMSE difference (QRF - QP):\n")
    cat("  Mean difference:", round(rmse_test$estimate[1] - rmse_test$estimate[2], 4), "\n")
    cat("  p-value:", round(rmse_test$p.value, 4), ifelse(rmse_test$p.value < 0.05, "*", ""), "\n")
    
    rf_r2_valid <- results_qrf_final_combined$R_squared[!is.na(results_qrf_final_combined$R_squared)]
    qp_r2_valid <- results_qp_combined$R_squared[!is.na(results_qp_combined$R_squared)]
    
    if (length(rf_r2_valid) > 10 && length(qp_r2_valid) > 10) {  
      r2_test <- t.test(rf_r2_valid, qp_r2_valid)
      cat("R² difference (QRF - QP):\n")
      cat("  Mean difference:", round(r2_test$estimate[1] - r2_test$estimate[2], 4), "\n")
      cat("  p-value:", round(r2_test$p.value, 4), ifelse(r2_test$p.value < 0.05, "*", ""), "\n")
    }
  }
}

#### OVERDISPERSION ANALYSIS ####
cat("\n--- OVERDISPERSION ANALYSIS ---\n")
qp_dispersion_summary <- results_qp_combined[, .(
  Mean_Dispersion = round(mean(Dispersion, na.rm = TRUE), 3),
  SD_Dispersion = round(sd(Dispersion, na.rm = TRUE), 3),
  Min_Dispersion = round(min(Dispersion, na.rm = TRUE), 3),
  Max_Dispersion = round(max(Dispersion, na.rm = TRUE), 3),
  Convergence_Rate = round(mean(Converged, na.rm = TRUE), 3)
)]

cat("Overdispersion analysis:\n")
print(qp_dispersion_summary)

if (qp_dispersion_summary$Mean_Dispersion > 1.5) {
  cat("\nNote: Mean dispersion > 1.5 indicates overdispersion in the data\n")
} else if (qp_dispersion_summary$Mean_Dispersion < 0.8) {
  cat("\nNote: Mean dispersion < 0.8 indicates potential underdispersion\n")
} else {
  cat("\nNote: Dispersion close to 1, suggesting Poisson assumption may be reasonable\n")
}

#### MODEL INTERPRETATION ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== QUASI-POISSON MODEL INTERPRETATION ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

if (!is.null(model_summaries_qp[["qp_model"]])) {
  cat("\n--- QUASI-POISSON COEFFICIENTS ---\n")
  
  qp_summary <- model_summaries_qp[["qp_model"]]$summary
  
  # Model coefficients
  cat("\nQUASI-POISSON MODEL COEFFICIENTS:\n")
  coefs <- qp_summary$coefficients
  coef_df <- data.frame(
    Variable = rownames(coefs),
    Estimate = round(coefs[, "Estimate"], 4),
    Std_Error = round(coefs[, "Std. Error"], 4),
    t_value = round(coefs[, "t value"], 3),
    p_value = round(coefs[, "Pr(>|t|)"], 4),
    Significant = ifelse(coefs[, "Pr(>|t|)"] < 0.05, "*", "")
  )
  print(coef_df)
  
  cat("\nInterpretation: Coefficients represent log rate ratios\n")
  cat("Positive coefficients increase expected count, negative decrease it\n")
  
  # Model fit statistics
  cat("\nMODEL FIT STATISTICS:\n")
  cat("Dispersion parameter:", round(model_summaries_qp[["qp_model"]]$dispersion, 3), "\n")
  cat("Residual deviance:", round(model_summaries_qp[["qp_model"]]$deviance, 2), "\n")
  cat("Null deviance:", round(model_summaries_qp[["qp_model"]]$null_deviance, 2), "\n")
  
  # Pseudo R-squared
  pseudo_r2 <- 1 - (model_summaries_qp[["qp_model"]]$deviance / model_summaries_qp[["qp_model"]]$null_deviance)
  cat("Pseudo R-squared:", round(pseudo_r2, 3), "\n")
  
  cat("AIC:", round(model_summaries_qp[["qp_model"]]$aic, 2), "\n")
}

# RF Variable Importance comparison
cat("\n--- RANDOM FOREST VARIABLE IMPORTANCE ---\n")
cat("Top 10 variables for presence prediction:\n")
print(head(presence_final_importance[, c("Variable", "Importance_pct")], 10))

if (nrow(abundance_final_importance) > 0) {
  cat("\nTop 10 variables for abundance prediction:\n")
  print(head(abundance_final_importance[, c("Variable", "Importance_pct")], 10))
}

#### PREDICTION COMPARISON PLOTS ####
cat("\n=== CREATING COMPARISON PLOTS ===\n")

if (nrow(predictions_qp_combined) > 0 && nrow(predictions_final_combined) > 0) {
  # Combine predictions for plotting (using original scale)
  rf_pred_plot <- predictions_final_combined[set == "test", .(obs, pred, model = "Random Forest")]
  qp_pred_plot <- predictions_qp_combined[set == "test", .(obs = obs_orig, pred, model = "Quasi-Poisson")]
  
  combined_pred_plot <- rbind(rf_pred_plot, qp_pred_plot)
  
  # Scatter plot comparison
  p_comparison <- ggplot(combined_pred_plot, aes(x = obs, y = pred, color = model)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~model) +
    labs(
      title = "Sanionia Abundance: Random Forest vs Quasi-Poisson",
      subtitle = "Test set predictions (perfect predictions would fall on diagonal line)",
      x = "Observed Abundance",
      y = "Predicted Abundance",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Random Forest" = "#EBCC2A", "Quasi-Poisson" = "#F21A00"))
  
  print(p_comparison)
  
  # Performance metrics comparison plot
  metrics_comparison <- rbind(
    data.table(Model = "Random Forest", 
               RMSE = results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]),
    data.table(Model = "Quasi-Poisson", 
               RMSE = results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)])
  )
  
  if (nrow(metrics_comparison) > 0) {
    p_metrics <- ggplot(metrics_comparison, aes(x = Model, y = RMSE, fill = Model)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = "RMSE Distribution: Random Forest vs Quasi-Poisson",
        subtitle = "Lower RMSE indicates better performance",
        y = "Root Mean Square Error"
      ) +
      theme_minimal() +
      scale_fill_manual(values = c("Random Forest" = "#EBCC2A", "Quasi-Poisson" = "#F21A00"))
    
    print(p_metrics)
  }
  
  # Residuals vs fitted plot for Quasi-Poisson
  if (nrow(predictions_qp_combined[set == "test"]) > 0) {
    qp_test_data <- predictions_qp_combined[set == "test"]
    qp_test_data[, residuals := obs_orig - pred]
    
    p_residuals <- ggplot(qp_test_data, aes(x = pred, y = residuals)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue") +
      labs(
        title = "Quasi-Poisson Model: Residuals vs Fitted Values",
        subtitle = "Test set predictions - should show random scatter around zero",
        x = "Fitted Values",
        y = "Residuals"
      ) +
      theme_minimal()
    
    print(p_residuals)
  }
}

#### SUMMARY AND RECOMMENDATIONS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== SUMMARY AND RECOMMENDATIONS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\nMODEL COMPARISON:\n")
cat("- Random Forest: Non-parametric, handles complex interactions, ensemble method\n")
cat("- Quasi-Poisson: Parametric, accounts for overdispersion, interpretable coefficients\n")

cat("\nQUASI-POISSON ADVANTAGES:\n")
cat("- Handles overdispersion in count data (dispersion parameter ≠ 1)\n")
cat("- More flexible than standard Poisson for overdispersed data\n")
cat("- Provides interpretable coefficients (log rate ratios)\n")
cat("- Well-established statistical foundation for count data\n")
cat("- Uses same mean-variance relationship as Poisson but allows different variance\n")

cat("\nPERFORMANCE SUMMARY:\n")
if (exists("presence_comparison")) {
  rf_better_acc <- presence_comparison$Mean_Accuracy[1] > presence_comparison$Mean_Accuracy[2]
  rf_better_auc <- presence_comparison$Mean_AUC[1] > presence_comparison$Mean_AUC[2]
  
  cat("Presence prediction:\n")
  cat("  - Better accuracy:", ifelse(rf_better_acc, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Better AUC:", ifelse(rf_better_auc, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("abundance_comparison") && nrow(abundance_comparison) > 0) {
  rf_better_rmse <- abundance_comparison$Mean_RMSE[1] < abundance_comparison$Mean_RMSE[2]
  rf_better_r2 <- abundance_comparison$Mean_R2[1] > abundance_comparison$Mean_R2[2]
  
  cat("Count prediction:\n")
  cat("  - Lower RMSE (better):", ifelse(rf_better_rmse, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Higher R² (better):", ifelse(rf_better_r2, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("qp_dispersion_summary")) {
  cat("Overdispersion:\n")
  cat("  - Mean dispersion parameter:", qp_dispersion_summary$Mean_Dispersion, "\n")
  if (qp_dispersion_summary$Mean_Dispersion > 1.2) {
    cat("  - Data shows overdispersion, Quasi-Poisson is appropriate\n")
  }
}

cat("\nRECOMMENDATION:\n")
cat("Choose Random Forest if: Complex non-linear relationships, prediction accuracy priority\n")
cat("Choose Quasi-Poisson if: Need interpretable coefficients, understand linear relationships\n")
cat("Consider Quasi-Poisson if: Count data with overdispersion, parametric approach preferred\n")






######################################## DESCHAMPSIA ###########################################################################################

rm(list = ls())

#### Deschampsia ####
### DESCHAMPSIA TRANSFORMATIONS FOR HURDLE MODEL 
library(randomForest)
library(quantregForest)
library(data.table)
library(dplyr)
library(rsample)
library(isotone)
library(pROC)
library(pdp)
library(clipr)
library(ggplot2)

select <- dplyr::select

# Load data
df <- read.csv("Deschampsia.csv")
df <- df %>% mutate(across(where(is.character), as.factor))
df$Presence <- as.factor(ifelse(df$Deschampsia > 0, 1, 0))

# Define transformation functions
transformations <- list(
  none = list(
    transform = function(x) x,
    inverse = function(x) x,
    name = "None"
  ),
  sqrt = list(
    transform = function(x) sqrt(pmax(x, 0)),
    inverse = function(x) pmax(x^2, 0),
    name = "Square Root"
  ),
  log = list(
    transform = function(x) log(pmax(x, 0.001)),  # Add small constant to avoid log(0)
    inverse = function(x) exp(x),
    name = "Log"
  ),
  asinh = list(
    transform = function(x) asinh(x),
    inverse = function(x) sinh(x),
    name = "Asinh"
  )
)

predictors <- c("PDIS_CLASS", "ASP", "SLOPE_CLASS", "MORF_CLASS", "CD_CLASS", "DISTCOAST", "EFAH",
                "VDEP", "TAS_08", "WEXP", "Bio9_81", "WEFF", "GDGFGD0", "PR_01", "HillSh", "FCF",
                "RSP", "SVF", "LSF", "MPI", "WSI", "Bio3_81", "MNCURV", "CSCURV", "VDISCN",
                "CVI", "FLAC", "PNCURV", "TXT", "TPW", "MBI", "Easting", "Northing")
factor_vars <- names(df[, predictors])[sapply(df[, predictors], is.factor)]

n_iter <- 100

# Initialize storage for all transformations
all_results_clf <- list()
all_results_qrf <- list()
varimp_list_clf <- list()
varimp_list_qrf <- list()
predictions_all <- list()
training_data_list <- list()
test_data_list <- list()

for (trans_name in names(transformations)) {
  trans_func <- transformations[[trans_name]]
  
  cat("Testing transformation:", trans_func$name, "\n")
  
  # Transform the response variable
  df_trans <- df
  df_trans$Deschampsia_trans <- trans_func$transform(df$Deschampsia)
  
  all_results_clf[[trans_name]] <- list()
  all_results_qrf[[trans_name]] <- list()
  varimp_list_clf[[trans_name]] <- list()
  varimp_list_qrf[[trans_name]] <- list()
  predictions_all[[trans_name]] <- list()
  training_data_list[[trans_name]] <- list()
  test_data_list[[trans_name]] <- list()
  
  for (i in 1:n_iter) {
    set.seed(i)
    repeat {
      split <- initial_split(df_trans, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      all_ok <- TRUE
      for (fac in factor_vars) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) {
          all_ok <- FALSE
          break
        }
      }
      if (all_ok) break
    }
    
    for (fac in factor_vars) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
    
    ## STEP 1: RF Classifier (same for all transformations)
    rf_clf <- randomForest(
      x = train[, predictors],
      y = train$Presence,
      ntree = 500,
      importance = TRUE,
      keep.forest = TRUE
    )
    
    pred_probs <- predict(rf_clf, newdata = test[, predictors], type = "prob")[, "1"]
    pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
    obs_class <- test$Presence
    
    # Metrics for RF classifier
    acc <- mean(pred_class == obs_class)
    sens <- mean(pred_class == "1" & obs_class == "1")
    spec <- mean(pred_class == "0" & obs_class == "0")
    logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
    brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
    auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
    
    spread_clf <- abs(0.5 - pred_probs)  # proxy for uncertainty
    
    all_results_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
      LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
      MedianSpread = median(spread_clf)
    )
    
    imp_clf <- importance(rf_clf)
    varimp_list_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Variable = rownames(imp_clf), 
      Importance = imp_clf[, 1], 
      Iteration = i, 
      Model = "RF_Presence"
    )
    
    ## STEP 2: QRF for predicted present cases (with transformation)
    test_presence_idx <- which(pred_class == "1")
    test_present_data <- test[test_presence_idx, ]
    train_present_data <- train[train$Presence == "1", ]
    
    if (nrow(test_present_data) > 10) {
      qrf <- quantregForest(
        x = train_present_data[, predictors],
        y = train_present_data$Deschampsia_trans,  # Use transformed response
        ntree = 500
      )
      
      preds_test_trans <- predict(qrf, newdata = test_present_data[, predictors], what = c(0.025, 0.5, 0.975))
      
      # Apply isotonic regression on transformed scale
      iso_model <- isoreg(preds_test_trans[, 2], test_present_data$Deschampsia_trans)
      iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
      calibrated_test_trans <- iso_fun(preds_test_trans[, 2])
      
      # Back-transform predictions and bounds
      preds_test_lower <- trans_func$inverse(preds_test_trans[, 1])
      preds_test_median <- trans_func$inverse(calibrated_test_trans)
      preds_test_upper <- trans_func$inverse(preds_test_trans[, 3])
      
      # Ensure predictions are within reasonable bounds
      preds_test_lower <- pmax(pmin(preds_test_lower, 100), 0)
      preds_test_median <- pmax(pmin(preds_test_median, 100), 0)
      preds_test_upper <- pmax(pmin(preds_test_upper, 100), 0)
      
      spread_test <- preds_test_upper - preds_test_lower
      resids_test <- preds_test_median - test_present_data$Deschampsia
      
      rmse <- sqrt(mean(resids_test^2))
      mae <- mean(abs(resids_test))
      bias <- mean(resids_test)
      r2 <- cor(test_present_data$Deschampsia, preds_test_median)^2
      nmae <- mae / mean(test_present_data$Deschampsia)
      coverage <- mean(test_present_data$Deschampsia >= preds_test_lower & test_present_data$Deschampsia <= preds_test_upper)
      median_spread <- median(spread_test)
      low_bias <- mean(preds_test_lower - test_present_data$Deschampsia)
      high_bias <- mean(preds_test_upper - test_present_data$Deschampsia)
      
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
        R_squared = r2, nMAE = nmae, Coverage_95PI = coverage,
        Median_Spread_95PI = median_spread,
        Lower_Q2.5_Bias = low_bias,
        Upper_Q97.5_Bias = high_bias
      )
      
      # Store predictions (back-transformed for training data too)
      train_preds_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.5)
      train_spread_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.975) -
        predict(qrf, newdata = train_present_data[, predictors], what = 0.025)
      
      predictions_all[[trans_name]][[i]] <- rbind(
        data.table(set = "test", obs = test_present_data$Deschampsia, pred = preds_test_median,
                   spread = spread_test, iteration = i, transformation = trans_func$name),
        data.table(set = "train", obs = train_present_data$Deschampsia,
                   pred = trans_func$inverse(train_preds_trans),
                   spread = trans_func$inverse(train_present_data$Deschampsia_trans + train_spread_trans/2) - 
                     trans_func$inverse(train_present_data$Deschampsia_trans - train_spread_trans/2),
                   iteration = i, transformation = trans_func$name)
      )
      
      imp_qrf <- importance(qrf)
      varimp_list_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Variable = rownames(imp_qrf), 
        Importance = imp_qrf[, 1], 
        Iteration = i, 
        Model = "QRF_Abundance"
      )
    } else {
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
        R_squared = NA, nMAE = NA, Coverage_95PI = NA,
        Median_Spread_95PI = NA,
        Lower_Q2.5_Bias = NA,
        Upper_Q97.5_Bias = NA
      )
      predictions_all[[trans_name]][[i]] <- NULL
      varimp_list_qrf[[trans_name]][[i]] <- NULL
    }
    
    training_data_list[[trans_name]][[i]] <- train
    test_data_list[[trans_name]][[i]] <- test
    
    if (i %% 10 == 0) cat("Transformation:", trans_func$name, "- Iteration", i, "\n")
  }
}


# Combine results across transformations
results_clf_combined <- rbindlist(lapply(all_results_clf, function(x) rbindlist(x)))
results_qrf_combined <- rbindlist(lapply(all_results_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_clf_combined <- rbindlist(lapply(varimp_list_clf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_qrf_combined <- rbindlist(lapply(varimp_list_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
predictions_combined <- rbindlist(lapply(predictions_all, function(x) rbindlist(x[!sapply(x, is.null)])))

# Summary statistics by transformation
cat("\n=== CLASSIFICATION RESULTS SUMMARY ===\n")
clf_summary <- results_clf_combined[, .(
  Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
  SD_Accuracy = sd(Accuracy, na.rm = TRUE),
  Mean_AUC = mean(AUC, na.rm = TRUE),
  SD_AUC = sd(AUC, na.rm = TRUE),
  Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
  Mean_Specificity = mean(Specificity, na.rm = TRUE)
), by = Transformation]
print(clf_summary)

cat("\n=== QUANTILE REGRESSION RESULTS SUMMARY ===\n")
qrf_summary <- results_qrf_combined[, .(
  Mean_RMSE = mean(RMSE, na.rm = TRUE),
  SD_RMSE = sd(RMSE, na.rm = TRUE),
  Mean_MAE = mean(MAE, na.rm = TRUE),
  SD_MAE = sd(MAE, na.rm = TRUE),
  Mean_R2 = mean(R_squared, na.rm = TRUE),
  SD_R2 = sd(R_squared, na.rm = TRUE),
  Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
  Mean_Bias = mean(Bias, na.rm = TRUE)
), by = Transformation]
print(qrf_summary)

# 1. PREPARE DATA FOR PLOTTING

# Spread plot data (prediction interval width vs observed)
spread_plot_data <- predictions_combined[set == "test", .(
  obs = obs,
  spread = spread,
  trans = transformation,
  iteration = iteration
)]

# Observed vs predicted data
pred_obs_data <- predictions_combined[set == "test", .(
  obs = obs,
  pred = pred,
  trans = transformation,
  iteration = iteration
)]

# Calculate residuals for each transformation
residuals_data <- predictions_combined[set == "test", .(
  pred = pred,
  residuals = obs - pred,
  trans = transformation,
  iteration = iteration
)]

# 2. SPREAD PLOT (Uncertainty vs Observed)
p1 <- ggplot(spread_plot_data, aes(x = obs, y = spread, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Observed Deschampsia Abundance", y = "Prediction Interval Width",
       title = "Uncertainty vs Observed (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# 3. VARIABLE IMPORTANCE ANALYSIS
cat("\n=== VARIABLE IMPORTANCE ANALYSIS ===\n")

# Calculate variable importance summary
varimp_summary <- varimp_qrf_combined %>%
  group_by(Transformation, Variable) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 4),
            .groups = "drop") %>%
  arrange(Transformation, desc(MedianImportance))

# Calculate percentage importance
varimp_summary <- varimp_summary %>%
  group_by(Transformation) %>%
  mutate(Importance_pct = 100 * MedianImportance / sum(MedianImportance)) %>%
  ungroup()

# For "None" transformation (equivalent to "raw")
none_imp <- varimp_summary %>%
  filter(Transformation == "None") %>%
  mutate(Cumulative_Importance = round(cumsum(Importance_pct), 1)) %>%
  arrange(desc(MedianImportance))

cat("Variable Importance for No Transformation (Raw):\n")
print(none_imp)
clipr::write_clip(none_imp)

# Variable importance comparison plot
p2 <- varimp_summary %>%
  group_by(Transformation) %>%
  slice_head(n = 10) %>%  # Top 10 variables per transformation
  ggplot(aes(x = reorder(Variable, MedianImportance), y = MedianImportance, fill = Transformation)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Variable", y = "Median Importance", 
       title = "Top 10 Variable Importance by Transformation") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  facet_wrap(~Transformation, scales = "free_y")

print(p2)

# 4. OBSERVED VS PREDICTED PLOT
p3 <- ggplot(pred_obs_data, aes(x = obs, y = pred, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed Deschampsia Abundance", y = "Predicted Deschampsia Abundance",
       title = "Observed vs Predicted (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# 5. RESIDUALS VS PREDICTED (for each transformation)
p4 <- ggplot(residuals_data, aes(x = pred, y = residuals)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~trans, scales = "free") +
  theme_minimal() +
  labs(title = "Residuals vs Predicted (by Transformation)",
       x = "Predicted Deschampsia Abundance",
       y = "Residuals (Observed - Predicted)")

print(p4)

# 6. INDIVIDUAL TRANSFORMATION PLOTS

# Function to create plots for a specific transformation
create_transformation_plots <- function(trans_name) {
  
  # Filter data for specific transformation
  trans_spread <- spread_plot_data[trans == trans_name]
  trans_pred_obs <- pred_obs_data[trans == trans_name]
  trans_residuals <- residuals_data[trans == trans_name]
  
  cat(paste("\n=== PLOTS FOR", trans_name, "TRANSFORMATION ===\n"))
  
  # Residuals vs Predicted for specific transformation
  p_resid <- ggplot(trans_residuals, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Residuals vs Predicted -", trans_name, "Transformation"),
         x = "Predicted Deschampsia Abundance",
         y = "Residuals (Observed - Predicted)")
  
  # Prediction Interval Width vs Observed for specific transformation
  p_spread <- ggplot(trans_spread, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    theme_minimal() +
    labs(title = paste("Prediction Interval Width vs Observed -", trans_name, "Transformation"),
         x = "Observed Deschampsia Abundance",
         y = "95% Prediction Interval Width")
  
  print(p_resid)
  print(p_spread)
  
  return(list(residuals = p_resid, spread = p_spread))
}

# Create individual plots for each transformation
transformation_plots <- list()
for (trans in unique(spread_plot_data$trans)) {
  transformation_plots[[trans]] <- create_transformation_plots(trans)
}

# 7. PERFORMANCE SUMMARY COMPARISON PLOT
perf_summary_long <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI, Bias) %>%
  group_by(Transformation) %>%
  summarize(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Mean_MAE = mean(MAE, na.rm = TRUE),
    Mean_R2 = mean(R_squared, na.rm = TRUE),
    Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
    Mean_Bias = abs(mean(Bias, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p5 <- ggplot(perf_summary_long, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_col() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Model Performance Comparison Across Transformations",
       x = "Transformation", y = "Mean Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p5)

# 8. BOXPLOT OF PERFORMANCE METRICS
perf_boxplot_data <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p6 <- ggplot(perf_boxplot_data, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Distribution of Performance Metrics Across Iterations",
       x = "Transformation", y = "Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p6)

# 9. CORRELATION BETWEEN OBSERVED AND PREDICTED BY TRANSFORMATION
cor_summary <- pred_obs_data %>%
  group_by(trans) %>%
  summarize(
    correlation = cor(obs, pred, use = "complete.obs"),
    rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    mae = mean(abs(obs - pred), na.rm = TRUE),
    .groups = "drop")
print(cor_summary)

#### EXTRACT TOP PREDICTORS FROM TRANSFORMATION ANALYSIS ####

# Extract top predictors from presence model (RF Classifier) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_presence <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

# Extract top predictors from abundance model (QRF) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_abundance <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

cat("\n=== TOP PREDICTORS EXTRACTED FROM RAW TRANSFORMATION (70% IMPORTANCE) ===\n")
cat("Top predictors for PRESENCE model (explaining 70% of importance):\n")
print(top_predictors_presence)
cat("Number of presence predictors:", length(top_predictors_presence), "\n")

cat("\nTop predictors for ABUNDANCE model (explaining 70% of importance):\n")
print(top_predictors_abundance)
cat("Number of abundance predictors:", length(top_predictors_abundance), "\n")

# Optional: Show the importance breakdown
cat("\n=== IMPORTANCE BREAKDOWN FOR SELECTED PREDICTORS ===\n")

# Presence model importance breakdown
presence_importance_breakdown <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_presence
  ) %>%
  filter(Selected)

cat("PRESENCE model - Selected variables and their importance:\n")
print(presence_importance_breakdown)

# Abundance model importance breakdown
abundance_importance_breakdown <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_abundance
  ) %>%
  filter(Selected)

cat("\nABUNDANCE model - Selected variables and their importance:\n")
print(abundance_importance_breakdown)


#### HURDLE FINAL MODEL USING TOP PREDICTORS ####

# Use the extracted top predictors for each model
predictors_presence <- c("DISTCOAST", "WEXP", "VDEP", "EFAH", "Easting", "MPI", "Northing", "WEFF", "LSF", "ASP", "HillSh", "PR_01", "MBI", "PDIS_CLASS", "RSP", "FLAC", "TAS_08")
predictors_abundance <- c("TPW", "LSF", "Northing", "MPI", "Easting", "ASP", "DISTCOAST", "EFAH", "MBI", "WEXP", "MNCURV", "FLAC", "SVF")

# Get factor variables for each predictor set
factor_vars_presence <- names(df[, predictors_presence])[sapply(df[, predictors_presence], is.factor)]
factor_vars_abundance <- names(df[, predictors_abundance])[sapply(df[, predictors_abundance], is.factor)]

n_iter <- 100
all_results_clf_final <- list()
all_results_qrf_final <- list()
varimp_list_clf_final <- list()
varimp_list_qrf_final <- list()
predictions_all_final <- list()
training_data_list_final <- list()
test_data_list_final <- list()

cat("\n=== RUNNING FINAL MODEL WITH TOP PREDICTORS ===\n")

for (i in 1:n_iter) {
  set.seed(i)
  repeat {
    split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
    train <- training(split)
    test <- testing(split)
    
    # Check factor levels for BOTH predictor sets
    all_ok <- TRUE
    for (fac in factor_vars_presence) {
      train_levels <- levels(train[[fac]])
      test_levels <- levels(factor(test[[fac]]))
      if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
    }
    if (all_ok) {
      for (fac in factor_vars_abundance) {
        train_levels <- levels(train[[fac]])
        test_levels <- levels(factor(test[[fac]]))
        if (!all(train_levels %in% test_levels)) { all_ok <- FALSE; break }
      }
    }
    if (all_ok) break
  }
  
  # Set factor levels
  for (fac in factor_vars_presence) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  for (fac in factor_vars_abundance) test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
  
  ## STEP 1: RF Classifier
  rf_clf <- randomForest(
    x = train[, predictors_presence],
    y = train$Presence,
    ntree = 500,
    importance = TRUE,
    keep.forest = TRUE
  )
  
  pred_probs <- predict(rf_clf, newdata = test[, predictors_presence], type = "prob")[, "1"]
  pred_class <- as.factor(ifelse(pred_probs >= 0.5, 1, 0))
  obs_class <- test$Presence
  
  # Metrics for classifier
  acc <- mean(pred_class == obs_class)
  sens <- mean(pred_class == "1" & obs_class == "1")
  spec <- mean(pred_class == "0" & obs_class == "0")
  logloss <- -mean(ifelse(obs_class == "1", log(pred_probs), log(1 - pred_probs)))
  brier <- mean((as.numeric(as.character(obs_class)) - pred_probs)^2)
  auc <- tryCatch(pROC::auc(obs_class, pred_probs), error = function(e) NA)
  spread_clf <- abs(0.5 - pred_probs)
  
  all_results_clf_final[[i]] <- data.table(
    Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
    LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
    MedianSpread = median(spread_clf)
  )
  
  varimp_list_clf_final[[i]] <- data.table(
    Variable = rownames(importance(rf_clf)),
    Importance = importance(rf_clf)[,1],
    Iteration = i, Model = "RF_Presence"
  )
  
  ## STEP 2: QRF for abundance
  test_presence_idx <- which(pred_class == "1")
  test_present_data <- test[test_presence_idx, ]
  train_present_data <- train[train$Presence == "1", ]
  
  if (nrow(test_present_data) > 10) {
    qrf <- quantregForest(
      x = train_present_data[, predictors_abundance],
      y = train_present_data$Deschampsia,
      ntree = 500
    )
    
    preds_test <- predict(qrf, newdata = test_present_data[, predictors_abundance], what = c(0.025, 0.5, 0.975))
    iso_model <- isoreg(preds_test[, 2], test_present_data$Deschampsia)
    iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
    calibrated_test <- pmax(pmin(iso_fun(preds_test[, 2]), 100), 0)
    
    spread_test <- preds_test[, 3] - preds_test[, 1]
    resids_test <- calibrated_test - test_present_data$Deschampsia
    
    rmse <- sqrt(mean(resids_test^2))
    mae <- mean(abs(resids_test))
    bias <- mean(resids_test)
    r2 <- cor(test_present_data$Deschampsia, calibrated_test)^2
    nrmse <- rmse / mean(test_present_data$Deschampsia, na.rm = TRUE)
    nrmse_median <- rmse / median(test_present_data$Deschampsia, na.rm = TRUE)
    nmae <- mae / mean(test_present_data$Deschampsia, na.rm = TRUE)
    nmae_median <- mae / median(test_present_data$Deschampsia, na.rm = TRUE)
    coverage <- mean(test_present_data$Deschampsia >= preds_test[,1] & test_present_data$Deschampsia <= preds_test[,3])
    median_spread <- median(spread_test)
    low_bias <- mean(preds_test[,1] - test_present_data$Deschampsia)
    high_bias <- mean(preds_test[,3] - test_present_data$Deschampsia)
    
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
      R_squared = r2,
      nRMSE = nrmse,
      nRMSE_Median = nrmse_median,
      nMAE = nmae,
      nMAE_Median = nmae_median,
      Coverage_95PI = coverage,
      Median_Spread_95PI = median_spread,
      Lower_Q2.5_Bias = low_bias,
      Upper_Q97.5_Bias = high_bias
    )
    
    predictions_all_final[[i]] <- rbind(
      data.table(set = "test", obs = test_present_data$Deschampsia, pred = calibrated_test,
                 spread = spread_test, iteration = i),
      data.table(set = "train", obs = train_present_data$Deschampsia,
                 pred = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.5),
                 spread = predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.975) -
                   predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.025),
                 iteration = i)
    )
    
    varimp_list_qrf_final[[i]] <- data.table(
      Variable = rownames(importance(qrf)),
      Importance = importance(qrf)[,1],
      Iteration = i, Model = "QRF_Abundance"
    )
  } else {
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
      R_squared = NA, nRMSE = NA, nRMSE_Median = NA, nMAE = NA, nMAE_Median = NA,
      Coverage_95PI = NA, Median_Spread_95PI = NA, Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA
    )
    predictions_all_final[[i]] <- NULL
    varimp_list_qrf_final[[i]] <- NULL
  }
  
  training_data_list_final[[i]] <- train
  test_data_list_final[[i]] <- test
  
  if (i %% 10 == 0) cat("Final model iteration", i, "\n")
}

# Combine results
results_clf_final_combined <- rbindlist(all_results_clf_final)
results_qrf_final_combined <- rbindlist(all_results_qrf_final[!sapply(all_results_qrf_final, is.null)])
varimp_clf_final_combined <- rbindlist(varimp_list_clf_final[!sapply(varimp_list_clf_final, is.null)])
varimp_qrf_final_combined <- rbindlist(varimp_list_qrf_final[!sapply(varimp_list_qrf_final, is.null)])
predictions_final_combined <- rbindlist(predictions_all_final[!sapply(predictions_all_final, is.null)])

# CLASSIFICATION SUMMARY
clf_final_summary <- results_clf_final_combined[, .(
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 3),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 3),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 3),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 3),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 3),
  SD_Sensitivity = round(sd(Sensitivity, na.rm = TRUE), 3),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 3),
  SD_Specificity = round(sd(Specificity, na.rm = TRUE), 3),
  Mean_LogLoss = round(mean(LogLoss, na.rm = TRUE), 3),
  SD_LogLoss = round(sd(LogLoss, na.rm = TRUE), 3),
  Mean_BrierScore = round(mean(BrierScore, na.rm = TRUE), 3),
  SD_BrierScore = round(sd(BrierScore, na.rm = TRUE), 3),
  Mean_MedianSpread = round(mean(MedianSpread, na.rm = TRUE), 3),
  SD_MedianSpread = round(sd(MedianSpread, na.rm = TRUE), 3)
)]
print(clf_final_summary)

qrf_final_summary <- results_qrf_final_combined[, .(
  Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
  SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
  Median_RMSE = round(median(RMSE, na.rm = TRUE), 3),
  Mean_MAE = round(mean(MAE, na.rm = TRUE), 3),
  SD_MAE = round(sd(MAE, na.rm = TRUE), 3),
  Median_MAE = round(median(MAE, na.rm = TRUE), 3),
  Mean_R2 = round(mean(R_squared, na.rm = TRUE), 3),
  SD_R2 = round(sd(R_squared, na.rm = TRUE), 3),
  Median_R2 = round(median(R_squared, na.rm = TRUE), 3),
  Median_nRMSE = round(median(nRMSE_Median, na.rm = TRUE), 3),
  IQR_nRMSE   = round(IQR(nRMSE_Median, na.rm = TRUE), 3),
  Mean_Bias = round(mean(Bias, na.rm = TRUE), 3),
  SD_Bias = round(sd(Bias, na.rm = TRUE), 3),
  Median_Bias = round(median(Bias, na.rm = TRUE), 3),
  Mean_Coverage_95PI = round(mean(Coverage_95PI, na.rm = TRUE), 3),
  SD_Coverage_95PI = round(sd(Coverage_95PI, na.rm = TRUE), 3),
  Median_Coverage_95PI = round(median(Coverage_95PI, na.rm = TRUE), 3),
  Mean_Median_Spread_95PI = round(mean(Median_Spread_95PI, na.rm = TRUE), 3),
  SD_Median_Spread_95PI = round(sd(Median_Spread_95PI, na.rm = TRUE), 3),
  Median_Spread_95PI = round(median(Median_Spread_95PI, na.rm = TRUE), 3),
  Mean_Lower_Q2.5_Bias = round(mean(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  SD_Lower_Q2.5_Bias = round(sd(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Median_Lower_Q2.5_Bias = round(median(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Mean_Upper_Q97.5_Bias = round(mean(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  SD_Upper_Q97.5_Bias = round(sd(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  Median_Upper_Q97.5_Bias = round(median(Upper_Q97.5_Bias, na.rm = TRUE), 3)
)]
print(qrf_final_summary)

cat("\nQUANTILE REGRESSION RESULTS (using top abundance predictors):\n")
qrf_final_summary <- results_qrf_final_combined[, .(
  Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
  SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
  Mean_MAE = round(mean(MAE, na.rm = TRUE), 3),
  SD_MAE = round(sd(MAE, na.rm = TRUE), 3),
  Mean_R2 = round(mean(R_squared, na.rm = TRUE), 3),
  SD_R2 = round(sd(R_squared, na.rm = TRUE), 3),
  Mean_nMAE = round(mean(nMAE, na.rm = TRUE), 3),
  SD_nMAE = round(sd(nMAE, na.rm = TRUE), 3),
  Mean_Bias = round(mean(Bias, na.rm = TRUE), 3),
  SD_Bias = round(sd(Bias, na.rm = TRUE), 3),
  Mean_Coverage_95PI = round(mean(Coverage_95PI, na.rm = TRUE), 3),
  SD_Coverage_95PI = round(sd(Coverage_95PI, na.rm = TRUE), 3),
  Mean_Median_Spread_95PI = round(mean(Median_Spread_95PI, na.rm = TRUE), 3),
  SD_Median_Spread_95PI = round(sd(Median_Spread_95PI, na.rm = TRUE), 3),
  Mean_Lower_Q2.5_Bias = round(mean(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  SD_Lower_Q2.5_Bias = round(sd(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Mean_Upper_Q97.5_Bias = round(mean(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  SD_Upper_Q97.5_Bias = round(sd(Upper_Q97.5_Bias, na.rm = TRUE), 3)
)]
print(qrf_final_summary)

#### VARIABLE IMPORTANCE ANALYSIS FOR FINAL MODELS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== VARIABLE IMPORTANCE FROM FINAL MODELS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create comprehensive importance tables from the final model runs
cat("\nTABLE 1: PRESENCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

presence_final_importance <- varimp_clf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(presence_final_importance)

cat("\nPRESENCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(presence_final_importance), "\n"))
cat(paste("- Top variable (", presence_final_importance$Variable[1], ") explains:", 
          presence_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", presence_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 2: ABUNDANCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

abundance_final_importance <- varimp_qrf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(abundance_final_importance)

cat("\nABUNDANCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(abundance_final_importance), "\n"))
cat(paste("- Top variable (", abundance_final_importance$Variable[1], ") explains:", 
          abundance_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", abundance_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

# Create a comparison table showing which variables are most important in each model
cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 3: COMPARISON - Top Variables in Each Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

# Get top 5 from each model for comparison
top_presence_vars <- presence_final_importance$Variable[1:min(5, nrow(presence_final_importance))]
top_abundance_vars <- abundance_final_importance$Variable[1:min(5, nrow(abundance_final_importance))]

# Create comparison dataframe
comparison_table <- data.frame(
  Rank = 1:5,
  Presence_Model = c(top_presence_vars, rep("", 5 - length(top_presence_vars)))[1:5],
  Presence_Importance = c(presence_final_importance$Importance_pct[1:min(5, nrow(presence_final_importance))], 
                          rep(NA, 5 - min(5, nrow(presence_final_importance))))[1:5],
  Abundance_Model = c(top_abundance_vars, rep("", 5 - length(top_abundance_vars)))[1:5],
  Abundance_Importance = c(abundance_final_importance$Importance_pct[1:min(5, nrow(abundance_final_importance))], 
                           rep(NA, 5 - min(5, nrow(abundance_final_importance))))[1:5]
)

print(comparison_table)

# Identify shared important variables
shared_vars <- intersect(top_presence_vars, top_abundance_vars)
if(length(shared_vars) > 0) {
  cat("\nShared important variables between models:", paste(shared_vars, collapse = ", "), "\n")
} else {
  cat("\nNo shared variables in top 5 between models\n")
}

#### PLOTTING FUNCTIONS - RUN THIS FIRST ####

# Define color palette for species
cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00")

library(ggplot2)
library(dplyr)

# Function to extract 70% importance variables from FINAL MODEL results
extract_final_70_percent <- function(varimp_final_combined, species_name, model_type) {
  
  # Calculate importance from final model runs and apply 70% threshold
  importance_data <- varimp_final_combined %>%
    group_by(Variable) %>%
    summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MedianImportance)) %>%
    mutate(
      Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
      Cumulative_pct = round(cumsum(Importance_pct), 2),
      Selected_70pct = Cumulative_pct <= 70,
      Rank = row_number(),
      Species = species_name,
      Model = model_type
    ) %>%
    filter(Selected_70pct) %>%  # Only keep variables that explain ≤70% cumulative importance
    select(Species, Model, Variable, Rank, Importance_pct, Cumulative_pct)
  
  return(importance_data)
}

create_individual_plot <- function(species_name, model_type, data_all) {
  
  data_subset <- data_all %>%
    filter(Species == species_name, Model == model_type) %>%
    arrange(desc(Rank))  # Most important at top
  
  if(nrow(data_subset) == 0) {
    cat(paste("No data found for", species_name, model_type, "\n"))
    return(NULL)
  }
  
  p <- ggplot(data_subset, aes(x = reorder(Variable, Rank), y = Importance_pct)) +
    geom_col(fill = cores[species_name], alpha = 0.8, color = "white", size = 0.5) +
    geom_text(aes(label = paste0(round(Importance_pct, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    labs(
      title = paste(species_name, "-", model_type, "Final Model (70% Variables)"),
      subtitle = paste("Variables explaining", round(max(data_subset$Cumulative_pct), 1), "% importance in FINAL model"),
      x = "Variables (Ranked by Final Model Importance)",
      y = "Individual Importance (%)",
      caption = paste("n =", nrow(data_subset), "variables from final model")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
      plot.subtitle = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

create_final_model_plots <- function(species_name, varimp_clf_final, varimp_qrf_final) {
  
  cat(paste("\n=== Creating plots for", species_name, "final model ===\n"))
  
  # Extract 70% variables from final models
  presence_data <- extract_final_70_percent(varimp_clf_final, species_name, "Presence")
  abundance_data <- extract_final_70_percent(varimp_qrf_final, species_name, "Abundance")
  
  # Combine data
  combined_data <- rbind(presence_data, abundance_data)
  
  # Create individual plots
  presence_plot <- create_individual_plot(species_name, "Presence", combined_data)
  abundance_plot <- create_individual_plot(species_name, "Abundance", combined_data)
  
  # Print plots
  if(!is.null(presence_plot)) {
    print(presence_plot)
    cat(paste("Created", species_name, "presence plot with", nrow(presence_data), "variables\n"))
  }
  
  if(!is.null(abundance_plot)) {
    print(abundance_plot) 
    cat(paste("Created", species_name, "abundance plot with", nrow(abundance_data), "variables\n"))
  }
  
  # Create comparison plot
  if(nrow(combined_data) > 0) {
    comparison_plot <- ggplot(combined_data, aes(x = reorder(Variable, Rank), y = Importance_pct, fill = Model)) +
      geom_col(position = "dodge", alpha = 0.8, color = "white", size = 0.5) +
      coord_flip() +
      labs(
        title = paste(species_name, "- Final Model: Presence vs Abundance (70% Variables)"),
        subtitle = "Variable importance from final models (70% threshold applied)",
        x = "Variables",
        y = "Individual Importance (%)",
        fill = "Model Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
        plot.subtitle = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      ) +
      scale_fill_manual(values = c("Presence" = cores[species_name], "Abundance" = adjustcolor(cores[species_name], alpha.f = 0.6))) +
      facet_wrap(~Model, scales = "free_y", ncol = 2)
    
    print(comparison_plot)
    cat(paste("Created", species_name, "comparison plot\n"))
  }
  
  return(combined_data)
}
#### INSERT THE PLOTTING CODE HERE ####
# Create plots using the final model importance results
Deschampsia_data <- create_final_model_plots("Deschampsia", varimp_clf_final_combined, varimp_qrf_final_combined)


# SAVE FINAL MODEL OUTPUTS
saveRDS(all_results_clf_final, "hurdle_rf_classifier_results_Deschampsia_top_predictors.rds")
saveRDS(all_results_qrf_final, "hurdle_qrf_abundance_results_Deschampsia_top_predictors.rds")
saveRDS(varimp_list_clf_final, "hurdle_varimp_rf_presence_Deschampsia_top_predictors.rds")
saveRDS(varimp_list_qrf_final, "hurdle_varimp_qrf_abundance_Deschampsia_top_predictors.rds")
saveRDS(predictions_all_final, "hurdle_predictions_all_Deschampsia_top_predictors.rds")
saveRDS(training_data_list_final, "hurdle_training_data_list_Deschampsia_top_predictors.rds")
saveRDS(test_data_list_final, "hurdle_test_data_list_Deschampsia_top_predictors.rds")


#### QUASI-POISSON MODEL COMPARISON ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RUNNING QUASI-POISSON MODEL FOR COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

library(MASS)     # For stepwise selection
library(car)      # For VIF
library(broom)    # For tidy model output
library(pROC)     # For AUC calculation

# Use top predictor explaining 70% importance used in final model
predictors_presence <- c("DISTCOAST", "WEXP", "VDEP", "EFAH", "Easting", "MPI", "Northing", "WEFF", "LSF", "ASP", "HillSh", "PR_01", "MBI", "PDIS_CLASS", "RSP", "FLAC", "TAS_08")
predictors_abundance <- c("TPW", "LSF", "Northing", "MPI", "Easting", "ASP", "DISTCOAST", "EFAH", "MBI", "WEXP", "MNCURV", "FLAC", "SVF")

# Initialize storage for Quasi-Poisson results
all_results_qp <- list()
predictions_all_qp <- list()
model_summaries_qp <- list()

# Helper function to calculate Quasi-Poisson performance metrics
calc_qp_metrics <- function(predicted, observed, predicted_probs = NULL) {
  # For count data, we'll calculate both classification metrics (presence/absence) and regression metrics
  
  # Convert to presence/absence for classification metrics
  pred_presence <- ifelse(predicted > 0, 1, 0)
  obs_presence <- ifelse(observed > 0, 1, 0)
  
  # Classification metrics (presence/absence)
  acc <- mean(pred_presence == obs_presence)
  
  # Calculate confusion matrix elements
  tp <- sum(pred_presence == 1 & obs_presence == 1)
  tn <- sum(pred_presence == 0 & obs_presence == 0)
  fp <- sum(pred_presence == 1 & obs_presence == 0)
  fn <- sum(pred_presence == 0 & obs_presence == 1)
  
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  
  # AUC for presence/absence prediction using predicted counts as proxy
  auc <- tryCatch(pROC::auc(obs_presence, predicted, quiet = TRUE), 
                  error = function(e) NA)
  
  # Regression metrics for count data
  resids <- predicted - observed
  rmse <- sqrt(mean(resids^2, na.rm = TRUE))
  mae <- mean(abs(resids), na.rm = TRUE)
  bias <- mean(resids, na.rm = TRUE)
  
  # R-squared for count data
  r2 <- tryCatch({
    cor(observed, predicted, use = "complete.obs")^2
  }, error = function(e) NA)
  
  # Normalized MAE
  nmae <- ifelse(mean(observed, na.rm = TRUE) > 0, 
                 mae / mean(observed, na.rm = TRUE), NA)
  
  # Zero prediction metrics
  obs_zeros <- sum(observed == 0)
  pred_zeros <- sum(predicted == 0)
  zero_accuracy <- mean((observed == 0) == (predicted == 0))
  
  return(list(
    # Classification metrics
    acc = acc, sens = sens, spec = spec, auc = as.numeric(auc),
    # Regression metrics
    rmse = rmse, mae = mae, bias = bias, r2 = r2, nmae = nmae,
    # Zero prediction metrics
    obs_zeros = obs_zeros, pred_zeros = pred_zeros, zero_accuracy = zero_accuracy
  ))
}

cat("Running Quasi-Poisson model with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  set.seed(i)
  
  # Use same train-test splits as RF model for fair comparison
  if (i <= length(training_data_list_final) && i <= length(test_data_list_final)) {
    train <- training_data_list_final[[i]]
    test <- test_data_list_final[[i]]
  } else {
    # Fallback: create new split
    repeat {
      split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      # Check factor levels
      all_ok <- TRUE
      for (fac in factor_vars_presence) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
      if (all_ok && length(unique(train$Presence)) >= 2) break
    }
    
    # Align factor levels
    for (fac in factor_vars_presence) {
      if (fac %in% names(test)) {
        test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
      }
    }
  }
  
  ## Quasi-Poisson Model
  tryCatch({
    # Prepare count data - convert abundance to integer counts
    # You may need to scale/round your abundance data appropriately
    train$Deschampsia_count <- round(train$Deschampsia)
    test$Deschampsia_count <- round(test$Deschampsia)
    
    # Ensure non-negative integers
    train$Deschampsia_count <- pmax(0, train$Deschampsia_count)
    test$Deschampsia_count <- pmax(0, test$Deschampsia_count)
    
    # Create formula combining both presence and abundance predictors
    all_predictors <- unique(c(predictors_presence, predictors_abundance))
    qp_formula <- as.formula(paste("Deschampsia_count ~", paste(all_predictors, collapse = " + ")))
    
    # Fit Quasi-Poisson model
    qp_model <- glm(qp_formula, data = train, family = quasipoisson())
    
    # Check for convergence issues
    if (!qp_model$converged) {
      warning(paste("Quasi-Poisson model did not converge in iteration", i))
    }
    
    # Predict on test set
    pred_counts_qp <- predict(qp_model, newdata = test, type = "response")
    
    # Calculate metrics
    metrics_qp <- calc_qp_metrics(pred_counts_qp, test$Deschampsia_count)
    
    # Get dispersion parameter
    dispersion <- summary(qp_model)$dispersion
    
    all_results_qp[[i]] <- data.table(
      Iteration = i,
      # Classification metrics (presence/absence)
      Accuracy = metrics_qp$acc,
      Sensitivity = metrics_qp$sens,
      Specificity = metrics_qp$spec,
      AUC = metrics_qp$auc,
      # Regression metrics (count prediction)
      RMSE = metrics_qp$rmse,
      MAE = metrics_qp$mae,
      Bias = metrics_qp$bias,
      R_squared = metrics_qp$r2,
      nMAE = metrics_qp$nmae,
      # Zero prediction metrics
      Obs_Zeros = metrics_qp$obs_zeros,
      Pred_Zeros = metrics_qp$pred_zeros,
      Zero_Accuracy = metrics_qp$zero_accuracy,
      # Model diagnostics
      Dispersion = dispersion,
      Converged = qp_model$converged,
      N_Test = nrow(test),
      N_Train = nrow(train)
    )
    
    # Store predictions
    train_pred_counts <- predict(qp_model, newdata = train, type = "response")
    
    predictions_all_qp[[i]] <- rbind(
      data.table(set = "test", obs = test$Deschampsia_count, pred = pred_counts_qp, 
                 obs_orig = test$Deschampsia, iteration = i),
      data.table(set = "train", obs = train$Deschampsia_count, pred = train_pred_counts, 
                 obs_orig = train$Deschampsia, iteration = i)
    )
    
    # Store model for first iteration
    if (i == 1) {
      model_summaries_qp[["qp_model"]] <- list(
        model = qp_model,
        summary = summary(qp_model),
        aic = AIC(qp_model),
        formula = qp_formula,
        dispersion = dispersion,
        deviance = qp_model$deviance,
        null_deviance = qp_model$null.deviance
      )
    }
    
  }, error = function(e) {
    cat("Quasi-Poisson Error in iteration", i, ":", e$message, "\n")
    
    # Store NA results for failed iterations
    all_results_qp[[i]] <- data.table(
      Iteration = i, Accuracy = NA, Sensitivity = NA, Specificity = NA, AUC = NA,
      RMSE = NA, MAE = NA, Bias = NA, R_squared = NA, nMAE = NA,
      Obs_Zeros = NA, Pred_Zeros = NA, Zero_Accuracy = NA,
      Dispersion = NA, Converged = FALSE, N_Test = NA, N_Train = NA
    )
    
    predictions_all_qp[[i]] <- NULL
  })
  
  if (i %% 10 == 0) cat("Quasi-Poisson iteration", i, "completed\n")
}

# Combine Quasi-Poisson results
results_qp_combined <- rbindlist(all_results_qp)
predictions_qp_combined <- rbindlist(predictions_all_qp[!sapply(predictions_all_qp, is.null)])

#### MODEL COMPARISON ANALYSIS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RANDOM FOREST vs QUASI-POISSON COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# PRESENCE PREDICTION COMPARISON
cat("\n--- PRESENCE PREDICTION COMPARISON ---\n")
cat("Random Forest vs Quasi-Poisson\n\n")

rf_clf_summary <- results_clf_final_combined[, .(
  Model = "Random Forest",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

qp_clf_summary <- results_qp_combined[, .(
  Model = "Quasi-Poisson",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

presence_comparison <- rbind(rf_clf_summary, qp_clf_summary)
print(presence_comparison)

# ABUNDANCE COMPARISON
cat("\n--- ABUNDANCE PREDICTION COMPARISON ---\n")
cat("Quantile Random Forest vs Quasi-Poisson\n\n")

if (nrow(results_qrf_final_combined) > 0 && nrow(results_qp_combined) > 0) {
  rf_abundance_summary <- results_qrf_final_combined[, .(
    Model = "Quantile RF",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  qp_abundance_summary <- results_qp_combined[, .(
    Model = "Quasi-Poisson",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  abundance_comparison <- rbind(rf_abundance_summary, qp_abundance_summary)
  print(abundance_comparison)
  
  # Statistical tests for abundance metrics
  cat("\n--- STATISTICAL SIGNIFICANCE TESTS ---\n")
  
  # Presence metrics
  accuracy_test <- t.test(results_clf_final_combined$Accuracy, results_qp_combined$Accuracy)
  auc_test <- t.test(results_clf_final_combined$AUC, results_qp_combined$AUC)
  
  cat("Accuracy difference (RF - QP):\n")
  cat("  Mean difference:", round(accuracy_test$estimate[1] - accuracy_test$estimate[2], 4), "\n")
  cat("  p-value:", round(accuracy_test$p.value, 4), ifelse(accuracy_test$p.value < 0.05, "*", ""), "\n")
  
  cat("AUC difference (RF - QP):\n")
  cat("  Mean difference:", round(auc_test$estimate[1] - auc_test$estimate[2], 4), "\n")
  cat("  p-value:", round(auc_test$p.value, 4), ifelse(auc_test$p.value < 0.05, "*", ""), "\n")
  
  # Abundance metrics
  rf_rmse_valid <- results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]
  qp_rmse_valid <- results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)]
  
  if (length(rf_rmse_valid) > 10 && length(qp_rmse_valid) > 10) {
    rmse_test <- t.test(rf_rmse_valid, qp_rmse_valid)
    cat("RMSE difference (QRF - QP):\n")
    cat("  Mean difference:", round(rmse_test$estimate[1] - rmse_test$estimate[2], 4), "\n")
    cat("  p-value:", round(rmse_test$p.value, 4), ifelse(rmse_test$p.value < 0.05, "*", ""), "\n")
    
    rf_r2_valid <- results_qrf_final_combined$R_squared[!is.na(results_qrf_final_combined$R_squared)]
    qp_r2_valid <- results_qp_combined$R_squared[!is.na(results_qp_combined$R_squared)]
    
    if (length(rf_r2_valid) > 10 && length(qp_r2_valid) > 10) {  
      r2_test <- t.test(rf_r2_valid, qp_r2_valid)
      cat("R² difference (QRF - QP):\n")
      cat("  Mean difference:", round(r2_test$estimate[1] - r2_test$estimate[2], 4), "\n")
      cat("  p-value:", round(r2_test$p.value, 4), ifelse(r2_test$p.value < 0.05, "*", ""), "\n")
    }
  }
}

#### OVERDISPERSION ANALYSIS ####
cat("\n--- OVERDISPERSION ANALYSIS ---\n")
qp_dispersion_summary <- results_qp_combined[, .(
  Mean_Dispersion = round(mean(Dispersion, na.rm = TRUE), 3),
  SD_Dispersion = round(sd(Dispersion, na.rm = TRUE), 3),
  Min_Dispersion = round(min(Dispersion, na.rm = TRUE), 3),
  Max_Dispersion = round(max(Dispersion, na.rm = TRUE), 3),
  Convergence_Rate = round(mean(Converged, na.rm = TRUE), 3)
)]

cat("Overdispersion analysis:\n")
print(qp_dispersion_summary)

if (qp_dispersion_summary$Mean_Dispersion > 1.5) {
  cat("\nNote: Mean dispersion > 1.5 indicates overdispersion in the data\n")
} else if (qp_dispersion_summary$Mean_Dispersion < 0.8) {
  cat("\nNote: Mean dispersion < 0.8 indicates potential underdispersion\n")
} else {
  cat("\nNote: Dispersion close to 1, suggesting Poisson assumption may be reasonable\n")
}

#### MODEL INTERPRETATION ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== QUASI-POISSON MODEL INTERPRETATION ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

if (!is.null(model_summaries_qp[["qp_model"]])) {
  cat("\n--- QUASI-POISSON COEFFICIENTS ---\n")
  
  qp_summary <- model_summaries_qp[["qp_model"]]$summary
  
  # Model coefficients
  cat("\nQUASI-POISSON MODEL COEFFICIENTS:\n")
  coefs <- qp_summary$coefficients
  coef_df <- data.frame(
    Variable = rownames(coefs),
    Estimate = round(coefs[, "Estimate"], 4),
    Std_Error = round(coefs[, "Std. Error"], 4),
    t_value = round(coefs[, "t value"], 3),
    p_value = round(coefs[, "Pr(>|t|)"], 4),
    Significant = ifelse(coefs[, "Pr(>|t|)"] < 0.05, "*", "")
  )
  print(coef_df)
  
  cat("\nInterpretation: Coefficients represent log rate ratios\n")
  cat("Positive coefficients increase expected count, negative decrease it\n")
  
  # Model fit statistics
  cat("\nMODEL FIT STATISTICS:\n")
  cat("Dispersion parameter:", round(model_summaries_qp[["qp_model"]]$dispersion, 3), "\n")
  cat("Residual deviance:", round(model_summaries_qp[["qp_model"]]$deviance, 2), "\n")
  cat("Null deviance:", round(model_summaries_qp[["qp_model"]]$null_deviance, 2), "\n")
  
  # Pseudo R-squared
  pseudo_r2 <- 1 - (model_summaries_qp[["qp_model"]]$deviance / model_summaries_qp[["qp_model"]]$null_deviance)
  cat("Pseudo R-squared:", round(pseudo_r2, 3), "\n")
  
  cat("AIC:", round(model_summaries_qp[["qp_model"]]$aic, 2), "\n")
}

# RF Variable Importance comparison
cat("\n--- RANDOM FOREST VARIABLE IMPORTANCE ---\n")
cat("Top 10 variables for presence prediction:\n")
print(head(presence_final_importance[, c("Variable", "Importance_pct")], 10))

if (nrow(abundance_final_importance) > 0) {
  cat("\nTop 10 variables for abundance prediction:\n")
  print(head(abundance_final_importance[, c("Variable", "Importance_pct")], 10))
}

#### PREDICTION COMPARISON PLOTS ####
cat("\n=== CREATING COMPARISON PLOTS ===\n")

if (nrow(predictions_qp_combined) > 0 && nrow(predictions_final_combined) > 0) {
  # Combine predictions for plotting (using original scale)
  rf_pred_plot <- predictions_final_combined[set == "test", .(obs, pred, model = "Random Forest")]
  qp_pred_plot <- predictions_qp_combined[set == "test", .(obs = obs_orig, pred, model = "Quasi-Poisson")]
  
  combined_pred_plot <- rbind(rf_pred_plot, qp_pred_plot)
  
  # Scatter plot comparison
  p_comparison <- ggplot(combined_pred_plot, aes(x = obs, y = pred, color = model)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~model) +
    labs(
      title = "Deschampsia Abundance: Random Forest vs Quasi-Poisson",
      subtitle = "Test set predictions (perfect predictions would fall on diagonal line)",
      x = "Observed Abundance",
      y = "Predicted Abundance",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Random Forest" = "#F21A00", "Quasi-Poisson" = "#F21A50"))
  
  print(p_comparison)
  
  # Performance metrics comparison plot
  metrics_comparison <- rbind(
    data.table(Model = "Random Forest", 
               RMSE = results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]),
    data.table(Model = "Quasi-Poisson", 
               RMSE = results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)])
  )
  
  if (nrow(metrics_comparison) > 0) {
    p_metrics <- ggplot(metrics_comparison, aes(x = Model, y = RMSE, fill = Model)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = "RMSE Distribution: Random Forest vs Quasi-Poisson",
        subtitle = "Lower RMSE indicates better performance",
        y = "Root Mean Square Error"
      ) +
      theme_minimal() +
      scale_fill_manual(values = c("Random Forest" = "#F21A00", "Quasi-Poisson" = "#F21A50"))
    
    print(p_metrics)
  }
  
  # Residuals vs fitted plot for Quasi-Poisson
  if (nrow(predictions_qp_combined[set == "test"]) > 0) {
    qp_test_data <- predictions_qp_combined[set == "test"]
    qp_test_data[, residuals := obs_orig - pred]
    
    p_residuals <- ggplot(qp_test_data, aes(x = pred, y = residuals)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue") +
      labs(
        title = "Quasi-Poisson Model: Residuals vs Fitted Values",
        subtitle = "Test set predictions - should show random scatter around zero",
        x = "Fitted Values",
        y = "Residuals"
      ) +
      theme_minimal()
    
    print(p_residuals)
  }
}

#### SUMMARY AND RECOMMENDATIONS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== SUMMARY AND RECOMMENDATIONS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\nMODEL COMPARISON:\n")
cat("- Random Forest: Non-parametric, handles complex interactions, ensemble method\n")
cat("- Quasi-Poisson: Parametric, accounts for overdispersion, interpretable coefficients\n")

cat("\nQUASI-POISSON ADVANTAGES:\n")
cat("- Handles overdispersion in count data (dispersion parameter ≠ 1)\n")
cat("- More flexible than standard Poisson for overdispersed data\n")
cat("- Provides interpretable coefficients (log rate ratios)\n")
cat("- Well-established statistical foundation for count data\n")
cat("- Uses same mean-variance relationship as Poisson but allows different variance\n")

cat("\nPERFORMANCE SUMMARY:\n")
if (exists("presence_comparison")) {
  rf_better_acc <- presence_comparison$Mean_Accuracy[1] > presence_comparison$Mean_Accuracy[2]
  rf_better_auc <- presence_comparison$Mean_AUC[1] > presence_comparison$Mean_AUC[2]
  
  cat("Presence prediction:\n")
  cat("  - Better accuracy:", ifelse(rf_better_acc, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Better AUC:", ifelse(rf_better_auc, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("abundance_comparison") && nrow(abundance_comparison) > 0) {
  rf_better_rmse <- abundance_comparison$Mean_RMSE[1] < abundance_comparison$Mean_RMSE[2]
  rf_better_r2 <- abundance_comparison$Mean_R2[1] > abundance_comparison$Mean_R2[2]
  
  cat("Count prediction:\n")
  cat("  - Lower RMSE (better):", ifelse(rf_better_rmse, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Higher R² (better):", ifelse(rf_better_r2, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("qp_dispersion_summary")) {
  cat("Overdispersion:\n")
  cat("  - Mean dispersion parameter:", qp_dispersion_summary$Mean_Dispersion, "\n")
  if (qp_dispersion_summary$Mean_Dispersion > 1.2) {
    cat("  - Data shows overdispersion, Quasi-Poisson is appropriate\n")
  }
}

cat("\nRECOMMENDATION:\n")
cat("Choose Random Forest if: Complex non-linear relationships, prediction accuracy priority\n")
cat("Choose Quasi-Poisson if: Need interpretable coefficients, understand linear relationships\n")
cat("Consider Quasi-Poisson if: Count data with overdispersion, parametric approach preferred\n")










rm(list = ls())

#### Colobanthus ####
### COLOBANTHUS TRANSFORMATIONS FOR HURDLE MODEL - FINAL FIXED VERSION
### Handles extremely imbalanced data (19 positives out of 541 total)

library(randomForest)
library(quantregForest)
library(data.table)
library(dplyr)
library(rsample)
library(isotone)
library(pROC)
library(pdp)
library(clipr)
library(ggplot2)


select <- dplyr::select

# Load data
df <- read.csv("Colobanthus.csv")
df <- df %>% mutate(across(where(is.character), as.factor))
df$Presence <- as.factor(ifelse(df$Colobanthus > 0, 1, 0))

# Define transformation functions
transformations <- list(
  none = list(
    transform = function(x) x,
    inverse = function(x) x,
    name = "None"
  ),
  sqrt = list(
    transform = function(x) sqrt(pmax(x, 0)),
    inverse = function(x) pmax(x^2, 0),
    name = "Square Root"
  ),
  log = list(
    transform = function(x) log(pmax(x, 0.001)),
    inverse = function(x) exp(x),
    name = "Log"
  ),
  asinh = list(
    transform = function(x) asinh(x),
    inverse = function(x) sinh(x),
    name = "Asinh"
  )
)

predictors <- c("PDIS_CLASS", "ASP", "SLOPE_CLASS", "MORF_CLASS", "CD_CLASS", "DISTCOAST", 
                "VDEP", "GDGFGD0", "FCF", "SVF", "LSF", "MPI", "FLAC", "TXT", "TPW", 
                "Easting", "Northing")

factor_vars <- names(df[, predictors])[sapply(df[, predictors], is.factor)]
n_iter <- 100

# CONFIGURATION OPTIONS - Choose your approach:
USE_OPTIMAL_THRESHOLD <- TRUE    # Use ROC-optimal threshold instead of 0.5
USE_LOWER_THRESHOLD <- FALSE     # Use fixed 0.3 threshold
MIN_QRF_SAMPLES <- 3            # Minimum samples needed for QRF (was 10)

cat("Dataset info: Total samples =", nrow(df), 
    ", Positive cases =", sum(df$Presence == "1"), 
    "(", round(100 * mean(df$Presence == "1"), 1), "%)\n")

# Initialize storage
all_results_clf <- list()
all_results_qrf <- list()
varimp_list_clf <- list()
varimp_list_qrf <- list()
predictions_all <- list()
training_data_list <- list()
test_data_list <- list()
threshold_used <- list()  # Track thresholds used

for (trans_name in names(transformations)) {
  trans_func <- transformations[[trans_name]]
  
  cat("Testing transformation:", trans_func$name, "\n")
  
  # Transform the response variable
  df_trans <- df
  df_trans$Colobanthus_trans <- trans_func$transform(df$Colobanthus)
  
  all_results_clf[[trans_name]] <- list()
  all_results_qrf[[trans_name]] <- list()
  varimp_list_clf[[trans_name]] <- list()
  varimp_list_qrf[[trans_name]] <- list()
  predictions_all[[trans_name]] <- list()
  training_data_list[[trans_name]] <- list()
  test_data_list[[trans_name]] <- list()
  threshold_used[[trans_name]] <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    set.seed(i)
    repeat {
      # Check if ALT_CLASS exists
      if ("ALT_CLASS" %in% names(df_trans)) {
        split <- initial_split(df_trans, prop = 0.8, strata = ALT_CLASS)
      } else {
        split <- initial_split(df_trans, prop = 0.8)
      }
      
      train <- training(split)
      test <- testing(split)
      
      # Check factor level consistency
      all_ok <- TRUE
      for (fac in factor_vars) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
      
      # Ensure both presence classes exist in training data
      if (length(unique(train$Presence)) < 2) {
        all_ok <- FALSE
      }
      
      if (all_ok) break
    }
    
    # Align factor levels
    for (fac in factor_vars) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
    
    ## STEP 1: RF Classifier
    rf_clf <- randomForest(
      x = train[, predictors],
      y = train$Presence,
      ntree = 500,
      importance = TRUE,
      keep.forest = TRUE
    )
    
    pred_probs <- predict(rf_clf, newdata = test[, predictors], type = "prob")[, "1"]
    
    # Determine threshold
    threshold <- 0.5  # default
    
    if (USE_OPTIMAL_THRESHOLD) {
      # Use ROC-optimal threshold
      roc_obj <- pROC::roc(test$Presence, pred_probs, quiet = TRUE)
      coords_result <- pROC::coords(roc_obj, "best", ret = "threshold", quiet = TRUE)
      threshold <- coords_result$threshold[1]  # Take first if multiple
    } else if (USE_LOWER_THRESHOLD) {
      threshold <- 0.3
    }
    
    threshold_used[[trans_name]][i] <- threshold
    
    # Make predictions with chosen threshold
    pred_class <- factor(ifelse(pred_probs >= threshold, 1, 0), levels = levels(test$Presence))
    obs_class <- test$Presence
    
    # Calculate metrics
    acc <- mean(pred_class == obs_class)
    
    pred_char <- as.character(pred_class)
    obs_char <- as.character(obs_class)
    
    sens <- ifelse(sum(obs_char == "1") > 0, 
                   mean(pred_char == "1" & obs_char == "1") / mean(obs_char == "1"), 
                   NA)
    spec <- ifelse(sum(obs_char == "0") > 0, 
                   mean(pred_char == "0" & obs_char == "0") / mean(obs_char == "0"), 
                   NA)
    
    logloss <- -mean(ifelse(obs_char == "1", log(pmax(pred_probs, 1e-15)), log(pmax(1 - pred_probs, 1e-15))))
    brier <- mean((as.numeric(obs_char) - pred_probs)^2)
    auc <- tryCatch(pROC::auc(obs_char, pred_probs, quiet = TRUE), error = function(e) NA)
    
    spread_clf <- abs(0.5 - pred_probs)
    
    all_results_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
      LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
      MedianSpread = median(spread_clf), Threshold = threshold
    )
    
    imp_clf <- importance(rf_clf)
    varimp_list_clf[[trans_name]][[i]] <- data.table(
      Transformation = trans_func$name,
      Variable = rownames(imp_clf), 
      Importance = imp_clf[, 1], 
      Iteration = i, 
      Model = "RF_Presence"
    )
    
    ## STEP 2: QRF for predicted present cases
    test_presence_idx <- which(pred_class == "1")
    test_present_data <- test[test_presence_idx, ]
    train_present_data <- train[train$Presence == "1", ]
    
    n_test_pos <- nrow(test_present_data)
    n_train_pos <- nrow(train_present_data)
    
    # More lenient condition for QRF training
    if (n_test_pos >= MIN_QRF_SAMPLES && n_train_pos >= 5) {
      tryCatch({
        qrf <- quantregForest(
          x = train_present_data[, predictors],
          y = train_present_data$Colobanthus_trans,
          ntree = 500
        )
        
        preds_test_trans <- predict(qrf, newdata = test_present_data[, predictors], what = c(0.025, 0.5, 0.975))
        
        # Apply isotonic regression on transformed scale
        iso_model <- isoreg(preds_test_trans[, 2], test_present_data$Colobanthus_trans)
        iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
        calibrated_test_trans <- iso_fun(preds_test_trans[, 2])
        
        # Back-transform predictions and bounds
        preds_test_lower <- trans_func$inverse(preds_test_trans[, 1])
        preds_test_median <- trans_func$inverse(calibrated_test_trans)
        preds_test_upper <- trans_func$inverse(preds_test_trans[, 3])
        
        # Ensure reasonable bounds
        preds_test_lower <- pmax(pmin(preds_test_lower, 100), 0)
        preds_test_median <- pmax(pmin(preds_test_median, 100), 0)
        preds_test_upper <- pmax(pmin(preds_test_upper, 100), 0)
        
        spread_test <- preds_test_upper - preds_test_lower
        resids_test <- preds_test_median - test_present_data$Colobanthus
        
        rmse <- sqrt(mean(resids_test^2, na.rm = TRUE))
        mae <- mean(abs(resids_test), na.rm = TRUE)
        bias <- mean(resids_test, na.rm = TRUE)
        r2 <- cor(test_present_data$Colobanthus, preds_test_median, use = "complete.obs")^2
        nmae <- mae / mean(test_present_data$Colobanthus, na.rm = TRUE)
        coverage <- mean(test_present_data$Colobanthus >= preds_test_lower & 
                           test_present_data$Colobanthus <= preds_test_upper, na.rm = TRUE)
        median_spread <- median(spread_test, na.rm = TRUE)
        low_bias <- mean(preds_test_lower - test_present_data$Colobanthus, na.rm = TRUE)
        high_bias <- mean(preds_test_upper - test_present_data$Colobanthus, na.rm = TRUE)
        
        all_results_qrf[[trans_name]][[i]] <- data.table(
          Transformation = trans_func$name,
          Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
          R_squared = r2, nMAE = nmae, Coverage_95PI = coverage,
          Median_Spread_95PI = median_spread,
          Lower_Q2.5_Bias = low_bias,
          Upper_Q97.5_Bias = high_bias,
          N_TestPos = n_test_pos, N_TrainPos = n_train_pos
        )
        
        # Store predictions
        train_preds_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.5)
        train_spread_trans <- predict(qrf, newdata = train_present_data[, predictors], what = 0.975) -
          predict(qrf, newdata = train_present_data[, predictors], what = 0.025)
        
        predictions_all[[trans_name]][[i]] <- rbind(
          data.table(set = "test", obs = test_present_data$Colobanthus, pred = preds_test_median,
                     spread = spread_test, iteration = i, transformation = trans_func$name),
          data.table(set = "train", obs = train_present_data$Colobanthus,
                     pred = trans_func$inverse(train_preds_trans),
                     spread = abs(trans_func$inverse(train_present_data$Colobanthus_trans + train_spread_trans/2) - 
                                    trans_func$inverse(train_present_data$Colobanthus_trans - train_spread_trans/2)),
                     iteration = i, transformation = trans_func$name)
        )
        
        imp_qrf <- importance(qrf)
        varimp_list_qrf[[trans_name]][[i]] <- data.table(
          Transformation = trans_func$name,
          Variable = rownames(imp_qrf), 
          Importance = imp_qrf[, 1], 
          Iteration = i, 
          Model = "QRF_Abundance"
        )
        
      }, error = function(e) {
        cat("QRF Error in iteration", i, ":", e$message, "\n")
        all_results_qrf[[trans_name]][[i]] <- data.table(
          Transformation = trans_func$name,
          Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
          R_squared = NA, nMAE = NA, Coverage_95PI = NA,
          Median_Spread_95PI = NA, Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA,
          N_TestPos = n_test_pos, N_TrainPos = n_train_pos
        )
        predictions_all[[trans_name]][[i]] <- NULL
        varimp_list_qrf[[trans_name]][[i]] <- NULL
      })
    } else {
      # Insufficient data for QRF
      all_results_qrf[[trans_name]][[i]] <- data.table(
        Transformation = trans_func$name,
        Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
        R_squared = NA, nMAE = NA, Coverage_95PI = NA,
        Median_Spread_95PI = NA, Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA,
        N_TestPos = n_test_pos, N_TrainPos = n_train_pos
      )
      predictions_all[[trans_name]][[i]] <- NULL
      varimp_list_qrf[[trans_name]][[i]] <- NULL
    }
    
    training_data_list[[trans_name]][[i]] <- train
    test_data_list[[trans_name]][[i]] <- test
    
    if (i %% 10 == 0) {
      cat("Transformation:", trans_func$name, "- Iteration", i, 
          "- Pred pos:", n_test_pos, "- Threshold:", round(threshold, 3), "\n")
    }
  }
}

# Combine results
results_clf_combined <- rbindlist(lapply(all_results_clf, function(x) rbindlist(x)))
results_qrf_combined <- rbindlist(lapply(all_results_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_clf_combined <- rbindlist(lapply(varimp_list_clf, function(x) rbindlist(x[!sapply(x, is.null)])))
varimp_qrf_combined <- rbindlist(lapply(varimp_list_qrf, function(x) rbindlist(x[!sapply(x, is.null)])))
predictions_combined <- rbindlist(lapply(predictions_all, function(x) rbindlist(x[!sapply(x, is.null)])))

# Diagnostic code to understand why QRF results are NaN

# 1. Check the structure of your results
cat("=== DIAGNOSTICS ===\n")
cat("Structure of predictions_combined:\n")
print(str(predictions_combined))

cat("\nFirst few rows of predictions_combined:\n")
print(head(predictions_combined))

cat("\nNumber of rows in predictions_combined:\n")
print(nrow(predictions_combined))

# 2. Check class distribution in your original data
cat("\nOriginal class distribution:\n")
print(table(df$Presence))

cat("\nColobanthus abundance distribution:\n")
print(summary(df$Colobanthus[df$Colobanthus > 0]))

# 3. Let's see how many positive predictions we're getting per iteration
if (nrow(predictions_combined) > 0) {
  cat("\nPredicted presence cases per iteration/transformation:\n")
  presence_counts <- predictions_combined[, .N, by = .(transformation, iteration)]
  print(summary(presence_counts$N))
} else {
  cat("\nNo predictions found - this explains the NaN values!\n")
}

# 4. Alternative: Check a single iteration manually to see what's happening
set.seed(1)
if ("ALT_CLASS" %in% names(df)) {
  split_test <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
} else {
  split_test <- initial_split(df, prop = 0.8)
}

train_test <- training(split_test)
test_test <- testing(split_test)

# Quick RF to see prediction distribution
rf_test <- randomForest(
  x = train_test[, predictors],
  y = train_test$Presence,
  ntree = 100
)

pred_probs_test <- predict(rf_test, newdata = test_test[, predictors], type = "prob")[, "1"]
pred_class_test <- factor(ifelse(pred_probs_test >= 0.5, 1, 0), levels = levels(test_test$Presence))

cat("\nSample iteration prediction distribution:\n")
print(table(pred_class_test))
cat("Number of predicted positives:", sum(pred_class_test == "1"), "\n")
cat("This explains why nrow(test_present_data) > 10 condition fails!\n")

# 5. Suggestion: Lower the threshold or change the condition
cat("\n=== SUGGESTIONS ===\n")
cat("1. Your model is very conservative - consider lowering the threshold from 0.5\n")
cat("2. Or change the condition from nrow(test_present_data) > 10 to a lower number like 5\n")
cat("3. The low sensitivity (24%) suggests many true positives are being missed\n")

# 6. Try with a lower threshold
pred_class_test_lower <- factor(ifelse(pred_probs_test >= 0.3, 1, 0), levels = levels(test_test$Presence))
cat("\nWith 0.3 threshold:\n")
print(table(pred_class_test_lower))
cat("Number of predicted positives with 0.3 threshold:", sum(pred_class_test_lower == "1"), "\n")

# 7. Check if the spread plot error is due to empty predictions_combined
if (nrow(predictions_combined) == 0) {
  cat("\nThe spread plot error is because predictions_combined is empty!\n")
  cat("This is because no QRF models were trained due to insufficient predicted positives.\n")
} else {
  # Try to identify the spread plot issue
  cat("\nDebugging spread plot data:\n")
  tryCatch({
    spread_plot_data <- predictions_combined[set == "test", .(
      obs = obs,
      spread = spread,
      trans = transformation,
      iteration = iteration
    )]
    cat("Spread plot data created successfully, nrow =", nrow(spread_plot_data), "\n")
  }, error = function(e) {
    cat("Error in creating spread plot data:", e$message, "\n")
    cat("Check column names in predictions_combined:\n")
    print(names(predictions_combined))
  })
}

# Summary statistics by transformation
cat("\n=== CLASSIFICATION RESULTS SUMMARY ===\n")
clf_summary <- results_clf_combined[, .(
  Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
  SD_Accuracy = sd(Accuracy, na.rm = TRUE),
  Mean_AUC = mean(AUC, na.rm = TRUE),
  SD_AUC = sd(AUC, na.rm = TRUE),
  Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
  Mean_Specificity = mean(Specificity, na.rm = TRUE)
), by = Transformation]
print(clf_summary)

cat("\n=== QUANTILE REGRESSION RESULTS SUMMARY ===\n")
qrf_summary <- results_qrf_combined[, .(
  Mean_RMSE = mean(RMSE, na.rm = TRUE),
  SD_RMSE = sd(RMSE, na.rm = TRUE),
  Mean_MAE = mean(MAE, na.rm = TRUE),
  SD_MAE = sd(MAE, na.rm = TRUE),
  Mean_R2 = mean(R_squared, na.rm = TRUE),
  SD_R2 = sd(R_squared, na.rm = TRUE),
  Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
  Mean_Bias = mean(Bias, na.rm = TRUE)
), by = Transformation]
print(qrf_summary)

# 1. PREPARE DATA FOR PLOTTING

# Spread plot data (prediction interval width vs observed)
spread_plot_data <- predictions_combined[set == "test", .(
  obs = obs,
  spread = spread,
  trans = transformation,
  iteration = iteration
)]

# Observed vs predicted data
pred_obs_data <- predictions_combined[set == "test", .(
  obs = obs,
  pred = pred,
  trans = transformation,
  iteration = iteration
)]

# Calculate residuals for each transformation
residuals_data <- predictions_combined[set == "test", .(
  pred = pred,
  residuals = obs - pred,
  trans = transformation,
  iteration = iteration
)]

# 2. SPREAD PLOT (Uncertainty vs Observed)
p1 <- ggplot(spread_plot_data, aes(x = obs, y = spread, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Observed Deschampsia Abundance", y = "Prediction Interval Width",
       title = "Uncertainty vs Observed (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# 3. VARIABLE IMPORTANCE ANALYSIS
cat("\n=== VARIABLE IMPORTANCE ANALYSIS ===\n")

# Calculate variable importance summary
varimp_summary <- varimp_qrf_combined %>%
  group_by(Transformation, Variable) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 4),
            .groups = "drop") %>%
  arrange(Transformation, desc(MedianImportance))

# Calculate percentage importance
varimp_summary <- varimp_summary %>%
  group_by(Transformation) %>%
  mutate(Importance_pct = 100 * MedianImportance / sum(MedianImportance)) %>%
  ungroup()

# For "None" transformation (equivalent to "raw")
none_imp <- varimp_summary %>%
  filter(Transformation == "None") %>%
  mutate(Cumulative_Importance = round(cumsum(Importance_pct), 1)) %>%
  arrange(desc(MedianImportance))

cat("Variable Importance for No Transformation (Raw):\n")
print(none_imp)
clipr::write_clip(none_imp)

# Variable importance comparison plot
p2 <- varimp_summary %>%
  group_by(Transformation) %>%
  slice_head(n = 10) %>%  # Top 10 variables per transformation
  ggplot(aes(x = reorder(Variable, MedianImportance), y = MedianImportance, fill = Transformation)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Variable", y = "Median Importance", 
       title = "Top 10 Variable Importance by Transformation") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  facet_wrap(~Transformation, scales = "free_y")

print(p2)

# 4. OBSERVED VS PREDICTED PLOT
p3 <- ggplot(pred_obs_data, aes(x = obs, y = pred, color = trans)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed Colobanthus Abundance", y = "Predicted Colobanthus Abundance",
       title = "Observed vs Predicted (All Transformations)") +
  scale_color_brewer(palette = "Set1", name = "Transformation") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# 5. RESIDUALS VS PREDICTED (for each transformation)
p4 <- ggplot(residuals_data, aes(x = pred, y = residuals)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~trans, scales = "free") +
  theme_minimal() +
  labs(title = "Residuals vs Predicted (by Transformation)",
       x = "Predicted Colobanthus Abundance",
       y = "Residuals (Observed - Predicted)")

print(p4)

# 6. INDIVIDUAL TRANSFORMATION PLOTS

# Function to create plots for a specific transformation
create_transformation_plots <- function(trans_name) {
  
  # Filter data for specific transformation
  trans_spread <- spread_plot_data[trans == trans_name]
  trans_pred_obs <- pred_obs_data[trans == trans_name]
  trans_residuals <- residuals_data[trans == trans_name]
  
  cat(paste("\n=== PLOTS FOR", trans_name, "TRANSFORMATION ===\n"))
  
  # Residuals vs Predicted for specific transformation
  p_resid <- ggplot(trans_residuals, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Residuals vs Predicted -", trans_name, "Transformation"),
         x = "Predicted Colobanthus Abundance",
         y = "Residuals (Observed - Predicted)")
  
  # Prediction Interval Width vs Observed for specific transformation
  p_spread <- ggplot(trans_spread, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    theme_minimal() +
    labs(title = paste("Prediction Interval Width vs Observed -", trans_name, "Transformation"),
         x = "Observed Colobanthus Abundance",
         y = "95% Prediction Interval Width")
  
  print(p_resid)
  print(p_spread)
  
  return(list(residuals = p_resid, spread = p_spread))
}

# Create individual plots for each transformation
transformation_plots <- list()
for (trans in unique(spread_plot_data$trans)) {
  transformation_plots[[trans]] <- create_transformation_plots(trans)
}

# 7. PERFORMANCE SUMMARY COMPARISON PLOT
perf_summary_long <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI, Bias) %>%
  group_by(Transformation) %>%
  summarize(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Mean_MAE = mean(MAE, na.rm = TRUE),
    Mean_R2 = mean(R_squared, na.rm = TRUE),
    Mean_Coverage = mean(Coverage_95PI, na.rm = TRUE),
    Mean_Bias = abs(mean(Bias, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p5 <- ggplot(perf_summary_long, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_col() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Model Performance Comparison Across Transformations",
       x = "Transformation", y = "Mean Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p5)

# 8. BOXPLOT OF PERFORMANCE METRICS
perf_boxplot_data <- results_qrf_combined %>%
  select(Transformation, RMSE, MAE, R_squared, Coverage_95PI) %>%
  tidyr::pivot_longer(cols = -Transformation, names_to = "Metric", values_to = "Value")

p6 <- ggplot(perf_boxplot_data, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y") +
  labs(title = "Distribution of Performance Metrics Across Iterations",
       x = "Transformation", y = "Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p6)

# 9. CORRELATION BETWEEN OBSERVED AND PREDICTED BY TRANSFORMATION
cor_summary <- pred_obs_data %>%
  group_by(trans) %>%
  summarize(
    correlation = cor(obs, pred, use = "complete.obs"),
    rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    mae = mean(abs(obs - pred), na.rm = TRUE),
    .groups = "drop")
print(cor_summary)

#### EXTRACT TOP PREDICTORS FROM TRANSFORMATION ANALYSIS ####

# Extract top predictors from presence model (RF Classifier) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_presence <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

# Extract top predictors from abundance model (QRF) - Raw transformation
# Using 70% cumulative importance threshold
top_predictors_abundance <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%  # Using raw transformation
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct)
  ) %>%
  filter(Cumulative_pct <= 70) %>%  # Variables explaining up to 70% of importance
  pull(Variable)

cat("\n=== TOP PREDICTORS EXTRACTED FROM RAW TRANSFORMATION (70% IMPORTANCE) ===\n")
cat("Top predictors for PRESENCE model (explaining 70% of importance):\n")
print(top_predictors_presence)
cat("Number of presence predictors:", length(top_predictors_presence), "\n")

cat("\nTop predictors for ABUNDANCE model (explaining 70% of importance):\n")
print(top_predictors_abundance)
cat("Number of abundance predictors:", length(top_predictors_abundance), "\n")

# Optional: Show the importance breakdown
cat("\n=== IMPORTANCE BREAKDOWN FOR SELECTED PREDICTORS ===\n")

# Presence model importance breakdown
presence_importance_breakdown <- varimp_clf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_presence
  ) %>%
  filter(Selected)

cat("PRESENCE model - Selected variables and their importance:\n")
print(presence_importance_breakdown)

# Abundance model importance breakdown
abundance_importance_breakdown <- varimp_qrf_combined %>%
  filter(Transformation == "None") %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Importance_pct = 100 * MedianImportance / sum(MedianImportance),
    Cumulative_pct = cumsum(Importance_pct),
    Selected = Variable %in% top_predictors_abundance
  ) %>%
  filter(Selected)

cat("\nABUNDANCE model - Selected variables and their importance:\n")
print(abundance_importance_breakdown)





#### HURDLE FINAL MODEL USING TOP PREDICTORS ####

# Define the predictor variables for each model
predictors_presence <- c("VDEP", "Northing", "LSF", "Easting", "FCF", "SVF", "MPI", "DISTCOAST", "TXT")
predictors_abundance <- c("DISTCOAST", "VDEP", "TPW", "FLAC", "ASP")

# Identify factor variables for presence model
factor_vars_presence <- names(df[, predictors_presence])[sapply(df[, predictors_presence], is.factor)]

# Identify factor variables for abundance model  
factor_vars_abundance <- names(df[, predictors_abundance])[sapply(df[, predictors_abundance], is.factor)]

n_iter <- 100
all_results_clf_final <- list()
all_results_qrf_final <- list()
varimp_list_clf_final <- list()
varimp_list_qrf_final <- list()
predictions_all_final <- list()
training_data_list_final <- list()
test_data_list_final <- list()

cat("\n=== RUNNING FINAL MODEL WITH TOP PREDICTORS ===\n")

for (i in 1:n_iter) {
  set.seed(i)
  repeat {
    split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
    train <- training(split)
    test <- testing(split)
    
    # FIXED: Check factor levels correctly - test levels must exist in training
    all_ok <- TRUE
    
    # Check presence predictors
    for (fac in factor_vars_presence) {
      if (fac %in% names(train) && fac %in% names(test)) {
        train_levels <- levels(droplevels(train[[fac]]))
        test_unique <- unique(as.character(test[[fac]]))
        # FIXED: Check if all TEST levels exist in TRAINING levels (not the other way around)
        if (!all(test_unique %in% train_levels)) {
          all_ok <- FALSE
          break
        }
      }
    }
    
    # Check abundance predictors
    if (all_ok) {
      for (fac in factor_vars_abundance) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          # FIXED: Same correction here
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
    }
    
    # Ensure both presence classes exist in training data
    if (length(unique(train$Presence)) < 2) {
      all_ok <- FALSE
    }
    
    if (all_ok) break
  }
  
  # Set factor levels for both predictor sets
  for (fac in factor_vars_presence) {
    if (fac %in% names(test)) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
  }
  for (fac in factor_vars_abundance) {
    if (fac %in% names(test)) {
      test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
    }
  }
  
  ## STEP 1: RF Classifier using TOP PRESENCE predictors
  rf_clf <- randomForest(
    x = train[, predictors_presence],
    y = train$Presence,
    ntree = 500,
    importance = TRUE,
    keep.forest = TRUE
  )
  
  pred_probs <- predict(rf_clf, newdata = test[, predictors_presence], type = "prob")[, "1"]
  
  # FIXED: Ensure pred_class has same factor levels as obs_class
  pred_class <- factor(ifelse(pred_probs >= 0.5, 1, 0), levels = levels(test$Presence))
  obs_class <- test$Presence
  
  # FIXED: Calculate metrics safely to avoid factor level issues
  acc <- mean(pred_class == obs_class)
  
  # Convert to character to avoid factor level mismatches
  pred_char <- as.character(pred_class)
  obs_char <- as.character(obs_class)
  
  # FIXED: Calculate sensitivity and specificity correctly
  sens <- ifelse(sum(obs_char == "1") > 0, 
                 mean(pred_char == "1" & obs_char == "1") / mean(obs_char == "1"), 
                 NA)
  spec <- ifelse(sum(obs_char == "0") > 0, 
                 mean(pred_char == "0" & obs_char == "0") / mean(obs_char == "0"), 
                 NA)
  
  # FIXED: Prevent log(0) errors
  logloss <- -mean(ifelse(obs_char == "1", 
                          log(pmax(pred_probs, 1e-15)), 
                          log(pmax(1 - pred_probs, 1e-15))))
  brier <- mean((as.numeric(obs_char) - pred_probs)^2)
  auc <- tryCatch(pROC::auc(obs_char, pred_probs, quiet = TRUE), error = function(e) NA)
  
  spread_clf <- abs(0.5 - pred_probs)
  
  all_results_clf_final[[i]] <- data.table(
    Iteration = i, Accuracy = acc, Sensitivity = sens, Specificity = spec,
    LogLoss = logloss, BrierScore = brier, AUC = as.numeric(auc),
    MedianSpread = median(spread_clf)
  )
  
  imp_clf <- importance(rf_clf)
  varimp_list_clf_final[[i]] <- data.table(
    Variable = rownames(imp_clf), 
    Importance = imp_clf[, 1], 
    Iteration = i, 
    Model = "RF_Presence"
  )
  
  ## STEP 2: QRF using TOP ABUNDANCE predictors for predicted present cases
  test_presence_idx <- which(as.character(pred_class) == "1")  # FIXED: Use character comparison
  test_present_data <- test[test_presence_idx, ]
  train_present_data <- train[train$Presence == "1", ]
  
  # FIXED: Lower threshold for more robust QRF training
  if (nrow(test_present_data) > 3) {  # Changed from 10 to 3
    tryCatch({
      qrf <- quantregForest(
        x = train_present_data[, predictors_abundance],
        y = train_present_data$Colobanthus,
        ntree = 500
      )
      
      preds_test <- predict(qrf, newdata = test_present_data[, predictors_abundance], what = c(0.025, 0.5, 0.975))
      iso_model <- isoreg(preds_test[, 2], test_present_data$Colobanthus)
      iso_fun <- with(iso_model, approxfun(x, y, rule = 2))
      calibrated_test <- pmax(pmin(iso_fun(preds_test[, 2]), 100), 0)
      
      spread_test <- preds_test[, 3] - preds_test[, 1]
      resids_test <- calibrated_test - test_present_data$Colobanthus
      
      rmse <- sqrt(mean(resids_test^2, na.rm = TRUE))
      mae <- mean(abs(resids_test), na.rm = TRUE)
      bias <- mean(resids_test, na.rm = TRUE)
      r2 <- cor(test_present_data$Colobanthus, calibrated_test, use = "complete.obs")^2
      
      # Calculate nMAE and nRMSE using the same normalization approach
      mean_obs <- mean(test_present_data$Colobanthus, na.rm = TRUE)
      nmae <- mae / mean_obs
      nrmse <- rmse / mean_obs  # Added nRMSE calculation same as nMAE
      
      coverage <- mean(test_present_data$Colobanthus >= preds_test[, 1] & 
                         test_present_data$Colobanthus <= preds_test[, 3], na.rm = TRUE)
      median_spread <- median(spread_test, na.rm = TRUE)
      low_bias <- mean(preds_test[, 1] - test_present_data$Colobanthus, na.rm = TRUE)
      high_bias <- mean(preds_test[, 3] - test_present_data$Colobanthus, na.rm = TRUE)
      
      all_results_qrf_final[[i]] <- data.table(
        Iteration = i, RMSE = rmse, MAE = mae, Bias = bias,
        R_squared = r2, nMAE = nmae, nRMSE = nrmse,  # Added nRMSE here
        Coverage_95PI = coverage,
        Median_Spread_95PI = median_spread,
        Lower_Q2.5_Bias = low_bias,
        Upper_Q97.5_Bias = high_bias
      )
      
      # FIXED: Add error handling for predictions
      train_preds <- predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.5)
      train_spread <- predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.975) -
        predict(qrf, newdata = train_present_data[, predictors_abundance], what = 0.025)
      
      predictions_all_final[[i]] <- rbind(
        data.table(set = "test", obs = test_present_data$Colobanthus, pred = calibrated_test,
                   spread = spread_test, iteration = i),
        data.table(set = "train", obs = train_present_data$Colobanthus,
                   pred = train_preds, spread = train_spread, iteration = i)
      )
      
      imp_qrf <- importance(qrf)
      varimp_list_qrf_final[[i]] <- data.table(
        Variable = rownames(imp_qrf), 
        Importance = imp_qrf[, 1], 
        Iteration = i, 
        Model = "QRF_Abundance"
      )
      
    }, error = function(e) {
      cat("QRF Error in iteration", i, ":", e$message, "\n")
      all_results_qrf_final[[i]] <- data.table(
        Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
        R_squared = NA, nMAE = NA, nRMSE = NA,  # Added nRMSE here too
        Coverage_95PI = NA, Median_Spread_95PI = NA, 
        Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA
      )
      predictions_all_final[[i]] <- NULL
      varimp_list_qrf_final[[i]] <- NULL
    })
  } else {
    cat("Iteration", i, "- Insufficient predicted positives:", nrow(test_present_data), "\n")
    all_results_qrf_final[[i]] <- data.table(
      Iteration = i, RMSE = NA, MAE = NA, Bias = NA,
      R_squared = NA, nMAE = NA, nRMSE = NA,  # Added nRMSE here too
      Coverage_95PI = NA, Median_Spread_95PI = NA, 
      Lower_Q2.5_Bias = NA, Upper_Q97.5_Bias = NA
    )
    predictions_all_final[[i]] <- NULL
    varimp_list_qrf_final[[i]] <- NULL
  }
  
  training_data_list_final[[i]] <- train
  test_data_list_final[[i]] <- test
  
  if (i %% 10 == 0) cat("Final model iteration", i, "\n")
}

# Print summary of predictors used in final model
cat("\n=== PREDICTORS USED IN FINAL MODEL ===\n")
cat("Presence model predictors:", paste(predictors_presence, collapse = ", "), "\n")
cat("Abundance model predictors:", paste(predictors_abundance, collapse = ", "), "\n")

# Combine final results
results_clf_final_combined <- rbindlist(all_results_clf_final)
results_qrf_final_combined <- rbindlist(all_results_qrf_final[!sapply(all_results_qrf_final, is.null)])
varimp_clf_final_combined <- rbindlist(varimp_list_clf_final[!sapply(varimp_list_clf_final, is.null)])
varimp_qrf_final_combined <- rbindlist(varimp_list_qrf_final[!sapply(varimp_list_qrf_final, is.null)])
predictions_final_combined <- rbindlist(predictions_all_final[!sapply(predictions_all_final, is.null)])

# Final model performance summary - COMPLETE VERSION WITH ALL METRICS INCLUDING nRMSE
cat("\n=== FINAL MODEL PERFORMANCE SUMMARY ===\n")
cat("CLASSIFICATION RESULTS (using top presence predictors):\n")
clf_final_summary <- results_clf_final_combined[, .(
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 3),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 3),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 3),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 3),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 3),
  SD_Sensitivity = round(sd(Sensitivity, na.rm = TRUE), 3),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 3),
  SD_Specificity = round(sd(Specificity, na.rm = TRUE), 3),
  Mean_LogLoss = round(mean(LogLoss, na.rm = TRUE), 3),
  SD_LogLoss = round(sd(LogLoss, na.rm = TRUE), 3),
  Mean_BrierScore = round(mean(BrierScore, na.rm = TRUE), 3),
  SD_BrierScore = round(sd(BrierScore, na.rm = TRUE), 3),
  Mean_MedianSpread = round(mean(MedianSpread, na.rm = TRUE), 3),
  SD_MedianSpread = round(sd(MedianSpread, na.rm = TRUE), 3)
)]
print(clf_final_summary)

cat("\nQUANTILE REGRESSION RESULTS (using top abundance predictors):\n")
qrf_final_summary <- results_qrf_final_combined[, .(
  Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 3),
  SD_RMSE = round(sd(RMSE, na.rm = TRUE), 3),
  Mean_MAE = round(mean(MAE, na.rm = TRUE), 3),
  SD_MAE = round(sd(MAE, na.rm = TRUE), 3),
  Mean_R2 = round(mean(R_squared, na.rm = TRUE), 3),
  SD_R2 = round(sd(R_squared, na.rm = TRUE), 3),
  Mean_nMAE = round(mean(nMAE, na.rm = TRUE), 3),
  SD_nMAE = round(sd(nMAE, na.rm = TRUE), 3),
  Mean_nRMSE = round(mean(nRMSE, na.rm = TRUE), 3),  # Added nRMSE summary
  SD_nRMSE = round(sd(nRMSE, na.rm = TRUE), 3),      # Added nRMSE summary
  Mean_Bias = round(mean(Bias, na.rm = TRUE), 3),
  SD_Bias = round(sd(Bias, na.rm = TRUE), 3),
  Mean_Coverage_95PI = round(mean(Coverage_95PI, na.rm = TRUE), 3),
  SD_Coverage_95PI = round(sd(Coverage_95PI, na.rm = TRUE), 3),
  Mean_Median_Spread_95PI = round(mean(Median_Spread_95PI, na.rm = TRUE), 3),
  SD_Median_Spread_95PI = round(sd(Median_Spread_95PI, na.rm = TRUE), 3),
  Mean_Lower_Q2.5_Bias = round(mean(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  SD_Lower_Q2.5_Bias = round(sd(Lower_Q2.5_Bias, na.rm = TRUE), 3),
  Mean_Upper_Q97.5_Bias = round(mean(Upper_Q97.5_Bias, na.rm = TRUE), 3),
  SD_Upper_Q97.5_Bias = round(sd(Upper_Q97.5_Bias, na.rm = TRUE), 3)
)]
print(qrf_final_summary)


#### VARIABLE IMPORTANCE ANALYSIS FOR FINAL MODELS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== VARIABLE IMPORTANCE FROM FINAL MODELS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create comprehensive importance tables from the final model runs
cat("\nTABLE 1: PRESENCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

presence_final_importance <- varimp_clf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(presence_final_importance)

cat("\nPRESENCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(presence_final_importance), "\n"))
cat(paste("- Top variable (", presence_final_importance$Variable[1], ") explains:", 
          presence_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", presence_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 2: ABUNDANCE MODEL - Variable Importance from Final Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

abundance_final_importance <- varimp_qrf_final_combined %>%
  group_by(Variable) %>%
  summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MedianImportance)) %>%
  mutate(
    Rank = row_number(),
    Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
    Cumulative_pct = round(cumsum(Importance_pct), 2)
  ) %>%
  select(Rank, Variable, MedianImportance, Importance_pct, Cumulative_pct)

print(abundance_final_importance)

cat("\nABUNDANCE MODEL SUMMARY:\n")
cat(paste("- Total variables in final model:", nrow(abundance_final_importance), "\n"))
cat(paste("- Top variable (", abundance_final_importance$Variable[1], ") explains:", 
          abundance_final_importance$Importance_pct[1], "% of importance\n"))
cat(paste("- Top 3 variables explain:", abundance_final_importance$Cumulative_pct[3], "% of importance\n"))
cat(paste("- All variables explain: 100% of importance\n"))

# Create a comparison table showing which variables are most important in each model
cat("\n", paste(rep("-", 80), collapse=""), "\n")
cat("TABLE 3: COMPARISON - Top Variables in Each Model\n")
cat(paste(rep("-", 80), collapse=""), "\n")

# Get top 5 from each model for comparison
top_presence_vars <- presence_final_importance$Variable[1:min(5, nrow(presence_final_importance))]
top_abundance_vars <- abundance_final_importance$Variable[1:min(5, nrow(abundance_final_importance))]

# Create comparison dataframe
comparison_table <- data.frame(
  Rank = 1:5,
  Presence_Model = c(top_presence_vars, rep("", 5 - length(top_presence_vars)))[1:5],
  Presence_Importance = c(presence_final_importance$Importance_pct[1:min(5, nrow(presence_final_importance))], 
                          rep(NA, 5 - min(5, nrow(presence_final_importance))))[1:5],
  Abundance_Model = c(top_abundance_vars, rep("", 5 - length(top_abundance_vars)))[1:5],
  Abundance_Importance = c(abundance_final_importance$Importance_pct[1:min(5, nrow(abundance_final_importance))], 
                           rep(NA, 5 - min(5, nrow(abundance_final_importance))))[1:5]
)

print(comparison_table)

# Identify shared important variables
shared_vars <- intersect(top_presence_vars, top_abundance_vars)
if(length(shared_vars) > 0) {
  cat("\nShared important variables between models:", paste(shared_vars, collapse = ", "), "\n")
} else {
  cat("\nNo shared variables in top 5 between models\n")
}

#### PLOTTING FUNCTIONS - RUN THIS FIRST ####

# Define color palette for species
cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00",
           "Colobanthus" = "#78B7C5")

library(ggplot2)
library(dplyr)

# Function to extract 70% importance variables from FINAL MODEL results
extract_final_70_percent <- function(varimp_final_combined, species_name, model_type) {
  
  # Calculate importance from final model runs and apply 70% threshold
  importance_data <- varimp_final_combined %>%
    group_by(Variable) %>%
    summarize(MedianImportance = median(Importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MedianImportance)) %>%
    mutate(
      Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 2),
      Cumulative_pct = round(cumsum(Importance_pct), 2),
      Selected_70pct = Cumulative_pct <= 70,
      Rank = row_number(),
      Species = species_name,
      Model = model_type
    ) %>%
    filter(Selected_70pct) %>%  # Only keep variables that explain ≤70% cumulative importance
    select(Species, Model, Variable, Rank, Importance_pct, Cumulative_pct)
  
  return(importance_data)
}

create_individual_plot <- function(species_name, model_type, data_all) {
  
  data_subset <- data_all %>%
    filter(Species == species_name, Model == model_type) %>%
    arrange(desc(Rank))  # Most important at top
  
  if(nrow(data_subset) == 0) {
    cat(paste("No data found for", species_name, model_type, "\n"))
    return(NULL)
  }
  
  p <- ggplot(data_subset, aes(x = reorder(Variable, Rank), y = Importance_pct)) +
    geom_col(fill = cores[species_name], alpha = 0.8, color = "white", size = 0.5) +
    geom_text(aes(label = paste0(round(Importance_pct, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    labs(
      title = paste(species_name, "-", model_type, "Final Model (70% Variables)"),
      subtitle = paste("Variables explaining", round(max(data_subset$Cumulative_pct), 1), "% importance in FINAL model"),
      x = "Variables (Ranked by Final Model Importance)",
      y = "Individual Importance (%)",
      caption = paste("n =", nrow(data_subset), "variables from final model")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
      plot.subtitle = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(p)
}

create_final_model_plots <- function(species_name, varimp_clf_final, varimp_qrf_final) {
  
  cat(paste("\n=== Creating plots for", species_name, "final model ===\n"))
  
  # Extract 70% variables from final models
  presence_data <- extract_final_70_percent(varimp_clf_final, species_name, "Presence")
  abundance_data <- extract_final_70_percent(varimp_qrf_final, species_name, "Abundance")
  
  # Combine data
  combined_data <- rbind(presence_data, abundance_data)
  
  # Create individual plots
  presence_plot <- create_individual_plot(species_name, "Presence", combined_data)
  abundance_plot <- create_individual_plot(species_name, "Abundance", combined_data)
  
  # Print plots
  if(!is.null(presence_plot)) {
    print(presence_plot)
    cat(paste("Created", species_name, "presence plot with", nrow(presence_data), "variables\n"))
  }
  
  if(!is.null(abundance_plot)) {
    print(abundance_plot) 
    cat(paste("Created", species_name, "abundance plot with", nrow(abundance_data), "variables\n"))
  }
  
  # Create comparison plot
  if(nrow(combined_data) > 0) {
    comparison_plot <- ggplot(combined_data, aes(x = reorder(Variable, Rank), y = Importance_pct, fill = Model)) +
      geom_col(position = "dodge", alpha = 0.8, color = "white", size = 0.5) +
      coord_flip() +
      labs(
        title = paste(species_name, "- Final Model: Presence vs Abundance (70% Variables)"),
        subtitle = "Variable importance from final models (70% threshold applied)",
        x = "Variables",
        y = "Individual Importance (%)",
        fill = "Model Type"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = cores[species_name]),
        plot.subtitle = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      ) +
      scale_fill_manual(values = c("Presence" = cores[species_name], "Abundance" = adjustcolor(cores[species_name], alpha.f = 0.6))) +
      facet_wrap(~Model, scales = "free_y", ncol = 2)
    
    print(comparison_plot)
    cat(paste("Created", species_name, "comparison plot\n"))
  }
  
  return(combined_data)
}
#### INSERT THE PLOTTING CODE HERE ####
# Create plots using the final model importance results
Colobanthus_data <- create_final_model_plots("Colobanthus", varimp_clf_final_combined, varimp_qrf_final_combined)



# SAVE FINAL MODEL OUTPUTS
saveRDS(all_results_clf_final, "hurdle_rf_classifier_results_Colobanthus_top_predictors.rds")
saveRDS(all_results_qrf_final, "hurdle_qrf_abundance_results_Colobanthus_top_predictors.rds")
saveRDS(varimp_list_clf_final, "hurdle_varimp_rf_presence_Colobanthus_top_predictors.rds")
saveRDS(varimp_list_qrf_final, "hurdle_varimp_qrf_abundance_Colobanthus_top_predictors.rds")
saveRDS(predictions_all_final, "hurdle_predictions_all_Colobanthus_top_predictors.rds")
saveRDS(training_data_list_final, "hurdle_training_data_list_Colobanthus_top_predictors.rds")
saveRDS(test_data_list_final, "hurdle_test_data_list_Colobanthus_top_predictors.rds")


#### QUASI-POISSON MODEL COMPARISON ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RUNNING QUASI-POISSON MODEL FOR COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

library(MASS)     # For stepwise selection
library(car)      # For VIF
library(broom)    # For tidy model output
library(pROC)     # For AUC calculation

# Initialize storage for Quasi-Poisson results
all_results_qp <- list()
predictions_all_qp <- list()
model_summaries_qp <- list()

# Helper function to calculate Quasi-Poisson performance metrics
calc_qp_metrics <- function(predicted, observed, predicted_probs = NULL) {
  # For count data, we'll calculate both classification metrics (presence/absence) and regression metrics
  
  # Convert to presence/absence for classification metrics
  pred_presence <- ifelse(predicted > 0, 1, 0)
  obs_presence <- ifelse(observed > 0, 1, 0)
  
  # Classification metrics (presence/absence)
  acc <- mean(pred_presence == obs_presence)
  
  # Calculate confusion matrix elements
  tp <- sum(pred_presence == 1 & obs_presence == 1)
  tn <- sum(pred_presence == 0 & obs_presence == 0)
  fp <- sum(pred_presence == 1 & obs_presence == 0)
  fn <- sum(pred_presence == 0 & obs_presence == 1)
  
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  
  # AUC for presence/absence prediction using predicted counts as proxy
  auc <- tryCatch(pROC::auc(obs_presence, predicted, quiet = TRUE), 
                  error = function(e) NA)
  
  # Regression metrics for count data
  resids <- predicted - observed
  rmse <- sqrt(mean(resids^2, na.rm = TRUE))
  mae <- mean(abs(resids), na.rm = TRUE)
  bias <- mean(resids, na.rm = TRUE)
  
  # R-squared for count data
  r2 <- tryCatch({
    cor(observed, predicted, use = "complete.obs")^2
  }, error = function(e) NA)
  
  # Normalized MAE
  nmae <- ifelse(mean(observed, na.rm = TRUE) > 0, 
                 mae / mean(observed, na.rm = TRUE), NA)
  
  # Zero prediction metrics
  obs_zeros <- sum(observed == 0)
  pred_zeros <- sum(predicted == 0)
  zero_accuracy <- mean((observed == 0) == (predicted == 0))
  
  return(list(
    # Classification metrics
    acc = acc, sens = sens, spec = spec, auc = as.numeric(auc),
    # Regression metrics
    rmse = rmse, mae = mae, bias = bias, r2 = r2, nmae = nmae,
    # Zero prediction metrics
    obs_zeros = obs_zeros, pred_zeros = pred_zeros, zero_accuracy = zero_accuracy
  ))
}

cat("Running Quasi-Poisson model with", n_iter, "iterations...\n")

for (i in 1:n_iter) {
  set.seed(i)
  
  # Use same train-test splits as RF model for fair comparison
  if (i <= length(training_data_list_final) && i <= length(test_data_list_final)) {
    train <- training_data_list_final[[i]]
    test <- test_data_list_final[[i]]
  } else {
    # Fallback: create new split
    repeat {
      split <- initial_split(df, prop = 0.8, strata = ALT_CLASS)
      train <- training(split)
      test <- testing(split)
      
      # Check factor levels
      all_ok <- TRUE
      for (fac in factor_vars_presence) {
        if (fac %in% names(train) && fac %in% names(test)) {
          train_levels <- levels(droplevels(train[[fac]]))
          test_unique <- unique(as.character(test[[fac]]))
          if (!all(test_unique %in% train_levels)) {
            all_ok <- FALSE
            break
          }
        }
      }
      if (all_ok && length(unique(train$Presence)) >= 2) break
    }
    
    # Align factor levels
    for (fac in factor_vars_presence) {
      if (fac %in% names(test)) {
        test[[fac]] <- factor(test[[fac]], levels = levels(train[[fac]]))
      }
    }
  }
  
  ## Quasi-Poisson Model
  tryCatch({
    # Prepare count data - convert abundance to integer counts
    # You may need to scale/round your abundance data appropriately
    train$Colobanthus_count <- round(train$Colobanthus)
    test$Colobanthus_count <- round(test$Colobanthus)
    
    # Ensure non-negative integers
    train$Colobanthus_count <- pmax(0, train$Colobanthus_count)
    test$Colobanthus_count <- pmax(0, test$Colobanthus_count)
    
    # Create formula combining both presence and abundance predictors
    all_predictors <- unique(c(predictors_presence, predictors_abundance))
    qp_formula <- as.formula(paste("Colobanthus_count ~", paste(all_predictors, collapse = " + ")))
    
    # Fit Quasi-Poisson model
    qp_model <- glm(qp_formula, data = train, family = quasipoisson())
    
    # Check for convergence issues
    if (!qp_model$converged) {
      warning(paste("Quasi-Poisson model did not converge in iteration", i))
    }
    
    # Predict on test set
    pred_counts_qp <- predict(qp_model, newdata = test, type = "response")
    
    # Calculate metrics
    metrics_qp <- calc_qp_metrics(pred_counts_qp, test$Colobanthus_count)
    
    # Get dispersion parameter
    dispersion <- summary(qp_model)$dispersion
    
    all_results_qp[[i]] <- data.table(
      Iteration = i,
      # Classification metrics (presence/absence)
      Accuracy = metrics_qp$acc,
      Sensitivity = metrics_qp$sens,
      Specificity = metrics_qp$spec,
      AUC = metrics_qp$auc,
      # Regression metrics (count prediction)
      RMSE = metrics_qp$rmse,
      MAE = metrics_qp$mae,
      Bias = metrics_qp$bias,
      R_squared = metrics_qp$r2,
      nMAE = metrics_qp$nmae,
      # Zero prediction metrics
      Obs_Zeros = metrics_qp$obs_zeros,
      Pred_Zeros = metrics_qp$pred_zeros,
      Zero_Accuracy = metrics_qp$zero_accuracy,
      # Model diagnostics
      Dispersion = dispersion,
      Converged = qp_model$converged,
      N_Test = nrow(test),
      N_Train = nrow(train)
    )
    
    # Store predictions
    train_pred_counts <- predict(qp_model, newdata = train, type = "response")
    
    predictions_all_qp[[i]] <- rbind(
      data.table(set = "test", obs = test$Colobanthus_count, pred = pred_counts_qp, 
                 obs_orig = test$Colobanthus, iteration = i),
      data.table(set = "train", obs = train$Colobanthus_count, pred = train_pred_counts, 
                 obs_orig = train$Colobanthus, iteration = i)
    )
    
    # Store model for first iteration
    if (i == 1) {
      model_summaries_qp[["qp_model"]] <- list(
        model = qp_model,
        summary = summary(qp_model),
        aic = AIC(qp_model),
        formula = qp_formula,
        dispersion = dispersion,
        deviance = qp_model$deviance,
        null_deviance = qp_model$null.deviance
      )
    }
    
  }, error = function(e) {
    cat("Quasi-Poisson Error in iteration", i, ":", e$message, "\n")
    
    # Store NA results for failed iterations
    all_results_qp[[i]] <- data.table(
      Iteration = i, Accuracy = NA, Sensitivity = NA, Specificity = NA, AUC = NA,
      RMSE = NA, MAE = NA, Bias = NA, R_squared = NA, nMAE = NA,
      Obs_Zeros = NA, Pred_Zeros = NA, Zero_Accuracy = NA,
      Dispersion = NA, Converged = FALSE, N_Test = NA, N_Train = NA
    )
    
    predictions_all_qp[[i]] <- NULL
  })
  
  if (i %% 10 == 0) cat("Quasi-Poisson iteration", i, "completed\n")
}

# Combine Quasi-Poisson results
results_qp_combined <- rbindlist(all_results_qp)
predictions_qp_combined <- rbindlist(predictions_all_qp[!sapply(predictions_all_qp, is.null)])

#### MODEL COMPARISON ANALYSIS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== RANDOM FOREST vs QUASI-POISSON COMPARISON ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# PRESENCE PREDICTION COMPARISON
cat("\n--- PRESENCE PREDICTION COMPARISON ---\n")
cat("Random Forest vs Quasi-Poisson\n\n")

rf_clf_summary <- results_clf_final_combined[, .(
  Model = "Random Forest",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

qp_clf_summary <- results_qp_combined[, .(
  Model = "Quasi-Poisson",
  Mean_Accuracy = round(mean(Accuracy, na.rm = TRUE), 4),
  SD_Accuracy = round(sd(Accuracy, na.rm = TRUE), 4),
  Mean_AUC = round(mean(AUC, na.rm = TRUE), 4),
  SD_AUC = round(sd(AUC, na.rm = TRUE), 4),
  Mean_Sensitivity = round(mean(Sensitivity, na.rm = TRUE), 4),
  Mean_Specificity = round(mean(Specificity, na.rm = TRUE), 4)
)]

presence_comparison <- rbind(rf_clf_summary, qp_clf_summary)
print(presence_comparison)

# ABUNDANCE COMPARISON
cat("\n--- ABUNDANCE PREDICTION COMPARISON ---\n")
cat("Quantile Random Forest vs Quasi-Poisson\n\n")

if (nrow(results_qrf_final_combined) > 0 && nrow(results_qp_combined) > 0) {
  rf_abundance_summary <- results_qrf_final_combined[, .(
    Model = "Quantile RF",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  qp_abundance_summary <- results_qp_combined[, .(
    Model = "Quasi-Poisson",
    Mean_RMSE = round(mean(RMSE, na.rm = TRUE), 4),
    SD_RMSE = round(sd(RMSE, na.rm = TRUE), 4),
    Mean_MAE = round(mean(MAE, na.rm = TRUE), 4),
    SD_MAE = round(sd(MAE, na.rm = TRUE), 4),
    Mean_R2 = round(mean(R_squared, na.rm = TRUE), 4),
    SD_R2 = round(sd(R_squared, na.rm = TRUE), 4),
    Mean_Bias = round(mean(Bias, na.rm = TRUE), 4),
    Valid_Models = sum(!is.na(RMSE))
  )]
  
  abundance_comparison <- rbind(rf_abundance_summary, qp_abundance_summary)
  print(abundance_comparison)
  
  # Statistical tests for abundance metrics
  cat("\n--- STATISTICAL SIGNIFICANCE TESTS ---\n")
  
  # Presence metrics
  accuracy_test <- t.test(results_clf_final_combined$Accuracy, results_qp_combined$Accuracy)
  auc_test <- t.test(results_clf_final_combined$AUC, results_qp_combined$AUC)
  
  cat("Accuracy difference (RF - QP):\n")
  cat("  Mean difference:", round(accuracy_test$estimate[1] - accuracy_test$estimate[2], 4), "\n")
  cat("  p-value:", round(accuracy_test$p.value, 4), ifelse(accuracy_test$p.value < 0.05, "*", ""), "\n")
  
  cat("AUC difference (RF - QP):\n")
  cat("  Mean difference:", round(auc_test$estimate[1] - auc_test$estimate[2], 4), "\n")
  cat("  p-value:", round(auc_test$p.value, 4), ifelse(auc_test$p.value < 0.05, "*", ""), "\n")
  
  # Abundance metrics
  rf_rmse_valid <- results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]
  qp_rmse_valid <- results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)]
  
  if (length(rf_rmse_valid) > 10 && length(qp_rmse_valid) > 10) {
    rmse_test <- t.test(rf_rmse_valid, qp_rmse_valid)
    cat("RMSE difference (QRF - QP):\n")
    cat("  Mean difference:", round(rmse_test$estimate[1] - rmse_test$estimate[2], 4), "\n")
    cat("  p-value:", round(rmse_test$p.value, 4), ifelse(rmse_test$p.value < 0.05, "*", ""), "\n")
    
    rf_r2_valid <- results_qrf_final_combined$R_squared[!is.na(results_qrf_final_combined$R_squared)]
    qp_r2_valid <- results_qp_combined$R_squared[!is.na(results_qp_combined$R_squared)]
    
    if (length(rf_r2_valid) > 10 && length(qp_r2_valid) > 10) {  
      r2_test <- t.test(rf_r2_valid, qp_r2_valid)
      cat("R² difference (QRF - QP):\n")
      cat("  Mean difference:", round(r2_test$estimate[1] - r2_test$estimate[2], 4), "\n")
      cat("  p-value:", round(r2_test$p.value, 4), ifelse(r2_test$p.value < 0.05, "*", ""), "\n")
    }
  }
}

#### OVERDISPERSION ANALYSIS ####
cat("\n--- OVERDISPERSION ANALYSIS ---\n")
qp_dispersion_summary <- results_qp_combined[, .(
  Mean_Dispersion = round(mean(Dispersion, na.rm = TRUE), 3),
  SD_Dispersion = round(sd(Dispersion, na.rm = TRUE), 3),
  Min_Dispersion = round(min(Dispersion, na.rm = TRUE), 3),
  Max_Dispersion = round(max(Dispersion, na.rm = TRUE), 3),
  Convergence_Rate = round(mean(Converged, na.rm = TRUE), 3)
)]

cat("Overdispersion analysis:\n")
print(qp_dispersion_summary)

if (qp_dispersion_summary$Mean_Dispersion > 1.5) {
  cat("\nNote: Mean dispersion > 1.5 indicates overdispersion in the data\n")
} else if (qp_dispersion_summary$Mean_Dispersion < 0.8) {
  cat("\nNote: Mean dispersion < 0.8 indicates potential underdispersion\n")
} else {
  cat("\nNote: Dispersion close to 1, suggesting Poisson assumption may be reasonable\n")
}

#### MODEL INTERPRETATION ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== QUASI-POISSON MODEL INTERPRETATION ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

if (!is.null(model_summaries_qp[["qp_model"]])) {
  cat("\n--- QUASI-POISSON COEFFICIENTS ---\n")
  
  qp_summary <- model_summaries_qp[["qp_model"]]$summary
  
  # Model coefficients
  cat("\nQUASI-POISSON MODEL COEFFICIENTS:\n")
  coefs <- qp_summary$coefficients
  coef_df <- data.frame(
    Variable = rownames(coefs),
    Estimate = round(coefs[, "Estimate"], 4),
    Std_Error = round(coefs[, "Std. Error"], 4),
    t_value = round(coefs[, "t value"], 3),
    p_value = round(coefs[, "Pr(>|t|)"], 4),
    Significant = ifelse(coefs[, "Pr(>|t|)"] < 0.05, "*", "")
  )
  print(coef_df)
  
  cat("\nInterpretation: Coefficients represent log rate ratios\n")
  cat("Positive coefficients increase expected count, negative decrease it\n")
  
  # Model fit statistics
  cat("\nMODEL FIT STATISTICS:\n")
  cat("Dispersion parameter:", round(model_summaries_qp[["qp_model"]]$dispersion, 3), "\n")
  cat("Residual deviance:", round(model_summaries_qp[["qp_model"]]$deviance, 2), "\n")
  cat("Null deviance:", round(model_summaries_qp[["qp_model"]]$null_deviance, 2), "\n")
  
  # Pseudo R-squared
  pseudo_r2 <- 1 - (model_summaries_qp[["qp_model"]]$deviance / model_summaries_qp[["qp_model"]]$null_deviance)
  cat("Pseudo R-squared:", round(pseudo_r2, 3), "\n")
  
  cat("AIC:", round(model_summaries_qp[["qp_model"]]$aic, 2), "\n")
}

# RF Variable Importance comparison
cat("\n--- RANDOM FOREST VARIABLE IMPORTANCE ---\n")
cat("Top 10 variables for presence prediction:\n")
print(head(presence_final_importance[, c("Variable", "Importance_pct")], 10))

if (nrow(abundance_final_importance) > 0) {
  cat("\nTop 10 variables for abundance prediction:\n")
  print(head(abundance_final_importance[, c("Variable", "Importance_pct")], 10))
}

#### PREDICTION COMPARISON PLOTS ####
cat("\n=== CREATING COMPARISON PLOTS ===\n")

if (nrow(predictions_qp_combined) > 0 && nrow(predictions_final_combined) > 0) {
  # Combine predictions for plotting (using original scale)
  rf_pred_plot <- predictions_final_combined[set == "test", .(obs, pred, model = "Random Forest")]
  qp_pred_plot <- predictions_qp_combined[set == "test", .(obs = obs_orig, pred, model = "Quasi-Poisson")]
  
  combined_pred_plot <- rbind(rf_pred_plot, qp_pred_plot)
  
  # Scatter plot comparison
  p_comparison <- ggplot(combined_pred_plot, aes(x = obs, y = pred, color = model)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~model) +
    labs(
      title = "Colobanthus Abundance: Random Forest vs Quasi-Poisson",
      subtitle = "Test set predictions (perfect predictions would fall on diagonal line)",
      x = "Observed Abundance",
      y = "Predicted Abundance",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Random Forest" = "#78B7C5", "Quasi-Poisson" = "#F21A00"))
  
  print(p_comparison)
  
  # Performance metrics comparison plot
  metrics_comparison <- rbind(
    data.table(Model = "Random Forest", 
               RMSE = results_qrf_final_combined$RMSE[!is.na(results_qrf_final_combined$RMSE)]),
    data.table(Model = "Quasi-Poisson", 
               RMSE = results_qp_combined$RMSE[!is.na(results_qp_combined$RMSE)])
  )
  
  if (nrow(metrics_comparison) > 0) {
    p_metrics <- ggplot(metrics_comparison, aes(x = Model, y = RMSE, fill = Model)) +
      geom_boxplot(alpha = 0.7) +
      labs(
        title = "RMSE Distribution: Random Forest vs Quasi-Poisson",
        subtitle = "Lower RMSE indicates better performance",
        y = "Root Mean Square Error"
      ) +
      theme_minimal() +
      scale_fill_manual(values = c("Random Forest" = "#78B7C5", "Quasi-Poisson" = "#F21A00"))
    
    print(p_metrics)
  }
  
  # Residuals vs fitted plot for Quasi-Poisson
  if (nrow(predictions_qp_combined[set == "test"]) > 0) {
    qp_test_data <- predictions_qp_combined[set == "test"]
    qp_test_data[, residuals := obs_orig - pred]
    
    p_residuals <- ggplot(qp_test_data, aes(x = pred, y = residuals)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = TRUE, color = "blue") +
      labs(
        title = "Quasi-Poisson Model: Residuals vs Fitted Values",
        subtitle = "Test set predictions - should show random scatter around zero",
        x = "Fitted Values",
        y = "Residuals"
      ) +
      theme_minimal()
    
    print(p_residuals)
  }
}

#### SUMMARY AND RECOMMENDATIONS ####
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("=== SUMMARY AND RECOMMENDATIONS ===\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\nMODEL COMPARISON:\n")
cat("- Random Forest: Non-parametric, handles complex interactions, ensemble method\n")
cat("- Quasi-Poisson: Parametric, accounts for overdispersion, interpretable coefficients\n")

cat("\nQUASI-POISSON ADVANTAGES:\n")
cat("- Handles overdispersion in count data (dispersion parameter ≠ 1)\n")
cat("- More flexible than standard Poisson for overdispersed data\n")
cat("- Provides interpretable coefficients (log rate ratios)\n")
cat("- Well-established statistical foundation for count data\n")
cat("- Uses same mean-variance relationship as Poisson but allows different variance\n")

cat("\nPERFORMANCE SUMMARY:\n")
if (exists("presence_comparison")) {
  rf_better_acc <- presence_comparison$Mean_Accuracy[1] > presence_comparison$Mean_Accuracy[2]
  rf_better_auc <- presence_comparison$Mean_AUC[1] > presence_comparison$Mean_AUC[2]
  
  cat("Presence prediction:\n")
  cat("  - Better accuracy:", ifelse(rf_better_acc, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Better AUC:", ifelse(rf_better_auc, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("abundance_comparison") && nrow(abundance_comparison) > 0) {
  rf_better_rmse <- abundance_comparison$Mean_RMSE[1] < abundance_comparison$Mean_RMSE[2]
  rf_better_r2 <- abundance_comparison$Mean_R2[1] > abundance_comparison$Mean_R2[2]
  
  cat("Count prediction:\n")
  cat("  - Lower RMSE (better):", ifelse(rf_better_rmse, "Random Forest", "Quasi-Poisson"), "\n")
  cat("  - Higher R² (better):", ifelse(rf_better_r2, "Random Forest", "Quasi-Poisson"), "\n")
}

if (exists("qp_dispersion_summary")) {
  cat("Overdispersion:\n")
  cat("  - Mean dispersion parameter:", qp_dispersion_summary$Mean_Dispersion, "\n")
  if (qp_dispersion_summary$Mean_Dispersion > 1.2) {
    cat("  - Data shows overdispersion, Quasi-Poisson is appropriate\n")
  }
}

cat("\nRECOMMENDATION:\n")
cat("Choose Random Forest if: Complex non-linear relationships, prediction accuracy priority\n")
cat("Choose Quasi-Poisson if: Need interpretable coefficients, understand linear relationships\n")
cat("Consider Quasi-Poisson if: Count data with overdispersion, parametric approach preferred\n")









###### NICHE ANALYSIS USING RESULTS FROM MODEL ##############################################################################

### Niche analysis with species-specific variable selection
rm(list=ls())

# Load required libraries
library(dplyr)
library(ade4)
library(ggplot2)
library(patchwork)
library(data.table)

select <- dplyr::select

# Read data
dados_wide <- read.csv("Species_Filtered_VIF.csv",header=T,sep=",")
dados_tudo <-  read.csv("ABUND_DOMINANTT.csv",header=T,sep=",")

# SOLUTION: Remove existing species columns from dados_wide first, then merge
cat("Removing existing species columns from dados_wide to avoid conflicts...\n")
dados_wide_clean <- dados_wide %>%
  select(-c(Sanionia, Usnea, Colobanthus, Deschampsia))

# Now merge with the abundance data
dados_wide <- dados_wide_clean %>%
  left_join(dados_tudo %>% select(SITE_QUAD, Deschampsia, Sanionia, Usnea, Colobanthus), 
            by = c("SITE_FINAL" = "SITE_QUAD")) %>%
  mutate(LICHEN = Usnea,
         BRIO = Sanionia, 
         PLANT = Deschampsia + Colobanthus)

cat("Merge completed successfully\n")

# Prepare log variables
dados_wide <- dados_wide %>%
  mutate(LICHEN_LOG = log(LICHEN + 1),
         BRIO_LOG = log(BRIO + 1),
         PLANT_LOG = log(PLANT + 1)) %>%
  mutate_at(vars(LICHEN_LOG, BRIO_LOG, PLANT_LOG), ~ifelse(is.infinite(.), 0, .))

### ENVIRONMENTAL VARIABLES SELECTION ###
# Define all available environmental variables (expanded based on your species results)
all_env_vars <- c("WEXP", "DISTCOAST", "VDEP", "EFAH", "HURS_RNGE", "ALTITUDE", "WEFF", "Bio3_81", 
                  "ASP", "Easting", "Northing", "Bio15_81", "FLAC", "TRI", "SVF", "HURS_09",
                  "MPI", "LSF", "HillSh", "PR_01", "MBI", "PDIS_CLASS", "RSP", "TAS_08", "TPW", "MNCURV", "FCF", "MXCURV")

# Check if all environmental variables exist in the dataset
missing_vars <- setdiff(all_env_vars, names(dados_wide))
if(length(missing_vars) > 0) {
  cat("Missing environmental variables:", paste(missing_vars, collapse = ", "), "\n")
  # Continue with available variables instead of stopping
  all_env_vars <- intersect(all_env_vars, names(dados_wide))
}

### SPECIES-SPECIFIC VARIABLE IMPORTANCE ANALYSIS ###
cat("Using updated species-specific important variables...\n")

# Define species-specific important variables based on new requirements
species_important_vars <- list(
  # Usnea - updated variables
  Usnea = c("Easting", "Northing", "MXCURV", "WEFF", "WEXP", "ALTITUDE", "ASP", "VDEP", "DISTCOAST"),
  
  # Sanionia - updated variables  
  Sanionia = c("Easting", "Northing", "TRI", "SVF", "WEFF", "ALTITUDE", "ASP", "FLAC", "VDEP", "DISTCOAST"),
  
  # Deschampsia - updated variables
  Deschampsia = c("Easting", "Northing", "MNCURV", "MBI", "EFAH", "MPI", "SVF", "WEXP", "LSF", "ASP", "TPW", "FLAC", "DISTCOAST"),
  
  # Colobanthus - updated variables
  Colobanthus = c("ASP", "TPW", "FLAC", "VDEP", "DISTCOAST")
)

# Display species-specific important variables
cat("\n=== UPDATED SPECIES-SPECIFIC IMPORTANT VARIABLES ===\n")
for(species in names(species_important_vars)) {
  cat(sprintf("\n%s:\n", species))
  imp_vars <- species_important_vars[[species]]
  if(length(imp_vars) > 0) {
    for(i in 1:length(imp_vars)) {
      cat(sprintf("  %d. %s\n", i, imp_vars[i]))
    }
  } else {
    cat("  No variables specified\n")
  }
}

# Get union of all important variables
species_specific_vars <- unique(unlist(species_important_vars))

# Filter to only include variables that exist in your dataset
available_vars <- intersect(species_specific_vars, all_env_vars)
missing_vars <- setdiff(species_specific_vars, all_env_vars)

if(length(missing_vars) > 0) {
  cat("\nWarning: Some specified variables not available in dataset:\n")
  cat(paste(missing_vars, collapse = ", "), "\n")
  cat("These will be excluded from analysis.\n")
}

species_specific_vars <- available_vars

cat("\n=== FINAL SELECTED VARIABLES ===\n")
cat("Species-specific important variables:", paste(species_specific_vars, collapse = ", "), "\n")
cat("Number of variables selected:", length(species_specific_vars), "\n")

# Use species-specific variables for analysis
env_vars <- species_specific_vars

# Create environmental data frame with species and selected environmental variables
env_mod <- dados_wide %>%
  dplyr::select(Usnea, Sanionia, Deschampsia, Colobanthus, all_of(env_vars))

# Rename environmental variables for better readability (only if they exist)
rename_mapping <- c(
  "EFAH" = "EFFEC_AIR_FLW_HGHT",
  "VDEP" = "VALDEPTH", 
  "WEXP" = "WINDEXPO",
  "DISTCOAST" = "DISTCOAST",
  "HURS_RNGE" = "REL_HUMID",
  "ALTITUDE" = "ALTITUDE",
  "WEFF" = "WINDEFFECT",
  "Bio3_81" = "ISOTHERMALITY"
)

# Apply renaming only for variables that exist
for(old_name in names(rename_mapping)) {
  if(old_name %in% names(env_mod)) {
    names(env_mod)[names(env_mod) == old_name] <- rename_mapping[old_name]
  }
}

cat("Environmental variables being used in analysis:\n")
env_var_names <- names(env_mod)[5:ncol(env_mod)]
cat(paste(env_var_names, collapse = ", "), "\n")

#### PCA - SELECTED ENVIRONMENTAL VARIABLES ONLY ####
# Extract only the selected environmental variables for PCA
env_vars_only <- env_mod %>%
  select(-c(Usnea, Sanionia, Deschampsia, Colobanthus))

# Check data types and convert to numeric if needed
cat("Checking data types of environmental variables:\n")
for(var in names(env_vars_only)) {
  cat(sprintf("%s: %s\n", var, class(env_vars_only[[var]])))
}

# Convert all environmental variables to numeric
cat("Converting variables to numeric...\n")
env_vars_numeric <- env_vars_only %>%
  mutate_all(~ as.numeric(as.character(.)))

# Check for any variables that couldn't be converted
numeric_check <- sapply(env_vars_numeric, function(x) all(is.na(x)))
if(any(numeric_check)) {
  cat("Warning: These variables could not be converted to numeric:\n")
  cat(paste(names(numeric_check)[numeric_check], collapse = ", "), "\n")
  # Remove problematic variables
  env_vars_numeric <- env_vars_numeric[, !numeric_check]
}

# Scale the environmental variables
env_mod_scale <- scale(env_vars_numeric)

# Check for any missing values
if(any(is.na(env_mod_scale))) {
  cat("Warning: NA values found in environmental data\n")
  complete_rows <- complete.cases(env_mod_scale)
  env_mod_scale <- env_mod_scale[complete_rows, ]
  env_mod <- env_mod[complete_rows, ]
  dados_wide <- dados_wide[complete_rows, ]
  cat("Removed", sum(!complete_rows), "rows with missing environmental data\n")
}

# Perform PCA on selected environmental variables
cat("Performing PCA on species-specific selected environmental variables...\n")
pca.env_mod <- dudi.pca(env_mod_scale, scannf=FALSE)

# Display PCA results
cat("PCA completed. Variables included:\n")
cat(paste(colnames(env_mod_scale), collapse = ", "), "\n")
cat("Variance explained by first two axes:", 
    round(100 * sum(pca.env_mod$eig[1:2])/sum(pca.env_mod$eig), 2), "%\n")

#### OMI ANALYSIS ####
# Create species matrix
spe <- dados_wide[, c("Usnea", "Sanionia", "Deschampsia", "Colobanthus")]

# Handle NA values in species data
if(any(is.na(spe))) {
  cat("Replacing NA values with 0 in species data\n")
  spe[is.na(spe)] <- 0
}

# Ensure we have the same number of rows for PCA and species data
if(nrow(spe) != nrow(env_mod_scale)) {
  spe <- spe[complete.cases(env_mod_scale), ]
}

cat("Performing OMI analysis...\n")
omi1 <- niche(pca.env_mod, spe, scannf=FALSE)

# Display variance explained by individual OMI axes
var_explained <- 100 * omi1$eig/sum(omi1$eig)
cat("\n=== OMI VARIANCE EXPLAINED ===\n")
cat("OMI axis 1 variance explained:", round(var_explained[1], 1), "%\n")
cat("OMI axis 2 variance explained:", round(var_explained[2], 1), "%\n")
cat("Total variance explained by first two axes:", round(sum(var_explained[1:2]), 1), "%\n")

# Statistical tests
cat("\n=== OMI STATISTICAL TESTS ===\n")
rtest_result <- rtest(omi1, nrepet=1000)
print(rtest_result)

niche_params <- niche.param(omi1)
print(niche_params)

### Export OMI scores with coordinates ###
pontos <- omi1$ls
pontos$SITE_FINAL <- dados_wide$SITE_FINAL[complete.cases(env_mod_scale)]

# Add coordinates if they exist
if(all(c("Easting", "Northing") %in% names(dados_wide))) {
  coord_data <- dados_wide %>%
    select(SITE_FINAL, Easting, Northing)
  pontos <- merge(pontos, coord_data, by="SITE_FINAL", all.x=TRUE)
  cat("Coordinates added to output\n")
} else {
  cat("Coordinate columns not found\n")
}

# Save results with species-specific label
output_filename <- paste0("OMI_coord_updated_species_specific_", length(env_vars), "vars.csv")
write.csv(pontos, file = output_filename, row.names = FALSE)
cat("Results saved to:", output_filename, "\n")

#### Correlations between environmental variables and OMI axes ####
cat("\n=== CORRELATIONS WITH OMI AXES ===\n")
cat("Correlations with OMI axis 1:\n")
cor_omi1 <- cor(as.matrix(omi1$ls[,1]), env_mod_scale)
print(round(cor_omi1, 3))

cat("\nCorrelations with OMI axis 2:\n")
cor_omi2 <- cor(as.matrix(omi1$ls[,2]), env_mod_scale)
print(round(cor_omi2, 3))

## Niche analysis
par(mfrow=c(1,2))
OMI1 <- sco.distri(omi1$ls[,1], spe, clab=0.7, y.rank = TRUE)
OMI2 <- sco.distri(omi1$ls[,2], spe, clab=0.7)

### Prepare data for plotting ###
tabela <- data.table(rbind(OMI1, OMI2))
tabela <- tabela %>%
  mutate(OMI = rep(c(1,1,1,1,2,2,2,2), 1),
         VEG = rep(c("Usnea", "Deschampsia", "Colobanthus", "Sanionia"), 2))

tabela_B <- data.table(OMI1)
tabela_B$VEG <- c("Usnea", "Deschampsia", "Colobanthus", "Sanionia")
tabela_B <- tabela_B %>%
  rename(OMI1_mean = mean, OMI1_var = var)

tabela_B <- cbind(tabela_B, OMI2)
tabela_B <- tabela_B %>%
  rename(OMI2_mean = mean, OMI2_var = var)

#### PREPARE SITE CLASSIFICATION DATA FOR PLOTTING ####

# Get the species data aligned with the environmental data (after removing incomplete cases)
spe_for_plotting <- spe[complete.cases(env_mod_scale), ]

# Create site classification based on species presence/dominance
site_classification <- data.frame(
  site_id = 1:nrow(spe_for_plotting),
  CS1 = omi1$ls$CS1,
  CS2 = omi1$ls$CS2,
  PC1 = pca.env_mod$li$Axis1,
  PC2 = pca.env_mod$li$Axis2,
  Usnea = spe_for_plotting$Usnea,
  Sanionia = spe_for_plotting$Sanionia,
  Deschampsia = spe_for_plotting$Deschampsia,
  Colobanthus = spe_for_plotting$Colobanthus
)

# Determine dominant species for each site
site_classification$dominant_species <- apply(site_classification[, c("Usnea", "Sanionia", "Deschampsia", "Colobanthus")], 1, function(x) {
  if(all(x == 0)) return("No_species")
  max_abundance <- max(x)
  dominant <- names(x)[which(x == max_abundance)]
  if(length(dominant) > 1) return("Mixed") # Tie for dominance
  return(dominant[1])
})

# Create presence-based coloring
site_classification$species_present <- apply(site_classification[, c("Usnea", "Sanionia", "Deschampsia", "Colobanthus")], 1, function(x) {
  present_species <- names(x)[x > 0]
  if(length(present_species) == 0) return("No_species")
  if(length(present_species) == 1) return(present_species[1])
  if(length(present_species) >= 2) return("Multiple_species")
})

# Create presence indicators for each species (for alternative layered approach)
site_classification$Usnea_present <- site_classification$Usnea > 0
site_classification$Sanionia_present <- site_classification$Sanionia > 0
site_classification$Deschampsia_present <- site_classification$Deschampsia > 0
site_classification$Colobanthus_present <- site_classification$Colobanthus > 0

# Print site classification summary
cat("\n=== SITE CLASSIFICATION SUMMARY ===\n")
cat("Dominant species classification:\n")
print(table(site_classification$dominant_species))
cat("\nSpecies presence classification:\n")
print(table(site_classification$species_present))

#### PLOTTING SECTION ####

#### DEFINE COLORS FOR INDIVIDUAL SPECIES ####
cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00",
           "Colobanthus" = "#78B7C5")

# Extended color scheme for site classifications
cores_extended <- c("Usnea" = "#00A08A",
                    "Sanionia" = "#EBCC2A", 
                    "Deschampsia" = "#F21A00",
                    "Colobanthus" = "#78B7C5",
                    "Mixed" = "#999999",
                    "Multiple_species" = "#2d2d2d",
                    "No_species" = "#CCCCCC")

#### NEW LAYOUT PANELS ####

## PANEL A - PCA CIRCLE WITH ONLY AXIS 1 AND AXIS 2 ARROWS
panel_A <- ggplot() +
  # Add a larger circle
  annotate("path",
           x = cos(seq(0, 2*pi, length.out = 100)) * 1.5,
           y = sin(seq(0, 2*pi, length.out = 100)) * 1.5,
           color = "black", size = 1) +
  # Axis 1 arrow (horizontal, pointing right)
  geom_segment(aes(x = 0, xend = 1.5, y = 0, yend = 0), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               colour = "black", size = 1) +
  # Axis 2 arrow (vertical, pointing up)
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1.5), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               colour = "black", size = 1) +
  # Axis labels
  geom_text(aes(x = 1.7, y = 0), label = "Axis1", size = 4, fontface = "bold") +
  geom_text(aes(x = 0, y = 1.7), label = "Axis2", size = 4, fontface = "bold") +
  coord_fixed() +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  xlim(-2.2, 2.2) +
  ylim(-2.2, 2.2) +
  ggtitle("A")

## PANEL B - ABUNDANCE CURVES ALONG OMI1
# Prepare data for abundance curves
ind_env_mod <- omi1$ls
ind_env_mod$Usnea <- spe$Usnea
ind_env_mod$Deschampsia <- spe$Deschampsia
ind_env_mod$Colobanthus <- spe$Colobanthus
ind_env_mod$Sanionia <- spe$Sanionia

# Create loess curves for OMI1
g1 <- ggplot(ind_env_mod) + 
  stat_smooth(aes(x = CS1, y = Usnea), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS1, y = Deschampsia), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS1, y = Colobanthus), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS1, y = Sanionia), method = "loess", se = FALSE)

gg1 <- ggplot_build(g1)
df2 <- data.frame(x = gg1$data[[1]]$x,
                  Usnea = gg1$data[[1]]$y,
                  Deschampsia = gg1$data[[2]]$y,
                  Colobanthus = gg1$data[[3]]$y,
                  Sanionia = gg1$data[[4]]$y)

# Color scheme: Bryophytes (blue), Lichens (red), Vascular Plants (dark)
panel_B <- ggplot(df2) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  geom_area(aes(x = x, y = Usnea), fill = "#00A08A", alpha = 0.7) +     # Usnea
  geom_area(aes(x = x, y = Sanionia), fill = "#EBCC2A", alpha = 0.7) +  # Sanionia
  geom_area(aes(x = x, y = Deschampsia), fill = "#F21A00", alpha = 0.7) + # Deschampsia
  geom_area(aes(x = x, y = Colobanthus), fill = "#78B7C5", alpha = 0.7) + # Colobanthus
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  xlab("OMI1") +
  ylab("Abundance") +
  ggtitle("B")

## PANEL C - NICHE BREADTH OMI1
panel_C <- tabela_B %>%
  ggplot(aes(x = OMI1_mean, y = VEG, color = VEG)) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  geom_point(size = 3) +
  geom_linerange(aes(xmin = OMI1_mean - OMI1_var, xmax = OMI1_mean + OMI1_var), size = 1.2) +
  scale_color_manual(values = cores) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12, face = "bold")) +
  xlab("OMI1") +
  ggtitle("C")

## PANEL D - PCA BIPLOT WITHOUT SITE POINTS (ARROWS ONLY)
env_pca_vars <- pca.env_mod$c1 %>%
  tibble::rownames_to_column()

# Calculate positions respecting original loadings but with good visibility
env_pca_vars <- env_pca_vars %>%
  mutate(
    arrow_x = CS1 * 1.7,  
    arrow_y = CS2 * 1.7,
    text_x = CS1 * 1.8,   
    text_y = CS2 * 1.8,
    hjust_val = case_when(
      CS1 > 0.2 ~ 0,      
      CS1 < -0.2 ~ 1,     
      TRUE ~ 0.5          
    ),
    vjust_val = case_when(
      CS2 > 0.2 ~ 0,      
      CS2 < -0.2 ~ 1,     
      TRUE ~ 0.5          
    )
  )

panel_D <- ggplot() +
  geom_vline(xintercept = seq(-1, 1, 0.5), color = "grey90", size = 0.5) +
  geom_hline(yintercept = seq(-1, 1, 0.5), color = "grey90", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.8) +
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  # Environmental arrows only - no site points
  geom_segment(data = env_pca_vars, aes(x = 0, xend = arrow_x, y = 0, yend = arrow_y), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               colour = "black", size = 0.8) +
  geom_text(data = env_pca_vars, aes(x = text_x, y = text_y, label = rowname, 
                                     hjust = hjust_val, vjust = vjust_val), 
            size = 3, fontface = "bold", color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0, size = 14, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-1.3, 1.3), breaks = seq(-1, 1, 0.5)) +
  scale_x_continuous(limits = c(-1.3, 1.3), breaks = seq(-1, 1, 0.5)) +
  xlab(paste0("PC1 (", round(100 * pca.env_mod$eig[1]/sum(pca.env_mod$eig), 1), "%)")) +
  ylab(paste0("PC2 (", round(100 * pca.env_mod$eig[2]/sum(pca.env_mod$eig), 1), "%)")) +
  ggtitle("D")

## PANEL E - OMI ORDINATION
# Create data frame for ellipses
ellipse_data_base <- data.frame(
  CS1 = rep(omi1$ls$CS1, 4),
  CS2 = rep(omi1$ls$CS2, 4),
  Species = rep(c("Usnea", "Sanionia", "Deschampsia", "Colobanthus"), each = nrow(omi1$ls)),
  Abundance = c(spe$Usnea, spe$Sanionia, spe$Deschampsia, spe$Colobanthus)
)

# Only keep sites where species is present
ellipse_data <- ellipse_data_base[ellipse_data_base$Abundance > 0, ]

panel_E <- omi1$ls %>%
  ggplot(aes(x = CS1, y = CS2)) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  stat_ellipse(data = ellipse_data, aes(x = CS1, y = CS2, color = Species), 
               type = "norm", level = 0.68, size = 1, alpha = 0.7) +
  geom_point(data = omi1$ls, aes(x = CS1, y = CS2), color = "black", size = 2, alpha = 0.4) +
  geom_point(data = omi1$li, aes(x = Axis1, y = Axis2), 
             color = c("#00A08A", "#EBCC2A", "#F21A00", "#78B7C5"), size = 4, shape = 19) +
  scale_color_manual(values = cores) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-4, 7), breaks = c(-4, -2, 0, 2, 4, 6)) +
  scale_x_continuous(limits = c(-4, 7), breaks = c(-2, 0, 2, 4, 6)) +
  xlab("OMI1") +
  ylab("OMI2") +
  ggtitle("E")

## PANEL F - NICHE BREADTH OMI2
panel_F <- tabela_B %>%
  ggplot(aes(x = OMI2_mean, y = VEG, color = VEG)) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  geom_point(size = 3) +
  geom_linerange(aes(xmin = OMI2_mean - OMI2_var, xmax = OMI2_mean + OMI2_var), size = 1.2) +
  scale_color_manual(values = cores) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        legend.position = "right",
        axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12, face = "bold")) +
  scale_x_continuous(limits = c(-4, 7), breaks = c(-2, 0, 2, 4, 6)) +
  xlab("OMI2") +
  ggtitle("F") +
  coord_flip()

## PANEL G - ABUNDANCE CURVES ALONG OMI2
# Create loess curves for OMI2
g2 <- ggplot(ind_env_mod) + 
  stat_smooth(aes(x = CS2, y = Usnea), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS2, y = Deschampsia), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS2, y = Colobanthus), method = "loess", se = FALSE) +
  stat_smooth(aes(x = CS2, y = Sanionia), method = "loess", se = FALSE)

gg2 <- ggplot_build(g2)
df3 <- data.frame(x = gg2$data[[1]]$x,
                  Usnea = gg2$data[[1]]$y,
                  Deschampsia = gg2$data[[2]]$y,
                  Colobanthus = gg2$data[[3]]$y,
                  Sanionia = gg2$data[[4]]$y)

panel_G <- ggplot(df3) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_area(aes(x = x, y = Usnea), fill = "#00A08A", alpha = 0.7) +     # Usnea
  geom_area(aes(x = x, y = Sanionia), fill = "#EBCC2A", alpha = 0.7) +  # Sanionia
  geom_area(aes(x = x, y = Deschampsia), fill = "#F21A00", alpha = 0.7) + # Deschampsia
  geom_area(aes(x = x, y = Colobanthus), fill = "#78B7C5", alpha = 0.7) + # Colobanthus
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12, face = "bold")) +
  scale_x_continuous(limits = c(-4, 7), breaks = c(-4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  xlab("OMI2") +
  ggtitle("G") +
  coord_flip()

### CREATE FINAL LAYOUT ###

# Layout design to match the reference image exactly
layout_design <- "
ABBCD
ABBCD
EFFFG
EFFFG
"

# Apply consistent theme to all panels
bold_theme <- theme(
  plot.title = element_text(hjust = 0, size = 16, face = "bold"),
  axis.text = element_text(size = 14, face = "bold"),
  axis.title = element_text(size = 16, face = "bold"), 
  axis.line = element_line(color = "black", size = 1.5)
)

# Apply theme to panels (except panel A which needs special handling)
panel_A <- panel_A + theme(plot.title = element_text(hjust = 0, size = 16, face = "bold"))
panel_B <- panel_B + bold_theme
panel_C <- panel_C + bold_theme  
panel_D <- panel_D + bold_theme
panel_E <- panel_E + bold_theme
panel_F <- panel_F + bold_theme
panel_G <- panel_G + bold_theme

# Create final plot with new layout
final_plot <- panel_A + panel_B + panel_C + panel_D + panel_E + panel_F + panel_G +
  plot_layout(design = layout_design) & theme(legend.position = "none")

# Remove legend references and use original color scheme
legend_data <- data.frame(
  Species = c("Usnea", "Sanionia", "Deschampsia", "Colobanthus"),
  Color = c("#00A08A", "#EBCC2A", "#F21A00", "#78B7C5")
)

# Print the final plot
print(final_plot)






########################### Partial responses to the top variables ################################################################


# Load required libraries
library(randomForest)
library(quantregForest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Read data
dados_wide <- read.csv("Species_Filtered_VIF.csv", header=T, sep=",")

#### DATA PREPARATION ####

# Set the proportion for training and test data
train_prop <- 0.8
set.seed(42)

train_data <- dados_wide %>% 
  group_by(ALT_CLASS) %>%
  sample_frac(train_prop) %>%
  ungroup()

# Create presence/absence variables for each species
train_data$Usnea_Presence <- as.factor(ifelse(train_data$Usnea > 0, 1, 0))
train_data$Sanionia_Presence <- as.factor(ifelse(train_data$Sanionia > 0, 1, 0))
train_data$Deschampsia_Presence <- as.factor(ifelse(train_data$Deschampsia > 0, 1, 0))
train_data$Colobanthus_Presence <- as.factor(ifelse(train_data$Colobanthus > 0, 1, 0))

#### DEFINE VARIABLES FOR HURDLE MODELS ####

# Species-specific top variables from hurdle model analysis
variables_list <- list(
  "Usnea" = c("ALTITUDE", "WEXP", "Northing", "Easting", "HURS_09", "DISTCOAST", "VDEP"),
  "Sanionia" = c("ALTITUDE", "VDEP", "WEFF", "FLAC", "Easting", "WEXP", "Northing", "DISTCOAST"),
  "Deschampsia" = c("DISTCOAST", "TPW", "WEXP", "Northing", "Easting", "LSF", "VDEP"),
  "Colobanthus" = c("FLAC", "VDEP", "DISTCOAST", "FCF", "Northing", "LSF")
)

# Get all unique variables for plotting
all_variables_to_plot <- unique(unlist(variables_list))

# Variable label mapping
variable_labels <- c(
  "WEXP" = "Wind Exposition",
  "HURS_09" = "September humidity",
  "DISTCOAST" = "Distance to coast",
  "VDEP" = "Valley Depth",
  "WEFF" = "Wind Effect",
  "FLAC" = "Flow Accumulation",
  "TPW" = "Topographic Wetness",
  "FCF" = "Frost Change Frequency",
  "LSF" = "LS Factor",
  "ALTITUDE" = "Altitude",
  "Northing" = "Northing",
  "Easting" = "Easting",
)

# Species names
species_names <- c("Usnea", "Sanionia", "Deschampsia", "Colobanthus")

#### HURDLE MODEL TRAINING ####

cat("Training hurdle models...\n")

# Storage for models
hurdle_models <- list()

for(species in species_names) {
  cat("Training", species, "hurdle model...\n")
  
  # Get presence and abundance column names
  presence_col <- paste0(species, "_Presence")
  abundance_col <- species
  
  # Get species-specific variables
  species_variables <- variables_list[[species]]
  
  # STAGE 1: Presence/Absence Model (Random Forest Classifier)
  presence_cols <- c(presence_col, species_variables)
  presence_data <- train_data[, presence_cols, drop = FALSE]
  names(presence_data)[1] <- "Presence"
  
  presence_model <- randomForest(
    Presence ~ ., 
    data = presence_data,
    ntree = 500,
    importance = TRUE
  )
  
  # STAGE 2: Abundance Model (Quantile Regression Forest)
  present_data <- train_data[train_data[[presence_col]] == "1", ]
  abundance_cols <- c(abundance_col, species_variables)
  abundance_data <- present_data[, abundance_cols, drop = FALSE]
  names(abundance_data)[1] <- "Abundance"
  
  abundance_model <- NULL
  if(nrow(abundance_data) > 10) {
    abundance_model <- quantregForest(
      x = abundance_data[, species_variables, drop = FALSE],
      y = abundance_data$Abundance,
      ntree = 500
    )
  }
  
  hurdle_models[[species]] <- list(
    presence = presence_model,
    abundance = abundance_model,
    presence_data = presence_data,
    abundance_data = abundance_data,
    species_variables = species_variables
  )
}

#### DEFINE COLORS ####

cores <- c("Usnea" = "#00A08A",
           "Sanionia" = "#EBCC2A", 
           "Deschampsia" = "#F21A00",
           "Colobanthus" = "#78B7C5")

#### IMPROVED HURDLE MODEL PARTIAL DEPENDENCE FUNCTIONS ####

create_hurdle_partial_data <- function(hurdle_model_info, variable, grid_resolution = 50) {
  presence_model <- hurdle_model_info$presence
  abundance_model <- hurdle_model_info$abundance
  presence_data <- hurdle_model_info$presence_data
  species_variables <- hurdle_model_info$species_variables
  
  # Use the full range from training data (including zeros)
  var_range <- range(train_data[[variable]], na.rm = TRUE)
  
  # Check for non-finite values
  if(!is.finite(var_range[1]) || !is.finite(var_range[2])) {
    warning(paste("Variable", variable, "has non-finite values. Skipping."))
    return(NULL)
  }
  
  var_grid <- seq(var_range[1], var_range[2], length.out = grid_resolution)
  predictions <- numeric(length(var_grid))
  
  for(i in seq_along(var_grid)) {
    # Create modified dataset for presence prediction
    temp_presence_data <- presence_data
    temp_presence_data[[variable]] <- var_grid[i]
    
    # Stage 1: Predict presence probability
    presence_prob <- predict(presence_model, newdata = temp_presence_data, type = "prob")
    if(ncol(presence_prob) >= 2) {
      prob_present <- mean(presence_prob[, 2], na.rm = TRUE)
    } else {
      prob_present <- mean(presence_prob[, 1], na.rm = TRUE)
    }
    
    # Stage 2: Predict abundance given presence
    abundance_given_presence <- 0
    if(!is.null(abundance_model) && !is.null(hurdle_model_info$abundance_data)) {
      tryCatch({
        # Create modified dataset for abundance prediction
        temp_abundance_data <- hurdle_model_info$abundance_data
        temp_abundance_data[[variable]] <- var_grid[i]
        
        # Predict abundance using quantile regression (median)
        abundance_pred <- predict(abundance_model, 
                                  newdata = temp_abundance_data[, species_variables, drop = FALSE], 
                                  what = 0.5)
        abundance_given_presence <- mean(abundance_pred, na.rm = TRUE)
      }, error = function(e) {
        # If abundance prediction fails, use mean abundance from training data
        abundance_given_presence <- mean(hurdle_model_info$abundance_data$Abundance, na.rm = TRUE)
      })
    }
    
    # Combine: Expected abundance = P(presence) × E(abundance|presence)
    predictions[i] <- prob_present * abundance_given_presence
  }
  
  result <- data.frame(
    x = var_grid,
    yhat = predictions
  )
  names(result)[1] <- variable
  
  return(result)
}

create_combined_abundance_data <- function(variable_name) {
  combined_data <- data.frame()
  
  for(species in species_names) {
    model_info <- hurdle_models[[species]]
    species_variables <- variables_list[[species]]
    
    if(variable_name %in% species_variables && variable_name %in% names(train_data)) {
      tryCatch({
        partial_data <- create_hurdle_partial_data(model_info, variable_name)
        if(!is.null(partial_data)) {  # Check if partial_data is not NULL
          temp_data <- data.frame(
            x = partial_data[[variable_name]],
            y = partial_data$yhat,
            species = species
          )
          combined_data <- rbind(combined_data, temp_data)
        }
      }, error = function(e) {
        cat("Error creating abundance data for", species, variable_name, ":", e$message, "\n")
      })
    }
  }
  
  return(combined_data)
}

create_combined_abundance_plot <- function(variable_name) {
  plot_data <- create_combined_abundance_data(variable_name)
  
  if(nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = x, y = y, color = species)) +
    geom_line(linewidth = 1.5) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 11, face = "bold", margin = ggplot2::margin(t = 15)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "none",
      plot.title = element_blank()
    ) +
    scale_color_manual(values = cores) +
    labs(
      x = variable_labels[[variable_name]] %||% variable_name,
      y = NULL
    ) +
    ylim(0, max(plot_data$y) * 1.1)  # Dynamic y-axis based on data
  
  return(p)
}

#### GENERATE ABUNDANCE PLOTS ####

cat("Generating hurdle model partial dependence plots...\n")
combined_plots_abundance <- list()

for(var in all_variables_to_plot) {
  cat("Creating abundance plot for:", var, "\n")
  combined_plots_abundance[[var]] <- create_combined_abundance_plot(var)
}

combined_plots_abundance <- combined_plots_abundance[!sapply(combined_plots_abundance, is.null)]
valid_abundance_plots <- sapply(combined_plots_abundance, function(x) inherits(x, "ggplot"))
combined_plots_abundance <- combined_plots_abundance[valid_abundance_plots]

cat("Displaying abundance model plots...\n")
if(length(combined_plots_abundance) > 0) {
  n_plots <- length(combined_plots_abundance)
  ncols <- min(3, n_plots)
  nrows <- ceiling(n_plots / ncols)
  
  tryCatch({
    grid.arrange(grobs = combined_plots_abundance, 
                 ncol = ncols, 
                 nrow = nrows,
                 left = textGrob("Expected abundance", rot = 90, vjust = 1, 
                                 gp = gpar(fontsize = 14, fontface = "bold")))
  }, error = function(e) {
    cat("Error displaying abundance plots:", e$message, "\n")
    do.call(grid.arrange, c(combined_plots_abundance, ncol = ncols, nrow = nrows))
  })
}

cat("\n=== HURDLE MODEL PARTIAL DEPENDENCE PLOTS COMPLETED ===\n")
cat("Expected abundance plots generated:", length(combined_plots_abundance), "\n")
cat("These plots now show the full ecological response including zero abundance areas.\n")

