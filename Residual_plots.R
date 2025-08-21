# Load the saved results
all_results_clf <- readRDS("hurdle_rf_classifier_results_USNEA_raw.rds")
all_results_qrf <- readRDS("hurdle_qrf_abundance_results_USNEA_raw.rds")
varimp_list_clf <- readRDS("hurdle_varimp_rf_presence_USNEA_raw.rds")
varimp_list_qrf <- readRDS("hurdle_varimp_qrf_abundance_USNEA_raw.rds")
predictions_all <- readRDS("hurdle_predictions_all_USNEA_raw.rds")
training_data_list <- readRDS("hurdle_training_data_list_USNEA_raw.rds")
test_data_list <- readRDS("hurdle_test_data_list_USNEA_raw.rds")

# ===== CLASSIFIER METRICS SUMMARY =====
print("=== RF CLASSIFIER PERFORMANCE SUMMARY ===")
metrics_clf <- rbindlist(all_results_clf)
metrics_clf_summary <- metrics_clf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_clf_summary$Model <- "RF_Presence"
print(metrics_clf_summary)
clipr::write_clip(metrics_clf_summary)

# ===== QRF ABUNDANCE METRICS SUMMARY =====
print("=== QRF ABUNDANCE PERFORMANCE SUMMARY ===")
metrics_qrf <- rbindlist(all_results_qrf)
metrics_qrf_summary <- metrics_qrf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_qrf_summary$Model <- "QRF_Abundance_sqrt"
metrics_qrf_summary$Transformation <- "Square_Root"
print(metrics_qrf_summary)
clipr::write_clip(metrics_qrf_summary)

# ===== VARIABLE IMPORTANCE SUMMARY =====
print("=== VARIABLE IMPORTANCE SUMMARY ===")

# Combine classifier and QRF importance
varimp_clf <- rbindlist(varimp_list_clf)
varimp_qrf <- rbindlist(varimp_list_qrf[!sapply(varimp_list_qrf, is.null)])
varimp_all <- rbind(varimp_clf, varimp_qrf)

# Summarize by model type
varimp_summary <- varimp_all %>%
  group_by(Variable, Model) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 1), .groups = "drop") %>%
  arrange(Model, desc(MedianImportance)) %>%
  group_by(Model) %>%
  mutate(Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 1),
         Cumulative_pct = round(cumsum(100 * MedianImportance / sum(MedianImportance)), 1)) %>%
  ungroup()

print(varimp_summary)
clipr::write_clip(varimp_summary)

# ===== PREDICTION ANALYSIS (TEST DATA ONLY) =====
# Get test predictions only
predictions_test <- rbindlist(predictions_all[!sapply(predictions_all, is.null)])[set == "test"]

if (nrow(predictions_test) > 0) {
  print("=== PREDICTION ANALYSIS ===")
  
  # ===== COVERAGE PLOT BY ABUNDANCE BINS =====
  # Use first test dataset for binning analysis
  test_data_sample <- test_data_list[[1]]
  present_idx <- which(test_data_sample$Presence == "1")
  test_present_sample <- test_data_sample[present_idx, ]
  
  if (nrow(test_present_sample) > 10) {
    # Create bins for coverage analysis
    breaks <- seq(0, max(test_present_sample$Usnea, na.rm = TRUE), length.out = 6)
    if (length(breaks) < 2) breaks <- c(0, max(test_present_sample$Usnea, na.rm = TRUE))
    
    obs_binned <- cut(test_present_sample$Usnea, breaks = breaks, include.lowest = TRUE)
    
    # Calculate coverage using predictions from first iteration
    pred_sample <- predictions_test[iteration == 1]
    if (nrow(pred_sample) > 0) {
      # Calculate prediction intervals (assuming spread represents interval width)
      pred_lower <- pred_sample$pred - pred_sample$spread/2
      pred_upper <- pred_sample$pred + pred_sample$spread/2
      
      coverage_by_bin <- tapply(pred_sample$obs >= pred_lower & pred_sample$obs <= pred_upper, 
                                obs_binned[1:nrow(pred_sample)], mean, na.rm = TRUE) * 100
      
      # Create bin labels
      bin_labels <- paste0("[", round(breaks[-length(breaks)], 1), ",", round(breaks[-1], 1), ")")
      bin_labels[length(bin_labels)] <- paste0("[", round(breaks[length(breaks)-1], 1), ",", 
                                               round(breaks[length(breaks)], 1), "]")
      
      # Coverage plot with larger labels
      par(cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
      bp <- barplot(coverage_by_bin, col = "#00A08A", ylim = c(0, 100),
                    main = "95% Prediction Interval Coverage by Abundance Bins",
                    xlab = "Observed Usnea (Binned)", ylab = "Coverage (%)",
                    names.arg = "", xaxt = "n")
      text(bp, par("usr")[3] - 5, labels = bin_labels, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 1.1)
      abline(h = 95, col = "red", lty = 2, lwd = 2)
      # Reset par
      par(cex.main = 1, cex.lab = 1, cex.axis = 1)
    }
  }
  
  # ===== UNCERTAINTY VS ABUNDANCE PLOT =====
  spread_plot <- ggplot(predictions_test, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.3, color = "#00A08A") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(x = "Observed Usnea", y = "Prediction Interval Width",
         title = "Prediction Uncertainty vs Observed Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(spread_plot)
  
  # ===== OBSERVED VS PREDICTED PLOT =====
  obs_pred <- predictions_test %>%
    group_by(obs) %>%
    summarize(pred = median(pred, na.rm = TRUE), .groups = "drop")
  
  spearman_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "spearman", use = "complete.obs")^2
  pearson_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "pearson", use = "complete.obs")^2
  
  r2_label <- paste0("Spearman R² = ", round(spearman_r2, 3), "\nPearson R² = ", round(pearson_r2, 3))
  
  obs_vs_pred_plot <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.6, color = "#00A08A") +
    geom_smooth(method = "lm", se = FALSE, color = "#00A08A") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = -Inf, label = r2_label,
             hjust = 1.1, vjust = -0.5, size = 5, color = "blue") +
    labs(x = "Observed Usnea", y = "Predicted Usnea",
         title = "Observed vs Predicted Usnea (Test Data)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(obs_vs_pred_plot)
  
  # ===== RESIDUALS PLOT =====
  obs_pred$residuals <- obs_pred$obs - obs_pred$pred
  
  residuals_plot <- ggplot(obs_pred, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.6, color = "#00A08A") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "#00A06A") +
    labs(x = "Predicted Usnea", y = "Residuals (Observed - Predicted)",
         title = "Residuals vs Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(residuals_plot)
}



# Load the saved results
all_results_clf <- readRDS("hurdle_rf_classifier_results_Sanionia.rds")
all_results_qrf <- readRDS("hurdle_qrf_abundance_results_Sanionia.rds")
varimp_list_clf <- readRDS("hurdle_varimp_rf_presence_Sanionia.rds")
varimp_list_qrf <- readRDS("hurdle_varimp_qrf_abundance_Sanionia.rds")
predictions_all <- readRDS("hurdle_predictions_all_Sanionia.rds")
training_data_list <- readRDS("hurdle_training_data_list_Sanionia.rds")
test_data_list <- readRDS("hurdle_test_data_list_Sanionia.rds")

# ===== CLASSIFIER METRICS SUMMARY =====
print("=== RF CLASSIFIER PERFORMANCE SUMMARY ===")
metrics_clf <- rbindlist(all_results_clf)
metrics_clf_summary <- metrics_clf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_clf_summary$Model <- "RF_Presence"
print(metrics_clf_summary)
clipr::write_clip(metrics_clf_summary)

# ===== QRF ABUNDANCE METRICS SUMMARY =====
print("=== QRF ABUNDANCE PERFORMANCE SUMMARY ===")
metrics_qrf <- rbindlist(all_results_qrf)
metrics_qrf_summary <- metrics_qrf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_qrf_summary$Model <- "QRF_Abundance_sqrt"
metrics_qrf_summary$Transformation <- "Raw"
print(metrics_qrf_summary)
clipr::write_clip(metrics_qrf_summary)

# ===== VARIABLE IMPORTANCE SUMMARY =====
print("=== VARIABLE IMPORTANCE SUMMARY ===")

# Combine classifier and QRF importance
varimp_clf <- rbindlist(varimp_list_clf)
varimp_qrf <- rbindlist(varimp_list_qrf[!sapply(varimp_list_qrf, is.null)])
varimp_all <- rbind(varimp_clf, varimp_qrf)

# Summarize by model type
varimp_summary <- varimp_all %>%
  group_by(Variable, Model) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 1), .groups = "drop") %>%
  arrange(Model, desc(MedianImportance)) %>%
  group_by(Model) %>%
  mutate(Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 1),
         Cumulative_pct = round(cumsum(100 * MedianImportance / sum(MedianImportance)), 1)) %>%
  ungroup()

print(varimp_summary)
clipr::write_clip(varimp_summary)

# ===== PREDICTION ANALYSIS (TEST DATA ONLY) =====
# Get test predictions only
predictions_test <- rbindlist(predictions_all[!sapply(predictions_all, is.null)])[set == "test"]

if (nrow(predictions_test) > 0) {
  print("=== PREDICTION ANALYSIS ===")
  
  # ===== COVERAGE PLOT BY ABUNDANCE BINS =====
  # Use first test dataset for binning analysis
  test_data_sample <- test_data_list[[1]]
  present_idx <- which(test_data_sample$Presence == "1")
  test_present_sample <- test_data_sample[present_idx, ]
  
  if (nrow(test_present_sample) > 10) {
    # Create bins for coverage analysis
    breaks <- seq(0, max(test_present_sample$Sanionia, na.rm = TRUE), length.out = 6)
    if (length(breaks) < 2) breaks <- c(0, max(test_present_sample$Sanionia, na.rm = TRUE))
    
    obs_binned <- cut(test_present_sample$Sanionia, breaks = breaks, include.lowest = TRUE)
    
    # Calculate coverage using predictions from first iteration
    pred_sample <- predictions_test[iteration == 1]
    if (nrow(pred_sample) > 0) {
      # Calculate prediction intervals (assuming spread represents interval width)
      pred_lower <- pred_sample$pred - pred_sample$spread/2
      pred_upper <- pred_sample$pred + pred_sample$spread/2
      
      coverage_by_bin <- tapply(pred_sample$obs >= pred_lower & pred_sample$obs <= pred_upper, 
                                obs_binned[1:nrow(pred_sample)], mean, na.rm = TRUE) * 100
      
      # Create bin labels
      bin_labels <- paste0("[", round(breaks[-length(breaks)], 1), ",", round(breaks[-1], 1), ")")
      bin_labels[length(bin_labels)] <- paste0("[", round(breaks[length(breaks)-1], 1), ",", 
                                               round(breaks[length(breaks)], 1), "]")
      
      # Coverage plot with larger labels
      par(cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
      bp <- barplot(coverage_by_bin, col = "#EBCC2A", ylim = c(0, 100),
                    main = "95% Prediction Interval Coverage by Abundance Bins",
                    xlab = "Observed Sanionia (Binned)", ylab = "Coverage (%)",
                    names.arg = "", xaxt = "n")
      text(bp, par("usr")[3] - 5, labels = bin_labels, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 1.1)
      abline(h = 95, col = "red", lty = 2, lwd = 2)
      # Reset par
      par(cex.main = 1, cex.lab = 1, cex.axis = 1)
    }
  }
  
  # ===== UNCERTAINTY VS ABUNDANCE PLOT =====
  spread_plot <- ggplot(predictions_test, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.3, color = "#EBCC2A") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(x = "Observed Sanionia", y = "Prediction Interval Width",
         title = "Prediction Uncertainty vs Observed Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(spread_plot)
  
  # ===== OBSERVED VS PREDICTED PLOT =====
  obs_pred <- predictions_test %>%
    group_by(obs) %>%
    summarize(pred = median(pred, na.rm = TRUE), .groups = "drop")
  
  spearman_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "spearman", use = "complete.obs")^2
  pearson_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "pearson", use = "complete.obs")^2
  
  r2_label <- paste0("Spearman R² = ", round(spearman_r2, 3), "\nPearson R² = ", round(pearson_r2, 3))
  
  obs_vs_pred_plot <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.6, color = "#EBCC2A") +
    geom_smooth(method = "lm", se = FALSE, color = "#EBCC2A") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = -Inf, label = r2_label,
             hjust = 1.1, vjust = -0.5, size = 5, color = "blue") +
    labs(x = "Observed Sanionia", y = "Predicted Sanionia",
         title = "Observed vs Predicted Sanionia (Test Data)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(obs_vs_pred_plot)
  
  # ===== RESIDUALS PLOT =====
  obs_pred$residuals <- obs_pred$obs - obs_pred$pred
  
  residuals_plot <- ggplot(obs_pred, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.6, color = "#EBCC2A") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "#EBCC6A") +
    labs(x = "Predicted Sanionia", y = "Residuals (Observed - Predicted)",
         title = "Residuals vs Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(residuals_plot)
}



# Load the saved results
all_results_clf <- readRDS("hurdle_rf_classifier_results_Deschampsia.rds")
all_results_qrf <- readRDS("hurdle_qrf_abundance_results_Deschampsia.rds")
varimp_list_clf <- readRDS("hurdle_varimp_rf_presence_Deschampsia.rds")
varimp_list_qrf <- readRDS("hurdle_varimp_qrf_abundance_Deschampsia.rds")
predictions_all <- readRDS("hurdle_predictions_all_Deschampsia.rds")
training_data_list <- readRDS("hurdle_training_data_list_Deschampsia.rds")
test_data_list <- readRDS("hurdle_test_data_list_Deschampsia.rds")

# ===== CLASSIFIER METRICS SUMMARY =====
print("=== RF CLASSIFIER PERFORMANCE SUMMARY ===")
metrics_clf <- rbindlist(all_results_clf)
metrics_clf_summary <- metrics_clf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_clf_summary$Model <- "RF_Presence"
print(metrics_clf_summary)
clipr::write_clip(metrics_clf_summary)

# ===== QRF ABUNDANCE METRICS SUMMARY =====
print("=== QRF ABUNDANCE PERFORMANCE SUMMARY ===")
metrics_qrf <- rbindlist(all_results_qrf)
metrics_qrf_summary <- metrics_qrf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_qrf_summary$Model <- "QRF_Abundance"
metrics_qrf_summary$Transformation <- "Raw"
print(metrics_qrf_summary)
clipr::write_clip(metrics_qrf_summary)

# ===== VARIABLE IMPORTANCE SUMMARY =====
print("=== VARIABLE IMPORTANCE SUMMARY ===")

# Combine classifier and QRF importance
varimp_clf <- rbindlist(varimp_list_clf)
varimp_qrf <- rbindlist(varimp_list_qrf[!sapply(varimp_list_qrf, is.null)])
varimp_all <- rbind(varimp_clf, varimp_qrf)

# Summarize by model type
varimp_summary <- varimp_all %>%
  group_by(Variable, Model) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 1), .groups = "drop") %>%
  arrange(Model, desc(MedianImportance)) %>%
  group_by(Model) %>%
  mutate(Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 1),
         Cumulative_pct = round(cumsum(100 * MedianImportance / sum(MedianImportance)), 1)) %>%
  ungroup()
print(varimp_summary, n = Inf)
print(varimp_summary)
clipr::write_clip(varimp_summary)

# ===== PREDICTION ANALYSIS (TEST DATA ONLY) =====
# Get test predictions only
predictions_test <- rbindlist(predictions_all[!sapply(predictions_all, is.null)])[set == "test"]

if (nrow(predictions_test) > 0) {
  print("=== PREDICTION ANALYSIS ===")
  
  # ===== COVERAGE PLOT BY ABUNDANCE BINS =====
  # Use first test dataset for binning analysis
  test_data_sample <- test_data_list[[1]]
  present_idx <- which(test_data_sample$Presence == "1")
  test_present_sample <- test_data_sample[present_idx, ]
  
  if (nrow(test_present_sample) > 10) {
    # Create bins for coverage analysis
    breaks <- seq(0, max(test_present_sample$Deschampsia, na.rm = TRUE), length.out = 6)
    if (length(breaks) < 2) breaks <- c(0, max(test_present_sample$Deschampsia, na.rm = TRUE))
    
    obs_binned <- cut(test_present_sample$Deschampsia, breaks = breaks, include.lowest = TRUE)
    
    # Calculate coverage using predictions from first iteration
    pred_sample <- predictions_test[iteration == 1]
    if (nrow(pred_sample) > 0) {
      # Calculate prediction intervals (assuming spread represents interval width)
      pred_lower <- pred_sample$pred - pred_sample$spread/2
      pred_upper <- pred_sample$pred + pred_sample$spread/2
      
      coverage_by_bin <- tapply(pred_sample$obs >= pred_lower & pred_sample$obs <= pred_upper, 
                                obs_binned[1:nrow(pred_sample)], mean, na.rm = TRUE) * 100
      
      # Create bin labels
      bin_labels <- paste0("[", round(breaks[-length(breaks)], 1), ",", round(breaks[-1], 1), ")")
      bin_labels[length(bin_labels)] <- paste0("[", round(breaks[length(breaks)-1], 1), ",", 
                                               round(breaks[length(breaks)], 1), "]")
      
      # Coverage plot with larger labels
      par(cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
      bp <- barplot(coverage_by_bin, col = "#F21A00", ylim = c(0, 100),
                    main = "95% Prediction Interval Coverage by Abundance Bins",
                    xlab = "Observed Deschampsia (Binned)", ylab = "Coverage (%)",
                    names.arg = "", xaxt = "n")
      text(bp, par("usr")[3] - 5, labels = bin_labels, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 1.1)
      abline(h = 95, col = "red", lty = 2, lwd = 2)
      # Reset par
      par(cex.main = 1, cex.lab = 1, cex.axis = 1)
    }
  }
  
  # ===== UNCERTAINTY VS ABUNDANCE PLOT =====
  spread_plot <- ggplot(predictions_test, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.3, color = "#F21A00") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(x = "Observed Deschampsia", y = "Prediction Interval Width",
         title = "Prediction Uncertainty vs Observed Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(spread_plot)
  
  # ===== OBSERVED VS PREDICTED PLOT =====
  obs_pred <- predictions_test %>%
    group_by(obs) %>%
    summarize(pred = median(pred, na.rm = TRUE), .groups = "drop")
  
  spearman_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "spearman", use = "complete.obs")^2
  pearson_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "pearson", use = "complete.obs")^2
  
  r2_label <- paste0("Spearman R² = ", round(spearman_r2, 3), "\nPearson R² = ", round(pearson_r2, 3))
  
  obs_vs_pred_plot <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.6, color = "#F21A00") +
    geom_smooth(method = "lm", se = FALSE, color = "#F21A00") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = -Inf, label = r2_label,
             hjust = 1.1, vjust = -0.5, size = 5, color = "blue") +
    labs(x = "Observed Deschampsia", y = "Predicted Deschampsia",
         title = "Observed vs Predicted Deschampsia (Test Data)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(obs_vs_pred_plot)
  
  # ===== RESIDUALS PLOT =====
  obs_pred$residuals <- obs_pred$obs - obs_pred$pred
  
  residuals_plot <- ggplot(obs_pred, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.6, color = "#F21A00") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "#F26A00") +
    labs(x = "Predicted Deschampsia", y = "Residuals (Observed - Predicted)",
         title = "Residuals vs Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(residuals_plot)
}


# Load the saved results
all_results_clf <- readRDS("hurdle_rf_classifier_results_Colobanthus_top_predictors.rds")
all_results_qrf <- readRDS("hurdle_qrf_abundance_results_Colobanthus_top_predictors.rds")
varimp_list_clf <- readRDS("hurdle_varimp_rf_presence_Colobanthus_top_predictors.rds")
varimp_list_qrf <- readRDS("hurdle_varimp_qrf_abundance_Colobanthus_top_predictors.rds")
predictions_all <- readRDS("hurdle_predictions_all_Colobanthus_top_predictors.rds")
training_data_list <- readRDS("hurdle_training_data_list_Colobanthus_top_predictors.rds")
test_data_list <- readRDS("hurdle_test_data_list_Colobanthus_top_predictors.rds")

# ===== CLASSIFIER METRICS SUMMARY =====
print("=== RF CLASSIFIER PERFORMANCE SUMMARY ===")
metrics_clf <- rbindlist(all_results_clf)
metrics_clf_summary <- metrics_clf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_clf_summary$Model <- "RF_Presence"
print(metrics_clf_summary)
clipr::write_clip(metrics_clf_summary)

# ===== QRF ABUNDANCE METRICS SUMMARY =====
print("=== QRF ABUNDANCE PERFORMANCE SUMMARY ===")
metrics_qrf <- rbindlist(all_results_qrf)
metrics_qrf_summary <- metrics_qrf[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), .SDcols = -"Iteration"]
metrics_qrf_summary$Model <- "QRF_Abundance"
metrics_qrf_summary$Transformation <- "Raw"
print(metrics_qrf_summary)
clipr::write_clip(metrics_qrf_summary)

# ===== VARIABLE IMPORTANCE SUMMARY =====
print("=== VARIABLE IMPORTANCE SUMMARY ===")

# Combine classifier and QRF importance
varimp_clf <- rbindlist(varimp_list_clf)
varimp_qrf <- rbindlist(varimp_list_qrf[!sapply(varimp_list_qrf, is.null)])
varimp_all <- rbind(varimp_clf, varimp_qrf)

# Summarize by model type
varimp_summary <- varimp_all %>%
  group_by(Variable, Model) %>%
  summarize(MedianImportance = round(median(Importance, na.rm = TRUE), 1), .groups = "drop") %>%
  arrange(Model, desc(MedianImportance)) %>%
  group_by(Model) %>%
  mutate(Importance_pct = round(100 * MedianImportance / sum(MedianImportance), 1),
         Cumulative_pct = round(cumsum(100 * MedianImportance / sum(MedianImportance)), 1)) %>%
  ungroup()
print(varimp_summary, n = Inf)
print(varimp_summary)
clipr::write_clip(varimp_summary)

# ===== PREDICTION ANALYSIS (TEST DATA ONLY) =====
# Get test predictions only
predictions_test <- rbindlist(predictions_all[!sapply(predictions_all, is.null)])[set == "test"]

if (nrow(predictions_test) > 0) {
  print("=== PREDICTION ANALYSIS ===")
  
  # ===== COVERAGE PLOT BY ABUNDANCE BINS =====
  # Use first test dataset for binning analysis
  test_data_sample <- test_data_list[[1]]
  present_idx <- which(test_data_sample$Presence == "1")
  test_present_sample <- test_data_sample[present_idx, ]
  
  if (nrow(test_present_sample) > 10) {
    # Create bins for coverage analysis
    breaks <- seq(0, max(test_present_sample$Colobanthus, na.rm = TRUE), length.out = 6)
    if (length(breaks) < 2) breaks <- c(0, max(test_present_sample$Colobanthus, na.rm = TRUE))
    
    obs_binned <- cut(test_present_sample$Colobanthus, breaks = breaks, include.lowest = TRUE)
    
    # Calculate coverage using predictions from first iteration
    pred_sample <- predictions_test[iteration == 1]
    if (nrow(pred_sample) > 0) {
      # Calculate prediction intervals (assuming spread represents interval width)
      pred_lower <- pred_sample$pred - pred_sample$spread/2
      pred_upper <- pred_sample$pred + pred_sample$spread/2
      
      coverage_by_bin <- tapply(pred_sample$obs >= pred_lower & pred_sample$obs <= pred_upper, 
                                obs_binned[1:nrow(pred_sample)], mean, na.rm = TRUE) * 100
      
      # Create bin labels
      bin_labels <- paste0("[", round(breaks[-length(breaks)], 1), ",", round(breaks[-1], 1), ")")
      bin_labels[length(bin_labels)] <- paste0("[", round(breaks[length(breaks)-1], 1), ",", 
                                               round(breaks[length(breaks)], 1), "]")
      
      # Coverage plot with larger labels
      par(cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)
      bp <- barplot(coverage_by_bin, col = "#78B7C5", ylim = c(0, 100),
                    main = "95% Prediction Interval Coverage by Abundance Bins",
                    xlab = "Observed Colobanthus (Binned)", ylab = "Coverage (%)",
                    names.arg = "", xaxt = "n")
      text(bp, par("usr")[3] - 5, labels = bin_labels, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 1.1)
      abline(h = 95, col = "red", lty = 2, lwd = 2)
      # Reset par
      par(cex.main = 1, cex.lab = 1, cex.axis = 1)
    }
  }
  
  # ===== UNCERTAINTY VS ABUNDANCE PLOT =====
  spread_plot <- ggplot(predictions_test, aes(x = obs, y = spread)) +
    geom_point(alpha = 0.3, color = "#78B7C5") +
    geom_smooth(method = "loess", se = FALSE, color = "#78B7C5") +
    labs(x = "Observed Colobanthus", y = "Prediction Interval Width",
         title = "Prediction Uncertainty vs Observed Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(spread_plot)
  
  # ===== OBSERVED VS PREDICTED PLOT =====
  obs_pred <- predictions_test %>%
    group_by(obs) %>%
    summarize(pred = median(pred, na.rm = TRUE), .groups = "drop")
  
  spearman_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "spearman", use = "complete.obs")^2
  pearson_r2 <- cor(obs_pred$obs, obs_pred$pred, method = "pearson", use = "complete.obs")^2
  
  r2_label <- paste0("Spearman R² = ", round(spearman_r2, 3), "\nPearson R² = ", round(pearson_r2, 3))
  
  obs_vs_pred_plot <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.6, color = "#78B7C5") +
    geom_smooth(method = "lm", se = FALSE, color = "#78B7C5") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = -Inf, label = r2_label,
             hjust = 1.1, vjust = -0.5, size = 5, color = "blue") +
    labs(x = "Observed Colobanthus", y = "Predicted Colobanthus",
         title = "Observed vs Predicted Colobanthus (Test Data)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(obs_vs_pred_plot)
  
  # ===== RESIDUALS PLOT =====
  obs_pred$residuals <- obs_pred$obs - obs_pred$pred
  
  residuals_plot <- ggplot(obs_pred, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.6, color = "#78B7C5") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "#45B0C5") +
    labs(x = "Predicted Colobanthus", y = "Residuals (Observed - Predicted)",
         title = "Residuals vs Predicted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(residuals_plot)
}

