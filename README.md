# Antarctic Species Distribution Models

## Description
R scripts for modeling species distribution using VIF-filtered environmental variables, Quantile Random Forest hurdle models and niche analysis.

## Files
- `vif_analysis.R` - Variable selection using VIF thresholds and Spearman rank correlation to ensure variables highly correlated to abundance were prioritised.
- `Resdiual_plpts.R` - Used to asses model performance
- `Hurdle_models_v GLM.R` - Hurdle models for presence and abundance predictions followed by OMI niche analysis and partial response plots
- `Hurdle_Current_Predictions.R` - Prediction hurdle models for presence and abundance predictions 


## Requirements
- R 4.0+
- Required packages: dplyr, ade4, randomForest, terra, quantregForest, isotone, pROC, pdp, usdm

## Usage
1. Run VIF analysis first: `source("vif_analysis.R")`
2. Run the Quantile Random Forest (QRF) hurdle model using abundance data and VIF filtered environmental variables: `source("Hurdle_current_predictions.R")`
3. Generate niche plots using the important environmental variables used in the abundnace component of the hurdle model: `source("niche_analysis.R")`


