# Load required packages

library(usdm)
library(dplyr)

rm(list=ls()) 

env_data <- read.csv("Environmental_Variables.csv")
species_data <-  read.csv("Abundance_data.csv",header=T,sep=",")

names(env_data)
names(species_data)

a <- env_data %>%
  filter(!SITE_FINAL %in% species_data$SITE_QUAD)

a$SITE_FINAL


# --- STEP 0: Merge environmental and species data by SITE_FINAL ---
# Replace 'env_data' and 'species_data' with your actual data frame names
merged_data <- left_join(env_data, species_data, by = c("SITE_FINAL" = "SITE_QUAD"))
### check for cor with Usnea

names(merged_data)

# Load required packages
library(dplyr)

# --- USER INPUTS ---
response_var <- "Deschampsia"          # Your species column name
predictor_cols <- 3:104           # Adjust this to your env var columns
cor_method <- "spearman"         # "pearson" or "spearman"

str(merged_data)

# --- STEP 1: Extract relevant data ---
# Get numeric environmental variables
env_vars <- merged_data[, predictor_cols] %>% 
  select(where(is.numeric))

df <- merged_data%>%
  select(all_of(response_var), all_of(names(env_vars)))

# --- STEP 2: Compute correlations ---
cor_table <- sapply(names(env_vars), function(var) {
  cor(df[[response_var]], df[[var]], method = cor_method, use = "complete.obs")
})

# --- STEP 3: Format and display results ---
cor_df <- data.frame(
  Variable = names(cor_table),
  Spearman_r = round(as.numeric(cor_table), 3)
) %>%
  arrange(desc(abs(Spearman_r)))

# Print nicely
cat(paste0("Spearman R between ", response_var, " and environmental variables:\n\n"))
print(cor_df, row.names = FALSE)


### now VIF



# --- USER PARAMETERS ---
response_var <- "Species"           # Name of species abundance column
predictor_cols <- 3:104            # Adjust based on merged_data column positions
vif_threshold <- 5
cor_method <- "spearman"

# --- STEP 1: Select numeric environmental predictors ---
numeric_vars <- merged_data[, predictor_cols] %>%
  dplyr::select(where(is.numeric))

# --- STEP 2: Compute Spearman correlation with vegetation response ---
cor_vec <- sapply(numeric_vars, function(x)
  cor(x, merged_data[[response_var]], method = cor_method, use = "complete.obs"))

# Rank predictors by correlation with species
cor_rank <- sort(abs(cor_vec), decreasing = TRUE)
ordered_vars <- names(cor_rank)

# --- STEP 3: Iterative VIF-based selection ---
selected_vars <- c()

for (var in ordered_vars) {
  candidate_vars <- c(selected_vars, var)
  
  if (length(candidate_vars) == 1) {
    selected_vars <- candidate_vars
  } else {
    vif_df <- usdm::vif(numeric_vars[, candidate_vars, drop = FALSE])
    max_vif <- max(vif_df$VIF, na.rm = TRUE)
    
    if (max_vif <= vif_threshold) {
      selected_vars <- candidate_vars
    }
  }
}

str(merged_data)

final_vars <- c(
  "SITE_FINAL",                   # always keep
  names(merged_data)[!sapply(merged_data, is.numeric)],  # all non-numeric
  selected_vars,                 # VIF-selected numeric variables
  response_var                   # species abundance
) %>% unique()

Deschampsia_filtered <- merged_data %>%
  select(all_of(final_vars), -IMAGE)

## save
write.csv(Species_filtered, "Species_filtered_VIF.csv", row.names = FALSE)

# --- STEP 5: Report ---
cat("\n Final selected predictors (VIF ≤", vif_threshold, "):\n")
print(selected_vars)

cat("\n Spearman correlations with", response_var, ":\n")
print(cor_vec[selected_vars])

cat("\n Final VIF values:\n")
print(usdm::vif(numeric_vars[, selected_vars, drop = FALSE]))



#########################################################################################################################################################################



# ─────────────────────────────────────────────────────────────────────────────
# SANIONIA: IMPROVED SPEARMAN + PRIORITY-BASED VIF SELECTION
# ─────────────────────────────────────────────────────────────────────────────

# 0. Clean workspace
rm(list=ls())
library(dplyr)

# 1. File paths – EDIT THESE
env_path     <- "Environmental_Variables.csv"
species_path <- "Abundance_Data.csv"

# 2. Read & label
raw           <- read.csv(env_path, header=FALSE, stringsAsFactors=FALSE)
var_names     <- as.character(raw[2,])
var_type_row  <- setNames(as.character(raw[1,]), var_names)
colnames(raw) <- var_names
dados         <- raw[-c(1,2), ]

# 3. Convert numeric‐looking columns
numflag <- sapply(dados, function(x)
  !all(is.na(suppressWarnings(as.numeric(x)))))
for(col in names(dados)[numflag]) {
  dados[[col]] <- as.numeric(dados[[col]])
}

# 4. Merge species data
sp        <- read.csv(species_path, stringsAsFactors=FALSE)
site_col  <- if("SITE_FINAL" %in% names(dados)) "SITE_FINAL" else
  grep("SITE", names(dados), ignore.case=TRUE, value=TRUE)[1]
combined  <- merge(dados, sp,
                   by.x=site_col, by.y="SITE_QUAD", all=FALSE)

# 5. Split predictors
env_start <- 15
env_end   <- ncol(combined) - 5
env_block <- combined[, env_start:env_end]
is_num    <- sapply(env_block, is.numeric)
env_num   <- env_block[, is_num, drop=FALSE]
env_cat   <- env_block[, !is_num, drop=FALSE]

# drop zero/infinite‐variance
keep_num <- sapply(env_num, function(x){
  v <- var(x, na.rm=TRUE); is.finite(v) && v>0
})
env_num <- env_num[, keep_num, drop=FALSE]

# exclude IMAGE & response from predictors
env_num <- env_num[, setdiff(names(env_num), c("IMAGE", "Sanionia")), drop=FALSE]

# 6. Extract Sanionia abundance
sanionia <- combined$Sanionia

# 7. Step 1 – Spearman among env vars (no species) - Threshold for microscale-microscale correlations
env_cor_thr <- 0.85  # Higher threshold for microscale-microscale correlations only
env_cor_mat <- cor(env_num, method="spearman", use="pairwise.complete.obs")
pairs_idx   <- which(upper.tri(env_cor_mat) & abs(env_cor_mat) > env_cor_thr,
                     arr.ind=TRUE)

if(nrow(pairs_idx)>0) {
  env_pairs <- data.frame(
    Var1 = rownames(env_cor_mat)[pairs_idx[,1]],
    Var2 = colnames(env_cor_mat)[pairs_idx[,2]],
    rho  = env_cor_mat[pairs_idx],
    stringsAsFactors=FALSE
  )
} else {
  env_pairs <- data.frame(
    Var1=character(), Var2=character(), rho=numeric(),
    stringsAsFactors=FALSE
  )
}

cat("\nStep 1: Highly inter-correlated environmental variables (|ρ| > ",
    env_cor_thr, "):\n", sep="")
print(env_pairs)

# 8. Step 2 – Spearman vs Sanionia
spearman_all <- function(vars, df, y){
  out <- data.frame(Variable=character(), rho=numeric(), p=numeric(),
                    stringsAsFactors=FALSE)
  for(v in vars){
    x <- df[[v]]
    if(length(unique(na.omit(x)))>1 &&
       length(unique(na.omit(y)))>1){
      tst <- try(cor.test(x, y, method="spearman", exact=FALSE),
                 silent=TRUE)
      if(!inherits(tst, "try-error")){
        out <- rbind(out, data.frame(
          Variable=v,
          rho      = as.numeric(tst$estimate),
          p        = tst$p.value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  out$p_adj <- p.adjust(out$p, method="BH")
  out <- out[order(-abs(out$rho)), ]
  rownames(out) <- NULL
  return(out)
}

all_vars   <- names(env_num)
cor_df_all <- spearman_all(all_vars, env_num, sanionia)

cat("\nStep 2: Spearman correlations vs Sanionia:\n")
print(cor_df_all)

# 9. Step 3 – IMPROVED: Create comprehensive priority ranking
cat("\nStep 3: Creating comprehensive priority ranking...\n")

# Create a dataframe with all variables and their species correlations
priority_df <- data.frame(
  Variable = cor_df_all$Variable,
  abs_rho = abs(cor_df_all$rho),
  rho = cor_df_all$rho,
  stringsAsFactors = FALSE
)

# Identify variables involved in high correlations
correlated_vars <- unique(c(env_pairs$Var1, env_pairs$Var2))

# For each variable, determine if it should be deprioritized
priority_df$status <- "available"
priority_df$reason <- ""

for(i in seq_len(nrow(env_pairs))){
  v1 <- env_pairs$Var1[i]
  v2 <- env_pairs$Var2[i]
  
  # Get absolute correlations with species
  r1 <- priority_df$abs_rho[priority_df$Variable == v1]
  r2 <- priority_df$abs_rho[priority_df$Variable == v2]
  
  if(length(r1) > 0 && length(r2) > 0){
    # The variable with LOWER species correlation gets deprioritized
    if(r1 > r2){
      priority_df$status[priority_df$Variable == v2] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Lower correlation with species than ", v1)
    } else if(r2 > r1){
      priority_df$status[priority_df$Variable == v1] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v1] <- paste0("Lower correlation with species than ", v2)
    } else {
      # If equal, deprioritize the second one (arbitrary but consistent)
      priority_df$status[priority_df$Variable == v2] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Equal correlation with species as ", v1, " (tie-breaker)")
    }
  }
}

# Sort by priority: available variables first, then by absolute correlation
priority_df <- priority_df[order(priority_df$status, -priority_df$abs_rho), ]

cat("\nPriority ranking (top 15):\n")
print(head(priority_df, 15))

# 10. Compute initial VIF for all env predictors (for reference only)
manual_vif <- function(df){
  vars <- names(df)
  sapply(vars, function(v){
    rest <- setdiff(vars, v)
    if(length(rest)==0) return(1)
    
    # Handle cases where perfect multicollinearity exists
    tryCatch({
      fit <- lm(df[[v]] ~ ., data=df[,rest,drop=FALSE])
      R2  <- summary(fit)$r.squared
      
      # Cap R² at 0.9999 to avoid infinite VIF
      if(R2 >= 0.9999) R2 <- 0.9999
      
      1/(1-R2)
    }, error = function(e) {
      return(Inf)  # Return Inf if model fails
    })
  })
}

vif_all <- manual_vif(env_num)
cat("\nInitial VIF values (for reference):\n")
print(sort(vif_all, decreasing = TRUE))

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")
cat("Note: Being more flexible to include more meaningful variables\n")

# Don't pre-filter by VIF - let the selection process handle it
all_available_vars <- priority_df$Variable[priority_df$status == "available"]
all_deprioritized_vars <- priority_df$Variable[priority_df$status == "deprioritised"]

# Start with empty selection
selected <- character()

# Phase 1: Add available (prioritized) variables in order of species correlation
cat("\n--- Phase 1: Adding prioritized variables ---\n")
for(v in all_available_vars){
  if(length(selected) == 0){
    # First variable - always add it
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    # Check if adding this variable creates acceptable VIF
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Also check if any existing variables would exceed VIF threshold
    max_existing_vif <- max(vifs[selected])
    
    if(new_var_vif <= vif_thr && max_existing_vif <= vif_thr){
      selected <- trial
      cat("Added (available):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    }
  }
}

# Phase 2: Try to add deprioritized variables that don't violate VIF
cat("\n--- Phase 2: Adding deprioritized variables where possible ---\n")
for(v in all_deprioritized_vars){
  if(length(selected) == 0){
    # If somehow no variables selected yet, add this one
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Check if any variable (new or existing) would exceed VIF threshold
    max_all_vif <- max(vifs)
    
    if(new_var_vif <= vif_thr && max_all_vif <= vif_thr){
      selected <- trial
      cat("Added (deprioritized):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    }
  }
}

# Phase 3: Final check - if we still have very few variables, try a different approach
cat("\n--- Phase 3: Final optimization ---\n")
if(length(selected) < 5){
  cat("Only", length(selected), "variables selected. Trying alternative approach...\n")
  
  # Try starting with the most important variables and build around them
  top_vars <- head(priority_df$Variable, 15)  # Top 15 by species correlation
  
  # Test different combinations starting from different points
  best_set <- selected
  best_count <- length(selected)
  
  for(i in 1:min(5, length(top_vars))){  # Try top 5 as starting points
    start_var <- top_vars[i]
    temp_selected <- start_var
    
    # Add other variables one by one, prioritizing by species correlation
    remaining_vars <- setdiff(top_vars, start_var)
    for(next_var in remaining_vars){
      trial <- unique(c(temp_selected, next_var))
      if(length(trial) > 1){
        vifs <- manual_vif(env_num[, trial, drop=FALSE])
        new_var_vif <- vifs[next_var]
        
        if(new_var_vif <= vif_thr){
          temp_selected <- trial
        }
      }
    }
    
    # If this combination is better, use it
    if(length(temp_selected) > best_count){
      best_set <- temp_selected
      best_count <- length(temp_selected)
    }
  }
  
  selected <- best_set
  cat("Alternative approach yielded", length(selected), "variables\n")
}

# 12. Final results
final_vifs <- if(length(selected) > 0)
  manual_vif(env_num[, selected, drop=FALSE]) else numeric()

results_tbl <- data.frame(
  Variable = selected,
  VIF = round(as.numeric(final_vifs), 3),
  Spearman_rho = round(cor_df_all$rho[match(selected, cor_df_all$Variable)], 3),
  Abs_Spearman_rho = round(abs(cor_df_all$rho[match(selected, cor_df_all$Variable)]), 3),
  Priority_Status = priority_df$status[match(selected, priority_df$Variable)],
  Var_Type = var_type_row[selected],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Abs_Spearman_rho))

cat("\nStep 4: Final selected variables (VIF < ", vif_thr, "):\n", sep="")
print(results_tbl)

# Show summary
cat("\nSummary:\n")
cat("Total variables selected:", nrow(results_tbl), "\n")
cat("Available variables selected:", sum(results_tbl$Priority_Status == "available"), "\n")
cat("Deprioritized variables selected:", sum(results_tbl$Priority_Status == "deprioritised"), "\n")

# 13. Save final dataset
id_cols <- setdiff(names(combined)[1:(env_start-1)], "IMAGE")
final_dataset <- combined %>%
  select(all_of(id_cols)) %>%
  bind_cols(env_num[, selected, drop=FALSE], env_cat) %>%
  mutate(Sanionia = sanionia)

write.csv(results_tbl, "Sanionia_VIF_results_improved.csv", row.names=FALSE)
write.csv(final_dataset, "VIF_Sanionia_FINAL_improved.csv", row.names=FALSE)
write.csv(priority_df, "Variable_Priority_Ranking.csv", row.names=FALSE)

message("Done: retained ", nrow(results_tbl), " variables with VIF < ", vif_thr, ".")



# ─────────────────────────────────────────────────────────────────────────────
# USNEA: STEPWISE SPEARMAN + PRIORITY-BASED VIF SELECTION
# ─────────────────────────────────────────────────────────────────────────────

# 0. Clean workspace
rm(list=ls())
library(dplyr)

# 1. File paths – EDIT THESE
env_path     <- "Environmental_Variables.csv"
species_path <- "Abundance_Data.csv"

# 2. Read & label
raw           <- read.csv(env_path, header=FALSE, stringsAsFactors=FALSE)
var_names     <- as.character(raw[2,])
var_type_row  <- setNames(as.character(raw[1,]), var_names)
colnames(raw) <- var_names
dados         <- raw[-c(1,2), ]

# 3. Convert numeric‐looking columns
numflag <- sapply(dados, function(x)
  !all(is.na(suppressWarnings(as.numeric(x)))))
for(col in names(dados)[numflag]) {
  dados[[col]] <- as.numeric(dados[[col]])
}

# 4. Merge species data
sp        <- read.csv(species_path, stringsAsFactors=FALSE)
site_col  <- if("SITE_FINAL" %in% names(dados)) "SITE_FINAL" else
  grep("SITE", names(dados), ignore.case=TRUE, value=TRUE)[1]
combined  <- merge(dados, sp,
                   by.x=site_col, by.y="SITE_QUAD", all=FALSE)

# 5. Split predictors
env_start <- 15
env_end   <- ncol(combined) - 5
env_block <- combined[, env_start:env_end]
is_num    <- sapply(env_block, is.numeric)
env_num   <- env_block[, is_num, drop=FALSE]
env_cat   <- env_block[, !is_num, drop=FALSE]

# drop zero/infinite‐variance
keep_num <- sapply(env_num, function(x){
  v <- var(x, na.rm=TRUE); is.finite(v) && v>0
})
env_num <- env_num[, keep_num, drop=FALSE]

# exclude IMAGE & response from predictors
env_num <- env_num[, setdiff(names(env_num), c("IMAGE", "Usnea")), drop=FALSE]

# 6. Extract Usnea abundance
usnea <- combined$Usnea

# 7. Step 1 – Spearman among env vars (no species) - Threshold for microscale-microscale correlations
env_cor_thr <- 0.85  # Higher threshold for microscale-microscale correlations only
env_cor_mat <- cor(env_num, method="spearman", use="pairwise.complete.obs")
pairs_idx   <- which(upper.tri(env_cor_mat) & abs(env_cor_mat) > env_cor_thr,
                     arr.ind=TRUE)

if(nrow(pairs_idx)>0) {
  env_pairs <- data.frame(
    Var1 = rownames(env_cor_mat)[pairs_idx[,1]],
    Var2 = colnames(env_cor_mat)[pairs_idx[,2]],
    rho  = env_cor_mat[pairs_idx],
    stringsAsFactors=FALSE
  )
} else {
  env_pairs <- data.frame(
    Var1=character(), Var2=character(), rho=numeric(),
    stringsAsFactors=FALSE
  )
}

cat("\nStep 1: Highly inter-correlated environmental variables (|ρ| > ",
    env_cor_thr, "):\n", sep="")
print(env_pairs)

# 8. Step 2 – Spearman vs Usnea
spearman_all <- function(vars, df, y){
  out <- data.frame(Variable=character(), rho=numeric(), p=numeric(),
                    stringsAsFactors=FALSE)
  for(v in vars){
    x <- df[[v]]
    if(length(unique(na.omit(x)))>1 &&
       length(unique(na.omit(y)))>1){
      tst <- try(cor.test(x, y, method="spearman", exact=FALSE),
                 silent=TRUE)
      if(!inherits(tst, "try-error")){
        out <- rbind(out, data.frame(
          Variable=v,
          rho      = as.numeric(tst$estimate),
          p        = tst$p.value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  out$p_adj <- p.adjust(out$p, method="BH")
  out <- out[order(-abs(out$rho)), ]
  rownames(out) <- NULL
  return(out)
}

all_vars   <- names(env_num)
cor_df_all <- spearman_all(all_vars, env_num, usnea)

cat("\nStep 2: Spearman correlations vs Usnea:\n")
print(cor_df_all)

# 9. Step 3 – IMPROVED: Create comprehensive priority ranking
cat("\nStep 3: Creating comprehensive priority ranking...\n")

# Create a dataframe with all variables and their Usnea correlations
priority_df <- data.frame(
  Variable = cor_df_all$Variable,
  abs_rho = abs(cor_df_all$rho),
  rho = cor_df_all$rho,
  stringsAsFactors = FALSE
)

# Identify variables involved in high correlations
correlated_vars <- unique(c(env_pairs$Var1, env_pairs$Var2))

# For each variable, determine if it should be deprioritised
priority_df$status <- "available"
priority_df$reason <- ""

for(i in seq_len(nrow(env_pairs))){
  v1 <- env_pairs$Var1[i]
  v2 <- env_pairs$Var2[i]
  
  # Get absolute correlations with Usnea
  r1 <- priority_df$abs_rho[priority_df$Variable == v1]
  r2 <- priority_df$abs_rho[priority_df$Variable == v2]
  
  if(length(r1) > 0 && length(r2) > 0){
    # The variable with LOWER Usnea correlation gets deprioritised
    if(r1 > r2){
      priority_df$status[priority_df$Variable == v2] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Lower correlation with Usnea than ", v1)
    } else if(r2 > r1){
      priority_df$status[priority_df$Variable == v1] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v1] <- paste0("Lower correlation with Usnea than ", v2)
    } else {
      # If equal, deprioritize the second one (arbitrary but consistent)
      priority_df$status[priority_df$Variable == v2] <- "deprioritized"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Equal correlation with Usnea as ", v1, " (tie-breaker)")
    }
  }
}

# Sort by priority: available variables first, then by absolute Usnea correlation
priority_df <- priority_df[order(priority_df$status, -priority_df$abs_rho), ]

cat("\nPriority ranking (top 15):\n")
print(head(priority_df, 15))

# 10. Compute initial VIF for all env predictors (for reference only)
manual_vif <- function(df){
  vars <- names(df)
  sapply(vars, function(v){
    rest <- setdiff(vars, v)
    if(length(rest)==0) return(1)
    
    # Handle cases where perfect multicollinearity exists
    tryCatch({
      fit <- lm(df[[v]] ~ ., data=df[,rest,drop=FALSE])
      R2  <- summary(fit)$r.squared
      
      # Cap R² at 0.9999 to avoid infinite VIF
      if(R2 >= 0.9999) R2 <- 0.9999
      
      1/(1-R2)
    }, error = function(e) {
      return(Inf)  # Return Inf if model fails
    })
  })
}

vif_all <- manual_vif(env_num)
cat("\nInitial VIF values (for reference):\n")
print(sort(vif_all, decreasing = TRUE))

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")
cat("Note: Being more flexible to include more meaningful variables\n")

# Don't pre-filter by VIF - let the selection process handle it
all_available_vars <- priority_df$Variable[priority_df$status == "available"]
all_deprioritized_vars <- priority_df$Variable[priority_df$status == "deprioritised"]

# Start with empty selection
selected <- character()

# Phase 1: Add available (prioritized) variables in order of Usnea correlation
cat("\n--- Phase 1: Adding prioritized variables ---\n")
for(v in all_available_vars){
  if(length(selected) == 0){
    # First variable - always add it
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    # Check if adding this variable creates acceptable VIF
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Also check if any existing variables would exceed VIF threshold
    max_existing_vif <- max(vifs[selected])
    
    if(new_var_vif <= vif_thr && max_existing_vif <= vif_thr){
      selected <- trial
      cat("Added (available):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    }
  }
}

# Phase 2: Try to add deprioritized variables that don't violate VIF
cat("\n--- Phase 2: Adding deprioritized variables where possible ---\n")
for(v in all_deprioritized_vars){
  if(length(selected) == 0){
    # If somehow no variables selected yet, add this one
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Check if any variable (new or existing) would exceed VIF threshold
    max_all_vif <- max(vifs)
    
    if(new_var_vif <= vif_thr && max_all_vif <= vif_thr){
      selected <- trial
      cat("Added (deprioritized):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    }
  }
}

# Phase 3: Final check - if we still have very few variables, try a different approach
cat("\n--- Phase 3: Final optimization ---\n")
if(length(selected) < 5){
  cat("Only", length(selected), "variables selected. Trying alternative approach...\n")
  
  # Try starting with the most important variables and build around them
  top_vars <- head(priority_df$Variable, 15)  # Top 15 by species correlation
  
  # Test different combinations starting from different points
  best_set <- selected
  best_count <- length(selected)
  
  for(i in 1:min(5, length(top_vars))){  # Try top 5 as starting points
    start_var <- top_vars[i]
    temp_selected <- start_var
    
    # Add other variables one by one, prioritizing by species correlation
    remaining_vars <- setdiff(top_vars, start_var)
    for(next_var in remaining_vars){
      trial <- unique(c(temp_selected, next_var))
      if(length(trial) > 1){
        vifs <- manual_vif(env_num[, trial, drop=FALSE])
        new_var_vif <- vifs[next_var]
        
        if(new_var_vif <= vif_thr){
          temp_selected <- trial
        }
      }
    }
    
    # If this combination is better, use it
    if(length(temp_selected) > best_count){
      best_set <- temp_selected
      best_count <- length(temp_selected)
    }
  }
  
  selected <- best_set
  cat("Alternative approach yielded", length(selected), "variables\n")
}

# 12. Final results
final_vifs <- if(length(selected) > 0)
  manual_vif(env_num[, selected, drop=FALSE]) else numeric()

results_tbl <- data.frame(
  Variable = selected,
  VIF = round(as.numeric(final_vifs), 3),
  Spearman_rho = round(cor_df_all$rho[match(selected, cor_df_all$Variable)], 3),
  Abs_Spearman_rho = round(abs(cor_df_all$rho[match(selected, cor_df_all$Variable)]), 3),
  Priority_Status = priority_df$status[match(selected, priority_df$Variable)],
  Var_Type = var_type_row[selected],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Abs_Spearman_rho))

cat("\nStep 4: Final selected variables (VIF < ", vif_thr, "):\n", sep="")
print(results_tbl)

# Show summary
cat("\nSummary:\n")
cat("Total variables selected:", nrow(results_tbl), "\n")
cat("Available variables selected:", sum(results_tbl$Priority_Status == "available"), "\n")
cat("Deprioritized variables selected:", sum(results_tbl$Priority_Status == "deprioritised"), "\n")

# 13. Save final dataset
id_cols <- setdiff(names(combined)[1:(env_start-1)], "IMAGE")
final_dataset <- combined %>%
  select(all_of(id_cols)) %>%
  bind_cols(env_num[, selected, drop=FALSE], env_cat) %>%
  mutate(Usnea = usnea)

write.csv(results_tbl, "Usnea_VIF_results_improved.csv", row.names=FALSE)
write.csv(final_dataset, "VIF_Usnea_FINAL_improved.csv", row.names=FALSE)
write.csv(priority_df, "Variable_Priority_Ranking_Usnea.csv", row.names=FALSE)

message("Done: retained ", nrow(results_tbl), " variables with VIF < ", vif_thr, ".")






# ─────────────────────────────────────────────────────────────────────────────
# DESCHAMPSIA: STEPWISE SPEARMAN + PRIORITY-BASED VIF SELECTION
# ─────────────────────────────────────────────────────────────────────────────

# 0. Clean workspace
rm(list=ls())
library(dplyr)

# 1. File paths – EDIT THESE
env_path     <- "Environmental_Variables.csv"
species_path <- "Abundance_Data.csv"

# 2. Read & label
raw           <- read.csv(env_path, header=FALSE, stringsAsFactors=FALSE)
var_names     <- as.character(raw[2,])
var_type_row  <- setNames(as.character(raw[1,]), var_names)
colnames(raw) <- var_names
dados         <- raw[-c(1,2), ]

# 3. Convert numeric‐looking columns
numflag <- sapply(dados, function(x)
  !all(is.na(suppressWarnings(as.numeric(x)))))
for(col in names(dados)[numflag]) {
  dados[[col]] <- as.numeric(dados[[col]])
}

# 4. Merge species data
sp        <- read.csv(species_path, stringsAsFactors=FALSE)
site_col  <- if("SITE_FINAL" %in% names(dados)) "SITE_FINAL" else
  grep("SITE", names(dados), ignore.case=TRUE, value=TRUE)[1]
combined  <- merge(dados, sp,
                   by.x=site_col, by.y="SITE_QUAD", all=FALSE)

# 5. Split predictors
env_start <- 15
env_end   <- ncol(combined) - 5
env_block <- combined[, env_start:env_end]
is_num    <- sapply(env_block, is.numeric)
env_num   <- env_block[, is_num, drop=FALSE]
env_cat   <- env_block[, !is_num, drop=FALSE]

# drop zero/infinite‐variance
keep_num <- sapply(env_num, function(x){
  v <- var(x, na.rm=TRUE); is.finite(v) && v>0
})
env_num <- env_num[, keep_num, drop=FALSE]

# exclude IMAGE & response from predictors
env_num <- env_num[, setdiff(names(env_num), c("IMAGE", "Deschampsia")), drop=FALSE]

# 6. Extract Deschampsia abundance
deschampsia <- combined$Deschampsia

# 7. Step 1 – Spearman among env vars (no species) - Threshold for microscale-microscale correlations
env_cor_thr <- 0.85  # Higher threshold for microscale-microscale correlations only
env_cor_mat <- cor(env_num, method="spearman", use="pairwise.complete.obs")
pairs_idx   <- which(upper.tri(env_cor_mat) & abs(env_cor_mat) > env_cor_thr,
                     arr.ind=TRUE)

if(nrow(pairs_idx)>0) {
  env_pairs <- data.frame(
    Var1 = rownames(env_cor_mat)[pairs_idx[,1]],
    Var2 = colnames(env_cor_mat)[pairs_idx[,2]],
    rho  = env_cor_mat[pairs_idx],
    stringsAsFactors=FALSE
  )
} else {
  env_pairs <- data.frame(
    Var1=character(), Var2=character(), rho=numeric(),
    stringsAsFactors=FALSE
  )
}

cat("\nStep 1: Highly inter-correlated environmental variables (|ρ| > ",
    env_cor_thr, "):\n", sep="")
print(env_pairs)

# 8. Step 2 – Spearman vs Deschampsia
spearman_all <- function(vars, df, y){
  out <- data.frame(Variable=character(), rho=numeric(), p=numeric(),
                    stringsAsFactors=FALSE)
  for(v in vars){
    x <- df[[v]]
    if(length(unique(na.omit(x)))>1 &&
       length(unique(na.omit(y)))>1){
      tst <- try(cor.test(x, y, method="spearman", exact=FALSE),
                 silent=TRUE)
      if(!inherits(tst, "try-error")){
        out <- rbind(out, data.frame(
          Variable=v,
          rho      = as.numeric(tst$estimate),
          p        = tst$p.value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  out$p_adj <- p.adjust(out$p, method="BH")
  out <- out[order(-abs(out$rho)), ]
  rownames(out) <- NULL
  return(out)
}

all_vars   <- names(env_num)
cor_df_all <- spearman_all(all_vars, env_num, deschampsia)

cat("\nStep 2: Spearman correlations vs Deschampsia:\n")
print(cor_df_all)

# 9. Step 3 – IMPROVED: Create comprehensive priority ranking
cat("\nStep 3: Creating comprehensive priority ranking...\n")

# Create a dataframe with all variables and their Deschampsia correlations
priority_df <- data.frame(
  Variable = cor_df_all$Variable,
  abs_rho = abs(cor_df_all$rho),
  rho = cor_df_all$rho,
  stringsAsFactors = FALSE
)

# Identify variables involved in high correlations
correlated_vars <- unique(c(env_pairs$Var1, env_pairs$Var2))

# For each variable, determine if it should be deprioritized
priority_df$status <- "available"
priority_df$reason <- ""

for(i in seq_len(nrow(env_pairs))){
  v1 <- env_pairs$Var1[i]
  v2 <- env_pairs$Var2[i]
  
  # Get absolute correlations with Deschampsia
  r1 <- priority_df$abs_rho[priority_df$Variable == v1]
  r2 <- priority_df$abs_rho[priority_df$Variable == v2]
  
  if(length(r1) > 0 && length(r2) > 0){
    # The variable with LOWER Deschampsia correlation gets deprioritised
    if(r1 > r2){
      priority_df$status[priority_df$Variable == v2] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Lower correlation with Deschampsia than ", v1)
    } else if(r2 > r1){
      priority_df$status[priority_df$Variable == v1] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v1] <- paste0("Lower correlation with Deschampsia than ", v2)
    } else {
      # If equal, deprioritize the second one (arbitrary but consistent)
      priority_df$status[priority_df$Variable == v2] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Equal correlation with Deschampsia as ", v1, " (tie-breaker)")
    }
  }
}

# Sort by priority: available variables first, then by absolute Deschampsia correlation
priority_df <- priority_df[order(priority_df$status, -priority_df$abs_rho), ]

cat("\nPriority ranking (top 15):\n")
print(head(priority_df, 15))

# 10. Compute initial VIF for all env predictors (for reference only)
manual_vif <- function(df){
  vars <- names(df)
  sapply(vars, function(v){
    rest <- setdiff(vars, v)
    if(length(rest)==0) return(1)
    
    # Handle cases where perfect multicollinearity exists
    tryCatch({
      fit <- lm(df[[v]] ~ ., data=df[,rest,drop=FALSE])
      R2  <- summary(fit)$r.squared
      
      # Cap R² at 0.9999 to avoid infinite VIF
      if(R2 >= 0.9999) R2 <- 0.9999
      
      1/(1-R2)
    }, error = function(e) {
      return(Inf)  # Return Inf if model fails
    })
  })
}

vif_all <- manual_vif(env_num)
cat("\nInitial VIF values (for reference):\n")
print(sort(vif_all, decreasing = TRUE))

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")
cat("Note: Being more flexible to include more meaningful variables\n")

# Don't pre-filter by VIF - let the selection process handle it
all_available_vars <- priority_df$Variable[priority_df$status == "available"]
all_deprioritized_vars <- priority_df$Variable[priority_df$status == "deprioritised"]

# Start with empty selection
selected <- character()

# Phase 1: Add available (prioritized) variables in order of Deschampsia correlation
cat("\n--- Phase 1: Adding prioritized variables ---\n")
for(v in all_available_vars){
  if(length(selected) == 0){
    # First variable - always add it
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    # Check if adding this variable creates acceptable VIF
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Also check if any existing variables would exceed VIF threshold
    max_existing_vif <- max(vifs[selected])
    
    if(new_var_vif <= vif_thr && max_existing_vif <= vif_thr){
      selected <- trial
      cat("Added (available):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    }
  }
}

# Phase 2: Try to add deprioritized variables that don't violate VIF
cat("\n--- Phase 2: Adding deprioritized variables where possible ---\n")
for(v in all_deprioritized_vars){
  if(length(selected) == 0){
    # If somehow no variables selected yet, add this one
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Check if any variable (new or existing) would exceed VIF threshold
    max_all_vif <- max(vifs)
    
    if(new_var_vif <= vif_thr && max_all_vif <= vif_thr){
      selected <- trial
      cat("Added (deprioritized):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    }
  }
}

# Phase 3: Final check - if we still have very few variables, try a different approach
cat("\n--- Phase 3: Final optimization ---\n")
if(length(selected) < 5){
  cat("Only", length(selected), "variables selected. Trying alternative approach...\n")
  
  # Try starting with the most important variables and build around them
  top_vars <- head(priority_df$Variable, 15)  # Top 15 by species correlation
  
  # Test different combinations starting from different points
  best_set <- selected
  best_count <- length(selected)
  
  for(i in 1:min(5, length(top_vars))){  # Try top 5 as starting points
    start_var <- top_vars[i]
    temp_selected <- start_var
    
    # Add other variables one by one, prioritizing by species correlation
    remaining_vars <- setdiff(top_vars, start_var)
    for(next_var in remaining_vars){
      trial <- unique(c(temp_selected, next_var))
      if(length(trial) > 1){
        vifs <- manual_vif(env_num[, trial, drop=FALSE])
        new_var_vif <- vifs[next_var]
        
        if(new_var_vif <= vif_thr){
          temp_selected <- trial
        }
      }
    }
    
    # If this combination is better, use it
    if(length(temp_selected) > best_count){
      best_set <- temp_selected
      best_count <- length(temp_selected)
    }
  }
  
  selected <- best_set
  cat("Alternative approach yielded", length(selected), "variables\n")
}

# 12. Final results
final_vifs <- if(length(selected) > 0)
  manual_vif(env_num[, selected, drop=FALSE]) else numeric()

results_tbl <- data.frame(
  Variable = selected,
  VIF = round(as.numeric(final_vifs), 3),
  Spearman_rho = round(cor_df_all$rho[match(selected, cor_df_all$Variable)], 3),
  Abs_Spearman_rho = round(abs(cor_df_all$rho[match(selected, cor_df_all$Variable)]), 3),
  Priority_Status = priority_df$status[match(selected, priority_df$Variable)],
  Var_Type = var_type_row[selected],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Abs_Spearman_rho))

cat("\nStep 4: Final selected variables (VIF < ", vif_thr, "):\n", sep="")
print(results_tbl)

# Show summary
cat("\nSummary:\n")
cat("Total variables selected:", nrow(results_tbl), "\n")
cat("Available variables selected:", sum(results_tbl$Priority_Status == "available"), "\n")
cat("Deprioritized variables selected:", sum(results_tbl$Priority_Status == "deprioritised"), "\n")

# 13. Save final dataset
id_cols <- setdiff(names(combined)[1:(env_start-1)], "IMAGE")
final_dataset <- combined %>%
  select(all_of(id_cols)) %>%
  bind_cols(env_num[, selected, drop=FALSE], env_cat) %>%
  mutate(Deschampsia = deschampsia)

write.csv(results_tbl, "Deschampsia_VIF_results_improved.csv", row.names=FALSE)
write.csv(final_dataset, "VIF_Deschampsia_FINAL_improved.csv", row.names=FALSE)
write.csv(priority_df, "Variable_Priority_Ranking_Deschampsia.csv", row.names=FALSE)

message("Done: retained ", nrow(results_tbl), " variables with VIF < ", vif_thr, ".")



# ─────────────────────────────────────────────────────────────────────────────
# COLOBANTHUS: IMPROVED STEPWISE SPEARMAN + PRIORITY-BASED VIF SELECTION
# ─────────────────────────────────────────────────────────────────────────────

# 0. Clean workspace
rm(list=ls())
library(dplyr)

# 1. File paths 
env_path     <- "Environmental_Variables.csv"
species_path <- "Abundance_Data.csv"

# 2. Read & label
raw           <- read.csv(env_path, header=FALSE, stringsAsFactors=FALSE)
var_names     <- as.character(raw[2,])
var_type_row  <- setNames(as.character(raw[1,]), var_names)
colnames(raw) <- var_names
dados         <- raw[-c(1,2), ]

# 3. Convert numeric‐looking columns
numflag <- sapply(dados, function(x)
  !all(is.na(suppressWarnings(as.numeric(x)))))
for(col in names(dados)[numflag]) {
  dados[[col]] <- as.numeric(dados[[col]])
}

# 4. Merge species data
sp        <- read.csv(species_path, stringsAsFactors=FALSE)
site_col  <- if("SITE_FINAL" %in% names(dados)) "SITE_FINAL" else
  grep("SITE", names(dados), ignore.case=TRUE, value=TRUE)[1]
combined  <- merge(dados, sp,
                   by.x=site_col, by.y="SITE_QUAD", all=FALSE)

# 5. Split predictors
env_start <- 15
env_end   <- ncol(combined) - 5
env_block <- combined[, env_start:env_end]
is_num    <- sapply(env_block, is.numeric)
env_num   <- env_block[, is_num, drop=FALSE]
env_cat   <- env_block[, !is_num, drop=FALSE]

# drop zero/infinite‐variance
keep_num <- sapply(env_num, function(x){
  v <- var(x, na.rm=TRUE); is.finite(v) && v>0
})
env_num <- env_num[, keep_num, drop=FALSE]

# exclude IMAGE & response from predictors
env_num <- env_num[, setdiff(names(env_num), c("IMAGE", "Colobanthus")), drop=FALSE]

# 6. Extract Colobanthus abundance
colobanthus <- combined$Colobanthus

# 7. Step 1 – Spearman among env vars (no species) - Threshold for microscale-microscale correlations
env_cor_thr <- 0.85  # Higher threshold for microscale-microscale correlations only
env_cor_mat <- cor(env_num, method="spearman", use="pairwise.complete.obs")
pairs_idx   <- which(upper.tri(env_cor_mat) & abs(env_cor_mat) > env_cor_thr,
                     arr.ind=TRUE)

if(nrow(pairs_idx)>0) {
  env_pairs <- data.frame(
    Var1 = rownames(env_cor_mat)[pairs_idx[,1]],
    Var2 = colnames(env_cor_mat)[pairs_idx[,2]],
    rho  = env_cor_mat[pairs_idx],
    stringsAsFactors=FALSE
  )
} else {
  env_pairs <- data.frame(
    Var1=character(), Var2=character(), rho=numeric(),
    stringsAsFactors=FALSE
  )
}

cat("\nStep 1: Highly inter-correlated environmental variables (|ρ| > ",
    env_cor_thr, "):\n", sep="")
print(env_pairs)

# 8. Step 2 – Spearman vs Colobanthus
spearman_all <- function(vars, df, y){
  out <- data.frame(Variable=character(), rho=numeric(), p=numeric(),
                    stringsAsFactors=FALSE)
  for(v in vars){
    x <- df[[v]]
    if(length(unique(na.omit(x)))>1 &&
       length(unique(na.omit(y)))>1){
      tst <- try(cor.test(x, y, method="spearman", exact=FALSE),
                 silent=TRUE)
      if(!inherits(tst, "try-error")){
        out <- rbind(out, data.frame(
          Variable=v,
          rho      = as.numeric(tst$estimate),
          p        = tst$p.value,
          stringsAsFactors=FALSE
        ))
      }
    }
  }
  out$p_adj <- p.adjust(out$p, method="BH")
  out <- out[order(-abs(out$rho)), ]
  rownames(out) <- NULL
  return(out)
}

all_vars   <- names(env_num)
cor_df_all <- spearman_all(all_vars, env_num, colobanthus)

cat("\nStep 2: Spearman correlations vs Colobanthus:\n")
print(cor_df_all)

# 9. Step 3 – IMPROVED: Create comprehensive priority ranking
cat("\nStep 3: Creating comprehensive priority ranking...\n")

# Create a dataframe with all variables and their Colobanthus correlations
priority_df <- data.frame(
  Variable = cor_df_all$Variable,
  abs_rho = abs(cor_df_all$rho),
  rho = cor_df_all$rho,
  stringsAsFactors = FALSE
)

# Identify variables involved in high correlations
correlated_vars <- unique(c(env_pairs$Var1, env_pairs$Var2))

# For each variable, determine if it should be deprioritized
priority_df$status <- "available"
priority_df$reason <- ""

for(i in seq_len(nrow(env_pairs))){
  v1 <- env_pairs$Var1[i]
  v2 <- env_pairs$Var2[i]
  
  # Get absolute correlations with Colobanthus
  r1 <- priority_df$abs_rho[priority_df$Variable == v1]
  r2 <- priority_df$abs_rho[priority_df$Variable == v2]
  
  if(length(r1) > 0 && length(r2) > 0){
    # The variable with LOWER Colobanthus correlation gets deprioritized
    if(r1 > r2){
      priority_df$status[priority_df$Variable == v2] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Lower correlation with Colobanthus than ", v1)
    } else if(r2 > r1){
      priority_df$status[priority_df$Variable == v1] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v1] <- paste0("Lower correlation with Colobanthus than ", v2)
    } else {
      # If equal, deprioritize the second one (arbitrary but consistent)
      priority_df$status[priority_df$Variable == v2] <- "deprioritised"
      priority_df$reason[priority_df$Variable == v2] <- paste0("Equal correlation with Colobanthus as ", v1, " (tie-breaker)")
    }
  }
}

# Sort by priority: available variables first, then by absolute Colobanthus correlation
priority_df <- priority_df[order(priority_df$status, -priority_df$abs_rho), ]

cat("\nPriority ranking (top 15):\n")
print(head(priority_df, 15))

# 10. Compute initial VIF for all env predictors (for reference only)
manual_vif <- function(df){
  vars <- names(df)
  sapply(vars, function(v){
    rest <- setdiff(vars, v)
    if(length(rest)==0) return(1)
    
    # Handle cases where perfect multicollinearity exists
    tryCatch({
      fit <- lm(df[[v]] ~ ., data=df[,rest,drop=FALSE])
      R2  <- summary(fit)$r.squared
      
      # Cap R² at 0.9999 to avoid infinite VIF
      if(R2 >= 0.9999) R2 <- 0.9999
      
      1/(1-R2)
    }, error = function(e) {
      return(Inf)  # Return Inf if model fails
    })
  })
}

vif_all <- manual_vif(env_num)
cat("\nInitial VIF values (for reference):\n")
print(sort(vif_all, decreasing = TRUE))

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")

# 11. Step 4 – FLEXIBLE: Include more variables with smart priority-based selection
vif_thr <- 5
cat("\nStep 4: Flexible VIF selection with priority consideration...\n")
cat("Note: Being more flexible to include more meaningful variables\n")

# Don't pre-filter by VIF - let the selection process handle it
all_available_vars <- priority_df$Variable[priority_df$status == "available"]
all_deprioritized_vars <- priority_df$Variable[priority_df$status == "deprioritized"]

# Start with empty selection
selected <- character()

# Phase 1: Add available (prioritised) variables in order of Colobanthus correlation
cat("\n--- Phase 1: Adding prioritized variables ---\n")
for(v in all_available_vars){
  if(length(selected) == 0){
    # First variable - always add it
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    # Check if adding this variable creates acceptable VIF
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Also check if any existing variables would exceed VIF threshold
    max_existing_vif <- max(vifs[selected])
    
    if(new_var_vif <= vif_thr && max_existing_vif <= vif_thr){
      selected <- trial
      cat("Added (available):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max existing VIF:", round(max_existing_vif, 2), "\n")
    }
  }
}

# Phase 2: Try to add deprioritised variables that don't violate VIF
cat("\n--- Phase 2: Adding deprioritized variables where possible ---\n")
for(v in all_deprioritized_vars){
  if(length(selected) == 0){
    # If somehow no variables selected yet, add this one
    selected <- v
    cat("Added (first variable):", v, "- No VIF constraint for single variable\n")
  } else {
    trial <- unique(c(selected, v))
    vifs <- manual_vif(env_num[, trial, drop=FALSE])
    new_var_vif <- vifs[v]
    
    # Check if any variable (new or existing) would exceed VIF threshold
    max_all_vif <- max(vifs)
    
    if(new_var_vif <= vif_thr && max_all_vif <= vif_thr){
      selected <- trial
      cat("Added (deprioritized):", v, "- VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    } else {
      cat("Skipped (VIF constraint violated):", v, "- New VIF:", round(new_var_vif, 2), 
          "- Max all VIF:", round(max_all_vif, 2), "\n")
    }
  }
}

# Phase 3: Final check - if we still have very few variables, try a different approach
cat("\n--- Phase 3: Final optimization ---\n")
if(length(selected) < 5){
  cat("Only", length(selected), "variables selected. Trying alternative approach...\n")
  
  # Try starting with the most important variables and build around them
  top_vars <- head(priority_df$Variable, 15)  # Top 15 by species correlation
  
  # Test different combinations starting from different points
  best_set <- selected
  best_count <- length(selected)
  
  for(i in 1:min(5, length(top_vars))){  # Try top 5 as starting points
    start_var <- top_vars[i]
    temp_selected <- start_var
    
    # Add other variables one by one, prioritizing by species correlation
    remaining_vars <- setdiff(top_vars, start_var)
    for(next_var in remaining_vars){
      trial <- unique(c(temp_selected, next_var))
      if(length(trial) > 1){
        vifs <- manual_vif(env_num[, trial, drop=FALSE])
        new_var_vif <- vifs[next_var]
        
        if(new_var_vif <= vif_thr){
          temp_selected <- trial
        }
      }
    }
    
    # If this combination is better, use it
    if(length(temp_selected) > best_count){
      best_set <- temp_selected
      best_count <- length(temp_selected)
    }
  }
  
  selected <- best_set
  cat("Alternative approach yielded", length(selected), "variables\n")
}

# 12. Final results
final_vifs <- if(length(selected) > 0)
  manual_vif(env_num[, selected, drop=FALSE]) else numeric()

results_tbl <- data.frame(
  Variable = selected,
  VIF = round(as.numeric(final_vifs), 3),
  Spearman_rho = round(cor_df_all$rho[match(selected, cor_df_all$Variable)], 3),
  Abs_Spearman_rho = round(abs(cor_df_all$rho[match(selected, cor_df_all$Variable)]), 3),
  Priority_Status = priority_df$status[match(selected, priority_df$Variable)],
  Var_Type = var_type_row[selected],
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Abs_Spearman_rho))

cat("\nStep 4: Final selected variables (VIF < ", vif_thr, "):\n", sep="")
print(results_tbl)

# Show summary
cat("\nSummary:\n")
cat("Total variables selected:", nrow(results_tbl), "\n")
cat("Available variables selected:", sum(results_tbl$Priority_Status == "available"), "\n")
cat("Deprioritized variables selected:", sum(results_tbl$Priority_Status == "deprioritised"), "\n")

# 13. Save final dataset
id_cols <- setdiff(names(combined)[1:(env_start-1)], "IMAGE")
final_dataset <- combined %>%
  select(all_of(id_cols)) %>%
  bind_cols(env_num[, selected, drop=FALSE], env_cat) %>%
  mutate(Colobanthus = colobanthus)

write.csv(results_tbl, "Colobanthus_VIF_results_improved.csv", row.names=FALSE)
write.csv(final_dataset, "VIF_Colobanthus_FINAL_improved.csv", row.names=FALSE)
write.csv(priority_df, "Variable_Priority_Ranking_Colobanthus.csv", row.names=FALSE)

message("Done: retained ", nrow(results_tbl), " variables with VIF < ", vif_thr, ".")



############################# Combine the CSV files #########################################################
library(dplyr)

# Define the species names and corresponding file names
species_files <- c(
  "Colobanthus" = "VIF_Colobanthus_FINAL_improved.csv",
  "Deschampsia" = "VIF_Deschampsia_FINAL_improved.csv", 
  "Usnea" = "VIF_Usnea_FINAL_improved.csv",
  "Sanionia" = "VIF_Sanionia_FINAL_improved.csv"
)

# Read all files first
species_data_list <- list()

for (species in names(species_files)) {
  file_path <- species_files[species]
  
  if (file.exists(file_path)) {
    temp_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    species_data_list[[species]] <- temp_data
    cat("Read:", file_path, "- Rows:", nrow(temp_data), "- Columns:", ncol(temp_data), "\n")
  } else {
    cat("Warning: File not found:", file_path, "\n")
  }
}

if (length(species_data_list) > 0) {
  cat("\n=== MERGING WITH DUPLICATE COLUMN HANDLING ===\n")
  
  # Start with first dataset
  combined_data <- species_data_list[[1]]
  cat("Starting with", names(species_data_list)[1], "\n")
  
  # Merge each subsequent dataset
  for (i in 2:length(species_data_list)) {
    species_name <- names(species_data_list)[i]
    next_data <- species_data_list[[i]]
    
    cat("Merging", species_name, "...\n")
    
    # Find common columns (excluding SITE_FINAL)
    common_cols <- intersect(colnames(combined_data), colnames(next_data))
    common_cols <- common_cols[common_cols != "SITE_FINAL"]
    
    if (length(common_cols) > 0) {
      cat("  Common columns found:", paste(common_cols, collapse = ", "), "\n")
      
      # Remove common columns from the dataset being merged
      cols_to_keep <- c("SITE_FINAL", setdiff(colnames(next_data), common_cols))
      next_data <- next_data[, cols_to_keep, drop = FALSE]
      cat("  Keeping unique columns from", species_name, ":", 
          paste(setdiff(colnames(next_data), "SITE_FINAL"), collapse = ", "), "\n")
    }
    
    # Merge the datasets
    combined_data <- merge(combined_data, next_data, by = "SITE_FINAL", all = TRUE)
    
    cat("  After merge - Rows:", nrow(combined_data), "- Columns:", ncol(combined_data), "\n")
  }
  
  # Display final summary
  cat("\n=== FINAL DATASET SUMMARY ===\n")
  cat("Total rows:", nrow(combined_data), "\n")
  cat("Total columns:", ncol(combined_data), "\n")
  cat("Unique sites:", length(unique(combined_data$SITE_FINAL)), "\n")
  
  # Save the result
  write.csv(combined_data, "VIF_ALL_SPECIES_FINAL_combined.csv", row.names = FALSE)
  cat("\nCombined dataset saved as: VIF_ALL_SPECIES_FINAL_combined.csv\n")
  
  # Show all column names
  cat("\nAll columns in final dataset:\n")
  print(colnames(combined_data))
  
  # Check for missing data
  cat("\n=== MISSING DATA CHECK ===\n")
  missing_counts <- colSums(is.na(combined_data))
  cols_with_missing <- missing_counts[missing_counts > 0]
  
  if (length(cols_with_missing) > 0) {
    cat("Columns with missing data:\n")
    print(cols_with_missing)
  } else {
    cat("No missing data!\n")
  }
  
  # Show first few rows
  cat("\nFirst 3 rows (first 8 columns):\n")
  print(combined_data[1:3, 1:min(8, ncol(combined_data))])
  
} else {
  cat("Error: No files were read successfully.\n")
}