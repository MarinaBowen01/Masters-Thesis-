### Current predictions code for each species 
### All code have been set up to be used on a HPC 



######## Colobanthus predictions using the HPC ##########

library(terra)  
library(sf)     
library(dplyr)
library(randomForest)
library(quantregForest)

cat("=== COLOBANTHUS STRICT PROCESSING ===\n")

#### CONFIGURATION ####
# Extended predictor variables for Colobanthus - NO DUMMY DATA ALLOWED
predictors_presence <- c("VDEP", "Northing", "LSF", "Easting", "FCF", "SVF", "MPI", "DISTCOAST", "TXT")
predictors_abundance <- c("DISTCOAST", "VDEP", "TPW", "FLAC", "ASP")

base_dir <- "/rds/general/user/mrb24/ephemeral"  
output_dir <- "/rds/general/user/mrb24/ephemeral"   

# Create output subdirectory
colobanthus_output_dir <- file.path(output_dir, "Colobanthus_AllFiles")
if (!dir.exists(colobanthus_output_dir)) {
  dir.create(colobanthus_output_dir, recursive = TRUE)
  cat("Created output directory:", colobanthus_output_dir, "\n")
}

terraOptions(memfrac = 0.3, tempdir = tempdir())
options(warn = 1)

#### FUNCTIONS ####

# Function to find all DEM files
find_all_dem_files <- function(base_path) {
  dem_dir <- file.path(base_path, "DEMICEFREE")
  if (!dir.exists(dem_dir)) {
    stop(paste("DEM directory not found:", dem_dir))
  }
  
  all_dem_files <- list.files(dem_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  
  if (length(all_dem_files) == 0) {
    stop("No DEM TIF files found")
  }
  
  # Sort files naturally
  all_dem_files <- sort(all_dem_files)
  
  cat(paste("Found", length(all_dem_files), "DEM files to process\n"))
  return(all_dem_files)
}

# Function to find matching file by spatial extent with minimum overlap threshold
find_matching_file_by_extent <- function(reference_file, search_dirs, patterns, min_overlap = 0.8) {
  tryCatch({
    ref_rast <- rast(reference_file)
    ref_ext <- ext(ref_rast)
    
    for (i in seq_along(search_dirs)) {
      search_dir <- search_dirs[i]
      pattern <- patterns[i]
      
      if (!dir.exists(search_dir)) {
        next
      }
      
      # Find all files matching pattern
      all_files <- list.files(search_dir, pattern = paste0(pattern, ".*\\.tif$"), 
                              full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      
      # Test each file for spatial overlap
      for (file in all_files) {
        tryCatch({
          test_rast <- rast(file)
          test_ext <- ext(test_rast)
          
          # Calculate overlap
          ref_vec <- as.vector(ref_ext)
          test_vec <- as.vector(test_ext)
          
          x_overlap <- max(0, min(ref_vec[2], test_vec[2]) - max(ref_vec[1], test_vec[1]))
          y_overlap <- max(0, min(ref_vec[4], test_vec[4]) - max(ref_vec[3], test_vec[3]))
          
          if (x_overlap > 0 && y_overlap > 0) {
            overlap_area <- x_overlap * y_overlap
            ref_area <- (ref_vec[2] - ref_vec[1]) * (ref_vec[4] - ref_vec[3])
            overlap_pct <- overlap_area / ref_area
            
            if (overlap_pct >= min_overlap) {
              return(file)
            }
          }
        }, error = function(e) {
          # Continue to next file
        })
      }
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Function to find all required matching files - STRICT MODE (no dummy data)
find_required_files_for_dem <- function(dem_file, base_path) {
  cat(paste("Finding required files for:", basename(dem_file), "\n"))
  
  # Define search directories and patterns for all required variables
  search_configs <- list(
    vdep = list(dirs = c(file.path(base_path, "VDEP")), patterns = c("VDEP", "VDEPICEFREE")),
    lsf = list(dirs = c(base_path), patterns = c("LSF", "LSFICEFREE")),
    asp = list(dirs = c(base_path), patterns = c("ASP", "ASPICEFREE")),
    tpw = list(dirs = c(base_path), patterns = c("TPW", "TPWICEFREE")),
    flac = list(dirs = c(base_path), patterns = c("FLAC", "FLACICEFREE")),
    svf = list(dirs = c(base_path), patterns = c("SVF", "SVFICEFREE")),
    mpi = list(dirs = c(base_path), patterns = c("MPI", "MPIICEFREE")),
    txt = list(dirs = c(base_path), patterns = c("TXT", "TXTICEFREE"))
  )
  
  matches <- list()
  
  for (var_name in names(search_configs)) {
    config <- search_configs[[var_name]]
    
    cat(paste("  Searching for", toupper(var_name), "..."))
    
    match_file <- find_matching_file_by_extent(
      dem_file, 
      config$dirs, 
      config$patterns, 
      min_overlap = 0.7
    )
    
    if (!is.null(match_file)) {
      matches[[var_name]] <- match_file
      cat(paste(" ✓", basename(match_file), "\n"))
    } else {
      cat(" ✗ NOT FOUND\n")
    }
  }
  
  return(matches)
}

# Strict check for all required variables - NO DUMMY DATA ALLOWED
check_strict_requirements <- function(matches) {
  # All presence predictors except coordinates and FCF (handled separately) and DISTCOAST (calculated)
  required_presence_vars <- c("vdep", "lsf", "svf", "mpi", "txt")
  # All abundance predictors except coordinates and DISTCOAST (calculated)
  required_abundance_vars <- c("vdep", "tpw", "flac", "asp")
  
  # Combined requirements
  all_required <- unique(c(required_presence_vars, required_abundance_vars))
  
  missing <- all_required[!all_required %in% names(matches)]
  
  if (length(missing) > 0) {
    cat(paste("  ✗ Missing required variables:", paste(toupper(missing), collapse = ", "), "\n"))
    return(FALSE)
  }
  
  cat("  ✓ All required variables found\n")
  return(TRUE)
}

# Proper distance to coast calculation
calculate_distance_to_coast_proper <- function(raster_template, coastline_path) {
  cat("  Calculating proper distance to coast...\n")
  
  if (!file.exists(coastline_path)) {
    cat("  ✗ Coastline file not found, cannot calculate proper distances\n")
    return(NULL)
  }
  
  tryCatch({
    # Read coastline
    coastline <- st_read(coastline_path, quiet = TRUE)
    coastline <- st_transform(coastline, crs(raster_template))
    
    # Get valid pixel coordinates
    coords <- xyFromCell(raster_template, which(!is.na(values(raster_template))))
    
    if (nrow(coords) == 0) {
      cat("  ✗ No valid pixels found\n")
      return(NULL)
    }
    
    # Convert to sf points
    points <- st_as_sf(data.frame(coords), coords = c("x", "y"), crs = crs(raster_template))
    
    # Calculate minimum distance to coastline for each point
    cat(paste("  Computing distances for", nrow(points), "points...\n"))
    distances <- st_distance(points, coastline)
    min_distances <- apply(distances, 1, min)
    
    # Create distance raster
    dist_raster <- raster_template
    values(dist_raster) <- NA
    valid_cells <- which(!is.na(values(raster_template)))
    values(dist_raster)[valid_cells] <- as.numeric(min_distances)
    
    names(dist_raster) <- "DISTCOAST"
    cat("  ✓ Distance to coast calculated\n")
    return(dist_raster)
    
  }, error = function(e) {
    cat(paste("  ✗ Error calculating distance to coast:", e$message, "\n"))
    return(NULL)
  })
}

# Function to convert aspect to categorical
convert_aspect_to_categorical <- function(aspect_raster) {
  cat("  Converting aspect to categorical\n")
  
  aspect_vals <- values(aspect_raster, mat = FALSE)
  aspect_cats <- rep(NA_character_, length(aspect_vals))
  
  valid_vals <- !is.na(aspect_vals)
  aspect_cats[valid_vals & (aspect_vals >= 337.5 | aspect_vals < 22.5)] <- "N"
  aspect_cats[valid_vals & aspect_vals >= 22.5 & aspect_vals < 67.5] <- "NE"
  aspect_cats[valid_vals & aspect_vals >= 67.5 & aspect_vals < 112.5] <- "E"
  aspect_cats[valid_vals & aspect_vals >= 112.5 & aspect_vals < 157.5] <- "SE"
  aspect_cats[valid_vals & aspect_vals >= 157.5 & aspect_vals < 202.5] <- "S"
  aspect_cats[valid_vals & aspect_vals >= 202.5 & aspect_vals < 247.5] <- "SW"
  aspect_cats[valid_vals & aspect_vals >= 247.5 & aspect_vals < 292.5] <- "W"
  aspect_cats[valid_vals & aspect_vals >= 292.5 & aspect_vals < 337.5] <- "NW"
  
  aspect_categorical <- aspect_raster
  values(aspect_categorical) <- aspect_cats
  
  return(aspect_categorical)
}

# Calculate coordinates
calculate_coordinates <- function(raster_template) {
  coords <- xyFromCell(raster_template, 1:ncell(raster_template))
  easting_raster <- raster_template
  northing_raster <- raster_template
  values(easting_raster) <- coords[, 1]
  values(northing_raster) <- coords[, 2]
  names(easting_raster) <- "Easting"
  names(northing_raster) <- "Northing"
  return(c(easting_raster, northing_raster))
}

# Clean prediction data - strict mode
clean_prediction_data_strict <- function(prediction_data, predictors) {
  cat("  Strict data cleaning (no NA tolerance)...\n")
  
  # Remove rows with any NA in predictor variables
  complete_rows <- complete.cases(prediction_data[, predictors, drop = FALSE])
  n_removed <- sum(!complete_rows)
  
  if (n_removed > 0) {
    cat(paste("  Removing", n_removed, "rows with NA values\n"))
    prediction_data <- prediction_data[complete_rows, , drop = FALSE]
  }
  
  if (nrow(prediction_data) == 0) {
    cat("  ✗ No complete cases remaining\n")
    return(NULL)
  }
  
  # Handle factor variables
  if ("ASP" %in% predictors && "ASP" %in% names(prediction_data)) {
    prediction_data$ASP <- factor(as.character(prediction_data$ASP), 
                                  levels = c("FLAT", "N", "NE", "E", "SE", "S", "SW", "W", "NW"))
  }
  
  cat(paste("  Clean data:", nrow(prediction_data), "rows\n"))
  return(prediction_data)
}

# Process single DEM file - STRICT MODE
process_single_dem_strict <- function(dem_file, matches, coastline_path, fcf_file, rf_model, qrf_model) {
  tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
  cat(paste("Processing:", tile_id, "\n"))
  
  tryCatch({
    # Load DEM as template
    altitude_raster <- rast(dem_file)
    names(altitude_raster) <- "ALTITUDE"
    
    if (sum(!is.na(values(altitude_raster))) == 0) {
      cat("  ✗ No valid altitude data\n")
      return(NULL)
    }
    
    template_raster <- altitude_raster
    all_rasters <- list()
    
    # Load FCF layer
    cat("  Loading FCF layer...\n")
    fcf_raster <- rast(fcf_file)
    names(fcf_raster) <- "FCF"
    
    if (!compareGeom(fcf_raster, template_raster, stopOnError = FALSE)) {
      fcf_raster <- resample(fcf_raster, template_raster, method = "bilinear")
    }
    all_rasters$FCF <- fcf_raster
    
    # Load all required variables - NO DUMMY DATA
    variable_map <- list(
      VDEP = "vdep", LSF = "lsf", ASP = "asp", TPW = "tpw", 
      FLAC = "flac", SVF = "svf", MPI = "mpi", TXT = "txt"
    )
    
    for (var_name in names(variable_map)) {
      match_key <- variable_map[[var_name]]
      
      if (match_key %in% names(matches)) {
        cat(paste("  Loading", var_name, "from", basename(matches[[match_key]]), "\n"))
        
        var_raster <- rast(matches[[match_key]])
        names(var_raster) <- var_name
        
        if (!compareGeom(var_raster, template_raster, stopOnError = FALSE)) {
          var_raster <- resample(var_raster, template_raster, method = "bilinear")
        }
        
        # Special handling for ASP (convert to categorical)
        if (var_name == "ASP") {
          var_raster <- convert_aspect_to_categorical(var_raster)
        }
        
        all_rasters[[var_name]] <- var_raster
      } else {
        cat(paste("  ✗ Missing required variable:", var_name, "\n"))
        return(NULL)
      }
    }
    
    # Calculate coordinates
    cat("  Calculating coordinates...\n")
    coord_rasters <- calculate_coordinates(template_raster)
    all_rasters$Easting <- coord_rasters[[1]]
    all_rasters$Northing <- coord_rasters[[2]]
    
    # Calculate proper distance to coast
    cat("  Calculating distance to coast...\n")
    dist_raster <- calculate_distance_to_coast_proper(template_raster, coastline_path)
    
    if (is.null(dist_raster)) {
      cat("  ✗ Could not calculate distance to coast\n")
      return(NULL)
    }
    
    all_rasters$DISTCOAST <- dist_raster
    
    # Convert to data frame
    cat("  Converting to data frame...\n")
    combined_rasters <- rast(all_rasters)
    prediction_data <- as.data.frame(combined_rasters, xy = TRUE, na.rm = TRUE)
    
    if (nrow(prediction_data) == 0) {
      cat("  ✗ No valid data after conversion\n")
      return(NULL)
    }
    
    # Clean data strictly
    presence_data <- clean_prediction_data_strict(prediction_data, predictors_presence)
    if (is.null(presence_data)) return(NULL)
    
    # Make predictions
    cat(paste("  Making predictions on", nrow(presence_data), "pixels\n"))
    
    prob_pred <- predict(rf_model, newdata = presence_data, type = "prob")[, "1"]
    binary_pred <- ifelse(prob_pred >= 0.5, 1, 0)
    
    # Abundance predictions for present pixels
    abundance_pred <- rep(NA, nrow(presence_data))
    present_idx <- which(binary_pred == 1)
    
    if (length(present_idx) > 0) {
      cat(paste("  Calculating abundance for", length(present_idx), "present pixels\n"))
      abundance_data <- clean_prediction_data_strict(presence_data[present_idx, ], predictors_abundance)
      
      if (!is.null(abundance_data) && nrow(abundance_data) > 0) {
        abundance_predictions <- predict(qrf_model, newdata = abundance_data, what = 0.5)
        abundance_pred[present_idx[1:length(abundance_predictions)]] <- pmax(pmin(abundance_predictions, 100), 0)
      }
    }
    
    # Create output rasters
    prob_raster <- binary_raster <- abundance_raster <- template_raster
    values(prob_raster) <- values(binary_raster) <- values(abundance_raster) <- NA
    
    # Map predictions back to raster
    cell_nums <- cellFromXY(template_raster, presence_data[, c("x", "y")])
    valid_cells <- !is.na(cell_nums)
    
    if (sum(valid_cells) > 0) {
      values(prob_raster)[cell_nums[valid_cells]] <- prob_pred[valid_cells]
      values(binary_raster)[cell_nums[valid_cells]] <- binary_pred[valid_cells]
      values(abundance_raster)[cell_nums[valid_cells]] <- abundance_pred[valid_cells]
    }
    
    names(prob_raster) <- "Presence_Probability"
    names(binary_raster) <- "Presence_Binary"
    names(abundance_raster) <- "Abundance_Prediction"
    
    return(list(
      tile_id = tile_id,
      presence_prob = prob_raster,
      presence_binary = binary_raster,
      abundance = abundance_raster,
      stats = list(
        pixels_present = sum(binary_pred == 1, na.rm = TRUE),
        total_pixels = nrow(presence_data),
        mean_probability = round(mean(prob_pred, na.rm = TRUE), 3)
      )
    ))
    
  }, error = function(e) {
    cat(paste("✗ ERROR:", e$message, "\n"))
    return(NULL)
  })
}

#### LOAD MODELS ####
train_models <- function() {
  training_file <- file.path(base_dir, "hurdle_training_data_list_Colobanthus_top_predictors.rds")
  if (!file.exists(training_file)) stop("Training data file not found!")
  
  training_data_list <- readRDS(training_file)
  train_data <- training_data_list[[1]]
  
  cat("Training models with extended predictors...\n")
  
  rf_clf <- randomForest(
    x = train_data[, predictors_presence, drop = FALSE],
    y = as.factor(train_data$Presence),
    ntree = 300,
    importance = TRUE
  )
  
  present_data <- train_data[train_data$Presence == "1" | train_data$Presence == 1, ]
  qrf <- quantregForest(
    x = present_data[, predictors_abundance, drop = FALSE],
    y = as.numeric(present_data$Colobanthus),
    ntree = 300
  )
  
  return(list(rf_classifier = rf_clf, qrf_model = qrf))
}

#### MAIN PROCESSING ####
main_processing <- function() {
  cat("Loading models...\n")
  models <- train_models()
  rf_model <- models$rf_classifier
  qrf_model <- models$qrf_model
  
  # Find all DEM files
  dem_files <- find_all_dem_files(base_dir)
  
  # FCF file
  fcf_file <- file.path(base_dir, "FCF", "Colobanthus_FCF_Climate_Envelope_Values.tif")
  if (!file.exists(fcf_file)) {
    stop("FCF file not found")
  }
  
  # Coastline file
  coastline_file <- file.path(base_dir, "COAST", "add_coastline_high_res_line_v7_10.shp")
  if (!file.exists(coastline_file)) {
    stop("Coastline file not found")
  }
  
  successful_tiles <- c()
  failed_tiles <- c()
  
  cat(paste("\n=== PROCESSING", length(dem_files), "DEM FILES ===\n"))
  
  for (i in seq_along(dem_files)) {
    dem_file <- dem_files[i]
    tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
    
    cat(paste("\n=== PROCESSING", i, "OF", length(dem_files), ":", tile_id, "===\n"))
    
    # Find required matching files
    matches <- find_required_files_for_dem(dem_file, base_dir)
    
    # Check strict requirements
    if (!check_strict_requirements(matches)) {
      failed_tiles <- c(failed_tiles, tile_id)
      next
    }
    
    # Process the tile
    result <- process_single_dem_strict(dem_file, matches, coastline_file, fcf_file, rf_model, qrf_model)
    
    if (!is.null(result)) {
      # Save results
      tile_output_dir <- file.path(colobanthus_output_dir, paste0("tile_", result$tile_id))
      if (!dir.exists(tile_output_dir)) dir.create(tile_output_dir, recursive = TRUE)
      
      writeRaster(result$presence_prob, file.path(tile_output_dir, paste0("presence_prob_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$presence_binary, file.path(tile_output_dir, paste0("presence_binary_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$abundance, file.path(tile_output_dir, paste0("abundance_", result$tile_id, ".tif")), overwrite = TRUE)
      
      successful_tiles <- c(successful_tiles, result$tile_id)
      cat(paste("✓ SUCCESS:", result$tile_id, "- Present pixels:", result$stats$pixels_present, "\n"))
    } else {
      failed_tiles <- c(failed_tiles, tile_id)
    }
    
    gc(full = TRUE)
  }
  
  # Summary
  cat(paste("\n=== FINAL SUMMARY ===\n"))
  cat(paste("Total files:", length(dem_files), "\n"))
  cat(paste("Successful:", length(successful_tiles), "\n"))
  cat(paste("Failed:", length(failed_tiles), "\n"))
  cat(paste("Success rate:", round(length(successful_tiles)/length(dem_files) * 100, 1), "%\n"))
  
  return(list(successful = successful_tiles, failed = failed_tiles))
}

# Execute
tryCatch({
  result <- main_processing()
  cat("✓ Processing completed successfully\n")
}, error = function(e) {
  cat(paste("✗ Critical error:", e$message, "\n"))
  quit(status = 1)
})

cat("=== SCRIPT EXECUTION COMPLETE ===\n")



############# Usnea prediction code ###################################


library(terra)
library(sf)
library(dplyr)
library(randomForest)
library(quantregForest)

cat("=== USNEA STRICT PROCESSING  ===\n")

#### CONFIGURATION ####
# Updated predictor variables as specified
predictors_presence <- c("ALTITUDE", "WEXP", "Easting", "Northing", "HURS_09", "VDEP", "DISTCOAST", "WEFF", "ASP", "Bio18_81")
predictors_abundance <- c("Northing", "ALTITUDE", "WEXP", "ASP", "DISTCOAST", "Easting", "WEFF", "MXCURV", "VDEP")

base_dir <- "/rds/general/user/mrb24/ephemeral"
output_dir <- "/rds/general/user/mrb24/ephemeral"

# Create output directory
usnea_output_dir <- file.path(output_dir, "USNEA_Strict_Processing")
if (!dir.exists(usnea_output_dir)) {
  dir.create(usnea_output_dir, recursive = TRUE)
  cat("Created output directory:", usnea_output_dir, "\n")
}

terraOptions(memfrac = 0.3, tempdir = tempdir())
options(warn = 1)

#### FUNCTIONS ####

# Find all DEM files
find_all_dem_files <- function(base_path) {
  # Check multiple possible DEM directories
  dem_dirs <- c("DEMICEFREE", "DEM", "USNEA_DEM", "ALTITUDE")
  
  for (dem_dir in dem_dirs) {
    dem_path <- file.path(base_path, dem_dir)
    if (dir.exists(dem_path)) {
      all_dem_files <- list.files(dem_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
      if (length(all_dem_files) > 0) {
        cat(paste("Found", length(all_dem_files), "DEM files in", dem_dir, "\n"))
        return(sort(all_dem_files))
      }
    }
  }
  
  stop("No DEM files found in any expected directory")
}

# Spatial extent matching function
find_matching_file_by_extent <- function(reference_file, search_dirs, patterns, min_overlap = 0.7) {
  tryCatch({
    ref_rast <- rast(reference_file)
    ref_ext <- ext(ref_rast)
    
    for (i in seq_along(search_dirs)) {
      search_dir <- search_dirs[i]
      pattern <- patterns[i]
      
      if (!dir.exists(search_dir)) {
        next
      }
      
      # Find all files matching pattern
      all_files <- list.files(search_dir, pattern = paste0(pattern, ".*\\.tif$"), 
                              full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      
      # Test each file for spatial overlap
      for (file in all_files) {
        tryCatch({
          test_rast <- rast(file)
          test_ext <- ext(test_rast)
          
          # Calculate overlap
          ref_vec <- as.vector(ref_ext)
          test_vec <- as.vector(test_ext)
          
          x_overlap <- max(0, min(ref_vec[2], test_vec[2]) - max(ref_vec[1], test_vec[1]))
          y_overlap <- max(0, min(ref_vec[4], test_vec[4]) - max(ref_vec[3], test_vec[3]))
          
          if (x_overlap > 0 && y_overlap > 0) {
            overlap_area <- x_overlap * y_overlap
            ref_area <- (ref_vec[2] - ref_vec[1]) * (ref_vec[4] - ref_vec[3])
            overlap_pct <- overlap_area / ref_area
            
            if (overlap_pct >= min_overlap) {
              return(file)
            }
          }
        }, error = function(e) {
          # Continue to next file
        })
      }
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Find all required matching files - STRICT MODE
find_required_files_for_dem <- function(dem_file, base_path) {
  cat(paste("Finding required files for:", basename(dem_file), "\n"))
  
  # Define search configurations for all required variables
  search_configs <- list(
    wexp = list(dirs = c(file.path(base_path, "WEXP")), patterns = c("WEXP", "WEXPICEFREE")),
    vdep = list(dirs = c(file.path(base_path, "VDEP")), patterns = c("VDEP", "VDEPICEFREE")),
    asp = list(dirs = c(file.path(base_path, "ASP")), patterns = c("ASP", "ASPICEFREE")),
    weff = list(dirs = c(file.path(base_path, "WEFF")), patterns = c("WEFF", "WEFFICEFREE")),
    mxcurv = list(dirs = c(base_path), patterns = c("MXCURV", "MXCURVICEFREE")),
    hurs_09 = list(dirs = c(file.path(base_path, "CLIMATE"), file.path(base_path, "HURS")), patterns = c("HURS", "hurs")),
    bio18_81 = list(dirs = c(file.path(base_path, "CLIMATE"), file.path(base_path, "Bio18_81")), patterns = c("Bio18", "bio18"))
  )
  
  matches <- list()
  
  for (var_name in names(search_configs)) {
    config <- search_configs[[var_name]]
    
    cat(paste("  Searching for", toupper(var_name), "..."))
    
    match_file <- find_matching_file_by_extent(
      dem_file, 
      config$dirs, 
      config$patterns, 
      min_overlap = 0.7
    )
    
    if (!is.null(match_file)) {
      matches[[var_name]] <- match_file
      cat(paste(" ✓", basename(match_file), "\n"))
    } else {
      cat(" ✗ NOT FOUND\n")
    }
  }
  
  return(matches)
}

# Check strict requirements for USNEA - NO DUMMY DATA
check_strict_requirements_usnea <- function(matches) {
  # Core required variables for USNEA presence (excluding coordinates and DISTCOAST)
  required_presence_vars <- c("wexp", "vdep", "weff", "asp", "hurs_09")
  # Core required variables for abundance (excluding coordinates and DISTCOAST) 
  required_abundance_vars <- c("wexp", "asp", "weff", "mxcurv", "vdep")
  
  # Combined requirements
  all_required <- unique(c(required_presence_vars, required_abundance_vars))
  
  missing <- all_required[!all_required %in% names(matches)]
  
  if (length(missing) > 0) {
    cat(paste("  ✗ Missing required variables:", paste(toupper(missing), collapse = ", "), "\n"))
    return(FALSE)
  }
  
  cat("  ✓ All required variables found\n")
  return(TRUE)
}

# Proper distance to coast calculation
calculate_distance_to_coast_proper <- function(raster_template, coastline_path) {
  cat("  Calculating proper distance to coast...\n")
  
  if (!file.exists(coastline_path)) {
    cat("  ✗ Coastline file not found, cannot calculate proper distances\n")
    return(NULL)
  }
  
  tryCatch({
    # Read coastline
    coastline <- st_read(coastline_path, quiet = TRUE)
    coastline <- st_transform(coastline, crs(raster_template))
    
    # Get valid pixel coordinates
    coords <- xyFromCell(raster_template, which(!is.na(values(raster_template))))
    
    if (nrow(coords) == 0) {
      cat("  ✗ No valid pixels found\n")
      return(NULL)
    }
    
    # Convert to sf points
    points <- st_as_sf(data.frame(coords), coords = c("x", "y"), crs = crs(raster_template))
    
    # Calculate minimum distance to coastline for each point
    cat(paste("  Computing distances for", nrow(points), "points...\n"))
    distances <- st_distance(points, coastline)
    min_distances <- apply(distances, 1, min)
    
    # Create distance raster
    dist_raster <- raster_template
    values(dist_raster) <- NA
    valid_cells <- which(!is.na(values(raster_template)))
    values(dist_raster)[valid_cells] <- as.numeric(min_distances)
    
    names(dist_raster) <- "DISTCOAST"
    cat("  ✓ Distance to coast calculated\n")
    return(dist_raster)
    
  }, error = function(e) {
    cat(paste("  ✗ Error calculating distance to coast:", e$message, "\n"))
    return(NULL)
  })
}

# Calculate coordinates
calculate_coordinates <- function(raster_template) {
  coords <- xyFromCell(raster_template, 1:ncell(raster_template))
  easting_raster <- raster_template
  northing_raster <- raster_template
  values(easting_raster) <- coords[, 1]
  values(northing_raster) <- coords[, 2]
  names(easting_raster) <- "Easting"
  names(northing_raster) <- "Northing"
  return(c(easting_raster, northing_raster))
}

# Convert aspect to categorical
convert_aspect_to_categorical <- function(aspect_raster) {
  cat("  Converting aspect to categorical\n")
  
  aspect_vals <- values(aspect_raster, mat = FALSE)
  aspect_cats <- rep(NA_character_, length(aspect_vals))
  
  valid_vals <- !is.na(aspect_vals)
  aspect_cats[valid_vals & (aspect_vals >= 337.5 | aspect_vals < 22.5)] <- "N"
  aspect_cats[valid_vals & aspect_vals >= 22.5 & aspect_vals < 67.5] <- "NE"
  aspect_cats[valid_vals & aspect_vals >= 67.5 & aspect_vals < 112.5] <- "E"
  aspect_cats[valid_vals & aspect_vals >= 112.5 & aspect_vals < 157.5] <- "SE"
  aspect_cats[valid_vals & aspect_vals >= 157.5 & aspect_vals < 202.5] <- "S"
  aspect_cats[valid_vals & aspect_vals >= 202.5 & aspect_vals < 247.5] <- "SW"
  aspect_cats[valid_vals & aspect_vals >= 247.5 & aspect_vals < 292.5] <- "W"
  aspect_cats[valid_vals & aspect_vals >= 292.5 & aspect_vals < 337.5] <- "NW"
  
  aspect_categorical <- aspect_raster
  values(aspect_categorical) <- aspect_cats
  
  return(aspect_categorical)
}

# Strict data cleaning - no tolerance for missing data
clean_prediction_data_strict <- function(prediction_data, predictors) {
  cat("  Strict data cleaning (no NA tolerance)...\n")
  
  # Remove rows with any NA in predictor variables
  complete_rows <- complete.cases(prediction_data[, predictors, drop = FALSE])
  n_removed <- sum(!complete_rows)
  
  if (n_removed > 0) {
    cat(paste("  Removing", n_removed, "rows with NA values\n"))
    prediction_data <- prediction_data[complete_rows, , drop = FALSE]
  }
  
  if (nrow(prediction_data) == 0) {
    cat("  ✗ No complete cases remaining\n")
    return(NULL)
  }
  
  # Handle factor variables
  if ("ASP" %in% predictors && "ASP" %in% names(prediction_data)) {
    prediction_data$ASP <- factor(as.character(prediction_data$ASP), 
                                  levels = c("FLAT", "N", "NE", "E", "SE", "S", "SW", "W", "NW"))
  }
  
  cat(paste("  Clean data:", nrow(prediction_data), "rows\n"))
  return(prediction_data)
}

# Process single DEM file - STRICT MODE
process_usnea_tile_strict <- function(dem_file, matches, coastline_path, rf_model, qrf_model) {
  tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
  cat(paste("Processing:", tile_id, "\n"))
  
  tryCatch({
    # Load DEM as template
    altitude_raster <- rast(dem_file)
    names(altitude_raster) <- "ALTITUDE"
    
    if (sum(!is.na(values(altitude_raster))) == 0) {
      cat("  ✗ No valid altitude data\n")
      return(NULL)
    }
    
    template_raster <- altitude_raster
    all_rasters <- list(ALTITUDE = altitude_raster)
    
    # Load all required variables - NO DUMMY DATA
    variable_map <- list(
      WEXP = "wexp", VDEP = "vdep", ASP = "asp", WEFF = "weff", 
      MXCURV = "mxcurv", HURS_09 = "hurs_09", Bio18_81 = "bio18_81"
    )
    
    for (var_name in names(variable_map)) {
      match_key <- variable_map[[var_name]]
      
      if (match_key %in% names(matches)) {
        cat(paste("  Loading", var_name, "from", basename(matches[[match_key]]), "\n"))
        
        var_raster <- rast(matches[[match_key]])
        names(var_raster) <- var_name
        
        if (!compareGeom(var_raster, template_raster, stopOnError = FALSE)) {
          var_raster <- resample(var_raster, template_raster, method = "bilinear")
        }
        
        # Special handling for ASP (convert to categorical)
        if (var_name == "ASP") {
          var_raster <- convert_aspect_to_categorical(var_raster)
        }
        
        all_rasters[[var_name]] <- var_raster
      } else {
        # Check if this variable is required for the predictors we're using
        required_for_presence <- var_name %in% c("WEXP", "VDEP", "WEFF", "ASP", "HURS_09")
        required_for_abundance <- var_name %in% c("WEXP", "ASP", "WEFF", "MXCURV", "VDEP")
        
        if (required_for_presence || required_for_abundance) {
          cat(paste("  ✗ Missing required variable:", var_name, "\n"))
          return(NULL)
        }
      }
    }
    
    # Calculate coordinates
    cat("  Calculating coordinates...\n")
    coord_rasters <- calculate_coordinates(template_raster)
    all_rasters$Easting <- coord_rasters[[1]]
    all_rasters$Northing <- coord_rasters[[2]]
    
    # Calculate proper distance to coast
    cat("  Calculating distance to coast...\n")
    dist_raster <- calculate_distance_to_coast_proper(template_raster, coastline_path)
    
    if (is.null(dist_raster)) {
      cat("  ✗ Could not calculate distance to coast\n")
      return(NULL)
    }
    
    all_rasters$DISTCOAST <- dist_raster
    
    # Convert to data frame
    cat("  Converting to data frame...\n")
    combined_rasters <- rast(all_rasters)
    prediction_data <- as.data.frame(combined_rasters, xy = TRUE, na.rm = TRUE)
    
    if (nrow(prediction_data) == 0) {
      cat("  ✗ No valid data after conversion\n")
      return(NULL)
    }
    
    # Check which predictors are available
    available_presence_vars <- intersect(predictors_presence, names(prediction_data))
    available_abundance_vars <- intersect(predictors_abundance, names(prediction_data))
    
    cat(paste("  Available presence predictors:", length(available_presence_vars), "of", length(predictors_presence), "\n"))
    cat(paste("  Available abundance predictors:", length(available_abundance_vars), "of", length(predictors_abundance), "\n"))
    
    # Need minimum number of predictors
    if (length(available_presence_vars) < 5) {
      cat("  ✗ Insufficient presence predictors available\n")
      return(NULL)
    }
    
    # Clean data strictly
    presence_data <- clean_prediction_data_strict(prediction_data, available_presence_vars)
    if (is.null(presence_data)) return(NULL)
    
    # Make predictions
    cat(paste("  Making predictions on", nrow(presence_data), "pixels\n"))
    
    prob_pred <- predict(rf_model, newdata = presence_data, type = "prob")[, "1"]
    binary_pred <- ifelse(prob_pred >= 0.5, 1, 0)
    
    # Abundance predictions for present pixels
    abundance_pred <- rep(NA, nrow(presence_data))
    present_idx <- which(binary_pred == 1)
    
    if (length(present_idx) > 0 && length(available_abundance_vars) >= 5) {
      cat(paste("  Calculating abundance for", length(present_idx), "present pixels\n"))
      abundance_data <- clean_prediction_data_strict(presence_data[present_idx, ], available_abundance_vars)
      
      if (!is.null(abundance_data) && nrow(abundance_data) > 0) {
        abundance_predictions <- predict(qrf_model, newdata = abundance_data, what = 0.5)
        abundance_pred[present_idx[1:length(abundance_predictions)]] <- pmax(pmin(abundance_predictions, 100), 0)
      }
    }
    
    # Create output rasters
    prob_raster <- binary_raster <- abundance_raster <- template_raster
    values(prob_raster) <- values(binary_raster) <- values(abundance_raster) <- NA
    
    # Map predictions back to raster
    cell_nums <- cellFromXY(template_raster, presence_data[, c("x", "y")])
    valid_cells <- !is.na(cell_nums)
    
    if (sum(valid_cells) > 0) {
      values(prob_raster)[cell_nums[valid_cells]] <- prob_pred[valid_cells]
      values(binary_raster)[cell_nums[valid_cells]] <- binary_pred[valid_cells]
      values(abundance_raster)[cell_nums[valid_cells]] <- abundance_pred[valid_cells]
    }
    
    names(prob_raster) <- "Presence_Probability"
    names(binary_raster) <- "Presence_Binary"
    names(abundance_raster) <- "Abundance_Prediction"
    
    return(list(
      tile_id = tile_id,
      presence_prob = prob_raster,
      presence_binary = binary_raster,
      abundance = abundance_raster,
      stats = list(
        pixels_present = sum(binary_pred == 1, na.rm = TRUE),
        total_pixels = nrow(presence_data),
        mean_probability = round(mean(prob_pred, na.rm = TRUE), 3),
        available_presence_vars = length(available_presence_vars),
        available_abundance_vars = length(available_abundance_vars)
      )
    ))
    
  }, error = function(e) {
    cat(paste("✗ ERROR:", e$message, "\n"))
    return(NULL)
  })
}

#### LOAD MODELS ####
train_models <- function() {
  training_file <- file.path(base_dir, "hurdle_training_data_list_USNEA_top_predictors.rds")
  if (!file.exists(training_file)) stop("Training data file not found!")
  
  training_data_list <- readRDS(training_file)
  train_data <- training_data_list[[1]]
  
  cat("Training USNEA models with updated predictors...\n")
  
  # Use available predictors from training data
  available_presence <- intersect(predictors_presence, names(train_data))
  available_abundance <- intersect(predictors_abundance, names(train_data))
  
  cat(paste("Training presence model with:", length(available_presence), "variables\n"))
  cat(paste("Training abundance model with:", length(available_abundance), "variables\n"))
  
  rf_clf <- randomForest(
    x = train_data[, available_presence, drop = FALSE],
    y = as.factor(train_data$Presence),
    ntree = 300,
    importance = TRUE
  )
  
  present_data <- train_data[train_data$Presence == "1" | train_data$Presence == 1, ]
  qrf <- quantregForest(
    x = present_data[, available_abundance, drop = FALSE],
    y = as.numeric(present_data$Usnea),
    ntree = 300
  )
  
  return(list(rf_classifier = rf_clf, qrf_model = qrf))
}

#### MAIN PROCESSING ####
main_processing <- function() {
  cat("Loading models...\n")
  models <- train_models()
  rf_model <- models$rf_classifier
  qrf_model <- models$qrf_model
  
  # Find all DEM files
  dem_files <- find_all_dem_files(base_dir)
  
  # Coastline file
  coastline_file <- file.path(base_dir, "COAST", "add_coastline_high_res_line_v7_10.shp")
  if (!file.exists(coastline_file)) {
    stop("Coastline file not found - required for proper distance calculation")
  }
  
  successful_tiles <- c()
  failed_tiles <- c()
  
  cat(paste("\n=== PROCESSING", length(dem_files), "DEM FILES ===\n"))
  
  for (i in seq_along(dem_files)) {
    dem_file <- dem_files[i]
    tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
    
    cat(paste("\n=== PROCESSING", i, "OF", length(dem_files), ":", tile_id, "===\n"))
    
    # Find required matching files
    matches <- find_required_files_for_dem(dem_file, base_dir)
    
    # Check strict requirements
    if (!check_strict_requirements_usnea(matches)) {
      failed_tiles <- c(failed_tiles, tile_id)
      next
    }
    
    # Process the tile
    result <- process_usnea_tile_strict(dem_file, matches, coastline_file, rf_model, qrf_model)
    
    if (!is.null(result)) {
      # Save results
      tile_output_dir <- file.path(usnea_output_dir, paste0("tile_", result$tile_id))
      if (!dir.exists(tile_output_dir)) dir.create(tile_output_dir, recursive = TRUE)
      
      writeRaster(result$presence_prob, file.path(tile_output_dir, paste0("presence_prob_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$presence_binary, file.path(tile_output_dir, paste0("presence_binary_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$abundance, file.path(tile_output_dir, paste0("abundance_", result$tile_id, ".tif")), overwrite = TRUE)
      
      successful_tiles <- c(successful_tiles, result$tile_id)
      cat(paste("✓ SUCCESS:", result$tile_id, "- Present pixels:", result$stats$pixels_present, "\n"))
      cat(paste("  Presence vars:", result$stats$available_presence_vars, "Abundance vars:", result$stats$available_abundance_vars, "\n"))
    } else {
      failed_tiles <- c(failed_tiles, tile_id)
    }
    
    gc(full = TRUE)
  }
  
  # Summary
  cat(paste("\n=== FINAL SUMMARY ===\n"))
  cat(paste("Total files:", length(dem_files), "\n"))
  cat(paste("Successful:", length(successful_tiles), "\n"))
  cat(paste("Failed:", length(failed_tiles), "\n"))
  cat(paste("Success rate:", round(length(successful_tiles)/length(dem_files) * 100, 1), "%\n"))
  cat(paste("Output directory:", usnea_output_dir, "\n"))
  
  if (length(failed_tiles) > 0) {
    cat("\n✗ Failed tiles (missing required data):\n")
    for (tile in head(failed_tiles, 10)) {
      cat(paste("  -", tile, "\n"))
    }
    if (length(failed_tiles) > 10) {
      cat(paste("  ... and", length(failed_tiles) - 10, "more\n"))
    }
  }
  
  return(list(successful = successful_tiles, failed = failed_tiles))
}

# Execute main processing
tryCatch({
  result <- main_processing()
  cat("✓ USNEA processing completed successfully\n")
  cat("✓ All predictions use proper distance-to-coast calculations\n")
  cat("✓ No dummy data was used - only real environmental variables\n")
}, error = function(e) {
  cat(paste("✗ Critical error:", e$message, "\n"))
  traceback()
  quit(status = 1)
})

cat("=== USNEA STRICT PROCESSING COMPLETE ===\n")





########### Sanionia current predictions code #####################################


library(terra)
library(sf)
library(dplyr)
library(randomForest)
library(quantregForest)

cat("=== SANIONIA STRICT PROCESSING  ===\n")

#### CONFIGURATION ####
# Updated predictor variables as specified
predictors_presence <- c("ALTITUDE", "Easting", "WEFF", "WEXP", "DISTCOAST", "Northing", "ASP", "VDEP", "Bio15_81", "FLAC", "TRI")
predictors_abundance <- c("VDEP", "ASP", "DISTCOAST", "SVF", "WEFF", "Easting", "FLAC", "TRI", "ALTITUDE", "Northing")

base_dir <- "/rds/general/user/mrb24/ephemeral"
output_dir <- "/rds/general/user/mrb24/ephemeral"

# Create output directory
sanionia_output_dir <- file.path(output_dir, "SANIONIA_Strict_Processing")
if (!dir.exists(sanionia_output_dir)) {
  dir.create(sanionia_output_dir, recursive = TRUE)
  cat("Created output directory:", sanionia_output_dir, "\n")
}

terraOptions(memfrac = 0.3, tempdir = tempdir())
options(warn = 1)

#### FUNCTIONS ####

# Find all DEM files
find_all_dem_files <- function(base_path) {
  # Check multiple possible DEM directories
  dem_dirs <- c("DEMICEFREE", "DEM", "SANIONIA_DEM", "ALTITUDE")
  
  for (dem_dir in dem_dirs) {
    dem_path <- file.path(base_path, dem_dir)
    if (dir.exists(dem_path)) {
      all_dem_files <- list.files(dem_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
      if (length(all_dem_files) > 0) {
        cat(paste("Found", length(all_dem_files), "DEM files in", dem_dir, "\n"))
        return(sort(all_dem_files))
      }
    }
  }
  
  stop("No DEM files found in any expected directory")
}

# Spatial extent matching function
find_matching_file_by_extent <- function(reference_file, search_dirs, patterns, min_overlap = 0.7) {
  tryCatch({
    ref_rast <- rast(reference_file)
    ref_ext <- ext(ref_rast)
    
    for (i in seq_along(search_dirs)) {
      search_dir <- search_dirs[i]
      pattern <- patterns[i]
      
      if (!dir.exists(search_dir)) {
        next
      }
      
      # Find all files matching pattern
      all_files <- list.files(search_dir, pattern = paste0(pattern, ".*\\.tif$"), 
                              full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      
      # Test each file for spatial overlap
      for (file in all_files) {
        tryCatch({
          test_rast <- rast(file)
          test_ext <- ext(test_rast)
          
          # Calculate overlap
          ref_vec <- as.vector(ref_ext)
          test_vec <- as.vector(test_ext)
          
          x_overlap <- max(0, min(ref_vec[2], test_vec[2]) - max(ref_vec[1], test_vec[1]))
          y_overlap <- max(0, min(ref_vec[4], test_vec[4]) - max(ref_vec[3], test_vec[3]))
          
          if (x_overlap > 0 && y_overlap > 0) {
            overlap_area <- x_overlap * y_overlap
            ref_area <- (ref_vec[2] - ref_vec[1]) * (ref_vec[4] - ref_vec[3])
            overlap_pct <- overlap_area / ref_area
            
            if (overlap_pct >= min_overlap) {
              return(file)
            }
          }
        }, error = function(e) {
          # Continue to next file
        })
      }
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Find all required matching files - STRICT MODE
find_required_files_for_dem <- function(dem_file, base_path) {
  cat(paste("Finding required files for:", basename(dem_file), "\n"))
  
  # Define search configurations for all required variables
  search_configs <- list(
    wexp = list(dirs = c(file.path(base_path, "WEXP")), patterns = c("WEXP", "WEXPICEFREE")),
    weff = list(dirs = c(file.path(base_path, "WEFF")), patterns = c("WEFF", "WEFFICEFREE")),
    vdep = list(dirs = c(file.path(base_path, "VDEP")), patterns = c("VDEP", "VDEPICEFREE")),
    asp = list(dirs = c(base_path), patterns = c("ASP", "ASPICEFREE")),
    flac = list(dirs = c(base_path), patterns = c("FLAC", "FLACICEFREE")),
    tri = list(dirs = c(base_path), patterns = c("TRI", "TRIICEFREE")),
    svf = list(dirs = c(base_path), patterns = c("SVF", "SVFICEFREE")),
    bio15_81 = list(dirs = c(base_path), patterns = c("Sanionia_Bio15_Climate_Envelope_Values", "Bio15", "bio15"))
  )
  
  matches <- list()
  
  for (var_name in names(search_configs)) {
    config <- search_configs[[var_name]]
    
    cat(paste("  Searching for", toupper(var_name), "..."))
    
    match_file <- find_matching_file_by_extent(
      dem_file, 
      config$dirs, 
      config$patterns, 
      min_overlap = 0.7
    )
    
    if (!is.null(match_file)) {
      matches[[var_name]] <- match_file
      cat(paste(" ✓", basename(match_file), "\n"))
    } else {
      cat(" ✗ NOT FOUND\n")
    }
  }
  
  return(matches)
}

# Check strict requirements for SANIONIA - NO DUMMY DATA
check_strict_requirements_sanionia <- function(matches) {
  # Core required variables for SANIONIA presence (excluding coordinates and DISTCOAST)
  required_presence_vars <- c("wexp", "weff", "vdep", "asp", "flac", "tri", "bio15_81")
  # Core required variables for abundance (excluding coordinates and DISTCOAST)
  required_abundance_vars <- c("vdep", "asp", "svf", "weff", "flac", "tri")
  
  # Combined requirements
  all_required <- unique(c(required_presence_vars, required_abundance_vars))
  
  missing <- all_required[!all_required %in% names(matches)]
  
  if (length(missing) > 0) {
    cat(paste("  ✗ Missing required variables:", paste(toupper(missing), collapse = ", "), "\n"))
    return(FALSE)
  }
  
  cat("  ✓ All required variables found\n")
  return(TRUE)
}

# Proper distance to coast calculation
calculate_distance_to_coast_proper <- function(raster_template, coastline_path) {
  cat("  Calculating proper distance to coast...\n")
  
  if (!file.exists(coastline_path)) {
    cat("  ✗ Coastline file not found, cannot calculate proper distances\n")
    return(NULL)
  }
  
  tryCatch({
    # Read coastline
    coastline <- st_read(coastline_path, quiet = TRUE)
    coastline <- st_transform(coastline, crs(raster_template))
    
    # Get valid pixel coordinates
    coords <- xyFromCell(raster_template, which(!is.na(values(raster_template))))
    
    if (nrow(coords) == 0) {
      cat("  ✗ No valid pixels found\n")
      return(NULL)
    }
    
    # Convert to sf points
    points <- st_as_sf(data.frame(coords), coords = c("x", "y"), crs = crs(raster_template))
    
    # Calculate minimum distance to coastline for each point
    cat(paste("  Computing distances for", nrow(points), "points...\n"))
    distances <- st_distance(points, coastline)
    min_distances <- apply(distances, 1, min)
    
    # Create distance raster
    dist_raster <- raster_template
    values(dist_raster) <- NA
    valid_cells <- which(!is.na(values(raster_template)))
    values(dist_raster)[valid_cells] <- as.numeric(min_distances)
    
    names(dist_raster) <- "DISTCOAST"
    cat("  ✓ Distance to coast calculated\n")
    return(dist_raster)
    
  }, error = function(e) {
    cat(paste("  ✗ Error calculating distance to coast:", e$message, "\n"))
    return(NULL)
  })
}

# Calculate coordinates
calculate_coordinates <- function(raster_template) {
  coords <- xyFromCell(raster_template, 1:ncell(raster_template))
  easting_raster <- raster_template
  northing_raster <- raster_template
  values(easting_raster) <- coords[, 1]
  values(northing_raster) <- coords[, 2]
  names(easting_raster) <- "Easting"
  names(northing_raster) <- "Northing"
  return(c(easting_raster, northing_raster))
}

# Convert aspect to categorical
convert_aspect_to_categorical <- function(aspect_raster) {
  cat("  Converting aspect to categorical\n")
  
  aspect_vals <- values(aspect_raster, mat = FALSE)
  aspect_cats <- rep(NA_character_, length(aspect_vals))
  
  valid_vals <- !is.na(aspect_vals)
  aspect_cats[valid_vals & (aspect_vals >= 337.5 | aspect_vals < 22.5)] <- "N"
  aspect_cats[valid_vals & aspect_vals >= 22.5 & aspect_vals < 67.5] <- "NE"
  aspect_cats[valid_vals & aspect_vals >= 67.5 & aspect_vals < 112.5] <- "E"
  aspect_cats[valid_vals & aspect_vals >= 112.5 & aspect_vals < 157.5] <- "SE"
  aspect_cats[valid_vals & aspect_vals >= 157.5 & aspect_vals < 202.5] <- "S"
  aspect_cats[valid_vals & aspect_vals >= 202.5 & aspect_vals < 247.5] <- "SW"
  aspect_cats[valid_vals & aspect_vals >= 247.5 & aspect_vals < 292.5] <- "W"
  aspect_cats[valid_vals & aspect_vals >= 292.5 & aspect_vals < 337.5] <- "NW"
  
  aspect_categorical <- aspect_raster
  values(aspect_categorical) <- aspect_cats
  
  return(aspect_categorical)
}

# Strict data cleaning - no tolerance for missing data
clean_prediction_data_strict <- function(prediction_data, predictors) {
  cat("  Strict data cleaning (no NA tolerance)...\n")
  
  # Remove rows with any NA in predictor variables
  complete_rows <- complete.cases(prediction_data[, predictors, drop = FALSE])
  n_removed <- sum(!complete_rows)
  
  if (n_removed > 0) {
    cat(paste("  Removing", n_removed, "rows with NA values\n"))
    prediction_data <- prediction_data[complete_rows, , drop = FALSE]
  }
  
  if (nrow(prediction_data) == 0) {
    cat("  ✗ No complete cases remaining\n")
    return(NULL)
  }
  
  # Handle factor variables
  if ("ASP" %in% predictors && "ASP" %in% names(prediction_data)) {
    prediction_data$ASP <- factor(as.character(prediction_data$ASP), 
                                  levels = c("FLAT", "N", "NE", "E", "SE", "S", "SW", "W", "NW"))
  }
  
  cat(paste("  Clean data:", nrow(prediction_data), "rows\n"))
  return(prediction_data)
}

# Process single DEM file - STRICT MODE
process_sanionia_tile_strict <- function(dem_file, matches, coastline_path, rf_model, qrf_model) {
  tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
  cat(paste("Processing:", tile_id, "\n"))
  
  tryCatch({
    # Load DEM as template
    altitude_raster <- rast(dem_file)
    names(altitude_raster) <- "ALTITUDE"
    
    if (sum(!is.na(values(altitude_raster))) == 0) {
      cat("  ✗ No valid altitude data\n")
      return(NULL)
    }
    
    template_raster <- altitude_raster
    all_rasters <- list(ALTITUDE = altitude_raster)
    
    # Load Bio15 climate data first for masking
    cat("  Loading Bio15 climate envelope...\n")
    bio15_raster <- rast(matches$bio15_81)
    names(bio15_raster) <- "Bio15_81"
    
    if (!compareGeom(bio15_raster, template_raster, stopOnError = FALSE)) {
      bio15_raster <- resample(bio15_raster, template_raster, method = "bilinear")
    }
    all_rasters$Bio15_81 <- bio15_raster
    
    # Create climate mask - only areas with valid Bio15 data
    climate_mask <- !is.na(bio15_raster) & !is.na(altitude_raster)
    cat(paste("  Climate envelope covers", sum(values(climate_mask), na.rm = TRUE), "pixels\n"))
    
    # Check if there's any overlap
    if (sum(values(climate_mask), na.rm = TRUE) == 0) {
      cat("  No overlap with climate envelope - skipping\n")
      return(NULL)
    }
    
    # Apply climate mask to template
    template_raster <- mask(template_raster, climate_mask)
    all_rasters$ALTITUDE <- template_raster
    
    # Load all required variables - NO DUMMY DATA
    variable_map <- list(
      WEXP = "wexp", WEFF = "weff", VDEP = "vdep", ASP = "asp", 
      FLAC = "flac", TRI = "tri", SVF = "svf"
    )
    
    for (var_name in names(variable_map)) {
      match_key <- variable_map[[var_name]]
      
      if (match_key %in% names(matches)) {
        cat(paste("  Loading", var_name, "from", basename(matches[[match_key]]), "\n"))
        
        var_raster <- rast(matches[[match_key]])
        names(var_raster) <- var_name
        
        if (!compareGeom(var_raster, template_raster, stopOnError = FALSE)) {
          var_raster <- resample(var_raster, template_raster, method = "bilinear")
        }
        
        # Apply climate mask
        var_raster <- mask(var_raster, climate_mask)
        
        # Special handling for ASP (convert to categorical)
        if (var_name == "ASP") {
          var_raster <- convert_aspect_to_categorical(var_raster)
        }
        
        all_rasters[[var_name]] <- var_raster
      } else {
        # Check if this variable is required for the predictors we're using
        required_for_presence <- var_name %in% c("WEXP", "WEFF", "VDEP", "ASP", "FLAC", "TRI")
        required_for_abundance <- var_name %in% c("VDEP", "ASP", "SVF", "WEFF", "FLAC", "TRI")
        
        if (required_for_presence || required_for_abundance) {
          cat(paste("  ✗ Missing required variable:", var_name, "\n"))
          return(NULL)
        }
      }
    }
    
    # Calculate coordinates
    cat("  Calculating coordinates...\n")
    coord_rasters <- calculate_coordinates(template_raster)
    all_rasters$Easting <- mask(coord_rasters[[1]], climate_mask)
    all_rasters$Northing <- mask(coord_rasters[[2]], climate_mask)
    
    # Calculate proper distance to coast
    cat("  Calculating distance to coast...\n")
    dist_raster <- calculate_distance_to_coast_proper(template_raster, coastline_path)
    
    if (is.null(dist_raster)) {
      cat("  ✗ Could not calculate distance to coast\n")
      return(NULL)
    }
    
    all_rasters$DISTCOAST <- mask(dist_raster, climate_mask)
    
    # Convert to data frame
    cat("  Converting to data frame...\n")
    combined_rasters <- rast(all_rasters)
    prediction_data <- as.data.frame(combined_rasters, xy = TRUE, na.rm = TRUE)
    
    if (nrow(prediction_data) == 0) {
      cat("  ✗ No valid data after conversion\n")
      return(NULL)
    }
    
    # Check which predictors are available
    available_presence_vars <- intersect(predictors_presence, names(prediction_data))
    available_abundance_vars <- intersect(predictors_abundance, names(prediction_data))
    
    cat(paste("  Available presence predictors:", length(available_presence_vars), "of", length(predictors_presence), "\n"))
    cat(paste("  Available abundance predictors:", length(available_abundance_vars), "of", length(predictors_abundance), "\n"))
    
    # Need minimum number of predictors
    if (length(available_presence_vars) < 7) {
      cat("  ✗ Insufficient presence predictors available\n")
      return(NULL)
    }
    
    # Clean data strictly
    presence_data <- clean_prediction_data_strict(prediction_data, available_presence_vars)
    if (is.null(presence_data)) return(NULL)
    
    # Make predictions
    cat(paste("  Making predictions on", nrow(presence_data), "pixels\n"))
    
    prob_pred <- predict(rf_model, newdata = presence_data, type = "prob")[, "1"]
    binary_pred <- ifelse(prob_pred >= 0.5, 1, 0)
    
    # Abundance predictions for present pixels
    abundance_pred <- rep(NA, nrow(presence_data))
    present_idx <- which(binary_pred == 1)
    
    if (length(present_idx) > 0 && length(available_abundance_vars) >= 6) {
      cat(paste("  Calculating abundance for", length(present_idx), "present pixels\n"))
      abundance_data <- clean_prediction_data_strict(presence_data[present_idx, ], available_abundance_vars)
      
      if (!is.null(abundance_data) && nrow(abundance_data) > 0) {
        abundance_predictions <- predict(qrf_model, newdata = abundance_data, what = 0.5)
        abundance_pred[present_idx[1:length(abundance_predictions)]] <- pmax(pmin(abundance_predictions, 100), 0)
      }
    }
    
    # Create output rasters (using original template before masking for full extent)
    original_template <- rast(dem_file)
    prob_raster <- binary_raster <- abundance_raster <- original_template
    values(prob_raster) <- values(binary_raster) <- values(abundance_raster) <- NA
    
    # Map predictions back to raster
    cell_nums <- cellFromXY(original_template, presence_data[, c("x", "y")])
    valid_cells <- !is.na(cell_nums)
    
    if (sum(valid_cells) > 0) {
      values(prob_raster)[cell_nums[valid_cells]] <- prob_pred[valid_cells]
      values(binary_raster)[cell_nums[valid_cells]] <- binary_pred[valid_cells]
      values(abundance_raster)[cell_nums[valid_cells]] <- abundance_pred[valid_cells]
    }
    
    names(prob_raster) <- "Presence_Probability"
    names(binary_raster) <- "Presence_Binary"
    names(abundance_raster) <- "Abundance_Prediction"
    
    return(list(
      tile_id = tile_id,
      presence_prob = prob_raster,
      presence_binary = binary_raster,
      abundance = abundance_raster,
      climate_mask = climate_mask,
      stats = list(
        pixels_present = sum(binary_pred == 1, na.rm = TRUE),
        total_pixels = nrow(presence_data),
        mean_probability = round(mean(prob_pred, na.rm = TRUE), 3),
        climate_pixels = sum(values(climate_mask), na.rm = TRUE),
        available_presence_vars = length(available_presence_vars),
        available_abundance_vars = length(available_abundance_vars)
      )
    ))
    
  }, error = function(e) {
    cat(paste("✗ ERROR:", e$message, "\n"))
    return(NULL)
  })
}

#### LOAD MODELS ####
train_models <- function() {
  training_file <- file.path(base_dir, "hurdle_training_data_list_SANIONIA_top_predictors.rds")
  if (!file.exists(training_file)) stop("Training data file not found!")
  
  training_data_list <- readRDS(training_file)
  train_data <- training_data_list[[1]]
  
  cat("Training SANIONIA models with updated predictors...\n")
  
  # Use available predictors from training data
  available_presence <- intersect(predictors_presence, names(train_data))
  available_abundance <- intersect(predictors_abundance, names(train_data))
  
  cat(paste("Training presence model with:", length(available_presence), "variables\n"))
  cat(paste("Training abundance model with:", length(available_abundance), "variables\n"))
  
  rf_clf <- randomForest(
    x = train_data[, available_presence, drop = FALSE],
    y = as.factor(train_data$Presence),
    ntree = 300,
    importance = TRUE
  )
  
  present_data <- train_data[train_data$Presence == "1" | train_data$Presence == 1, ]
  qrf <- quantregForest(
    x = present_data[, available_abundance, drop = FALSE],
    y = as.numeric(present_data$Sanionia),
    ntree = 300
  )
  
  return(list(rf_classifier = rf_clf, qrf_model = qrf))
}

#### MAIN PROCESSING ####
main_processing <- function() {
  cat("Loading models...\n")
  models <- train_models()
  rf_model <- models$rf_classifier
  qrf_model <- models$qrf_model
  
  # Find all DEM files
  dem_files <- find_all_dem_files(base_dir)
  
  # Coastline file
  coastline_file <- file.path(base_dir, "COAST", "add_coastline_high_res_line_v7_10.shp")
  if (!file.exists(coastline_file)) {
    stop("Coastline file not found - required for proper distance calculation")
  }
  
  successful_tiles <- c()
  failed_tiles <- c()
  
  cat(paste("\n=== PROCESSING", length(dem_files), "DEM FILES ===\n"))
  
  for (i in seq_along(dem_files)) {
    dem_file <- dem_files[i]
    tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
    
    cat(paste("\n=== PROCESSING", i, "OF", length(dem_files), ":", tile_id, "===\n"))
    
    # Find required matching files
    matches <- find_required_files_for_dem(dem_file, base_dir)
    
    # Check strict requirements
    if (!check_strict_requirements_sanionia(matches)) {
      failed_tiles <- c(failed_tiles, tile_id)
      next
    }
    
    # Process the tile
    result <- process_sanionia_tile_strict(dem_file, matches, coastline_file, rf_model, qrf_model)
    
    if (!is.null(result)) {
      # Save results
      tile_output_dir <- file.path(sanionia_output_dir, paste0("tile_", result$tile_id))
      if (!dir.exists(tile_output_dir)) dir.create(tile_output_dir, recursive = TRUE)
      
      writeRaster(result$presence_prob, file.path(tile_output_dir, paste0("presence_prob_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$presence_binary, file.path(tile_output_dir, paste0("presence_binary_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$abundance, file.path(tile_output_dir, paste0("abundance_", result$tile_id, ".tif")), overwrite = TRUE)
      
      # Save climate mask for reference
      writeRaster(result$climate_mask, file.path(tile_output_dir, paste0("climate_mask_", result$tile_id, ".tif")), overwrite = TRUE)
      
      successful_tiles <- c(successful_tiles, result$tile_id)
      cat(paste("✓ SUCCESS:", result$tile_id, "- Present pixels:", result$stats$pixels_present, "\n"))
      cat(paste("  Climate pixels:", result$stats$climate_pixels, "Presence vars:", result$stats$available_presence_vars, "\n"))
    } else {
      failed_tiles <- c(failed_tiles, tile_id)
    }
    
    gc(full = TRUE)
  }
  
  # Summary
  cat(paste("\n=== FINAL SUMMARY ===\n"))
  cat(paste("Total files:", length(dem_files), "\n"))
  cat(paste("Successful:", length(successful_tiles), "\n"))
  cat(paste("Failed:", length(failed_tiles), "\n"))
  cat(paste("Success rate:", round(length(successful_tiles)/length(dem_files) * 100, 1), "%\n"))
  cat(paste("Output directory:", sanionia_output_dir, "\n"))
  
  if (length(failed_tiles) > 0) {
    cat("\n✗ Failed tiles (missing required data):\n")
    for (tile in head(failed_tiles, 10)) {
      cat(paste("  -", tile, "\n"))
    }
    if (length(failed_tiles) > 10) {
      cat(paste("  ... and", length(failed_tiles) - 10, "more\n"))
    }
  }
  
  return(list(successful = successful_tiles, failed = failed_tiles))
}

# Execute main processing
tryCatch({
  result <- main_processing()
  cat("✓ SANIONIA processing completed successfully\n")
  cat("✓ All predictions clipped to Bio15 climate envelope\n")
  cat("✓ Proper distance-to-coast calculations applied\n")
  cat("✓ No dummy data used - only real environmental variables\n")
}, error = function(e) {
  cat(paste("✗ Critical error:", e$message, "\n"))
  traceback()
  quit(status = 1)
})

cat("=== SANIONIA STRICT PROCESSING COMPLETE ===\n")




################ Deschampsia current predictions #############################


#!/usr/bin/env Rscript

# DESCHAMPSIA Strict Processing - No Dummy Values, Proper Distance Calculation
# Author: Generated for HPC processing
# Date: 2025

library(terra)
library(sf)
library(dplyr)
library(randomForest)
library(quantregForest)

cat("=== DESCHAMPSIA STRICT PROCESSING - NO DUMMY VALUES ===\n")

#### CONFIGURATION ####
# Updated predictor variables as specified
predictors_presence <- c("DISTCOAST", "WEXP", "VDEP", "EFAH", "Easting", "MPI", "Northing", "WEFF", "LSF", "ASP", "HillSh", "PR_01", "MBI", "PDIS_CLASS", "RSP", "FLAC", "TAS_08")
predictors_abundance <- c("TPW", "LSF", "Northing", "MPI", "Easting", "ASP", "DISTCOAST", "EFAH", "MBI", "WEXP", "MNCURV", "FLAC", "SVF")

base_dir <- "/rds/general/user/mrb24/ephemeral"
output_dir <- "/rds/general/user/mrb24/ephemeral"

# Create output directory
deschampsia_output_dir <- file.path(output_dir, "DESCHAMPSIA_Strict_Processing")
if (!dir.exists(deschampsia_output_dir)) {
  dir.create(deschampsia_output_dir, recursive = TRUE)
  cat("Created output directory:", deschampsia_output_dir, "\n")
}

terraOptions(memfrac = 0.3, tempdir = tempdir())
options(warn = 1)

#### FUNCTIONS ####

# Find all DEM files
find_all_dem_files <- function(base_path) {
  # Check multiple possible DEM directories
  dem_dirs <- c("DEMICEFREE", "DEM", "DESCHAMPSIA_DEM", "ALTITUDE")
  
  for (dem_dir in dem_dirs) {
    dem_path <- file.path(base_path, dem_dir)
    if (dir.exists(dem_path)) {
      all_dem_files <- list.files(dem_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
      if (length(all_dem_files) > 0) {
        cat(paste("Found", length(all_dem_files), "DEM files in", dem_dir, "\n"))
        return(sort(all_dem_files))
      }
    }
  }
  
  stop("No DEM files found in any expected directory")
}

# Spatial extent matching function
find_matching_file_by_extent <- function(reference_file, search_dirs, patterns, min_overlap = 0.7) {
  tryCatch({
    ref_rast <- rast(reference_file)
    ref_ext <- ext(ref_rast)
    
    for (i in seq_along(search_dirs)) {
      search_dir <- search_dirs[i]
      pattern <- patterns[i]
      
      if (!dir.exists(search_dir)) {
        next
      }
      
      # Find all files matching pattern
      all_files <- list.files(search_dir, pattern = paste0(pattern, ".*\\.tif$"), 
                              full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      
      # Test each file for spatial overlap
      for (file in all_files) {
        tryCatch({
          test_rast <- rast(file)
          test_ext <- ext(test_rast)
          
          # Calculate overlap
          ref_vec <- as.vector(ref_ext)
          test_vec <- as.vector(test_ext)
          
          x_overlap <- max(0, min(ref_vec[2], test_vec[2]) - max(ref_vec[1], test_vec[1]))
          y_overlap <- max(0, min(ref_vec[4], test_vec[4]) - max(ref_vec[3], test_vec[3]))
          
          if (x_overlap > 0 && y_overlap > 0) {
            overlap_area <- x_overlap * y_overlap
            ref_area <- (ref_vec[2] - ref_vec[1]) * (ref_vec[4] - ref_vec[3])
            overlap_pct <- overlap_area / ref_area
            
            if (overlap_pct >= min_overlap) {
              return(file)
            }
          }
        }, error = function(e) {
          # Continue to next file
        })
      }
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Find all required matching files - STRICT MODE
find_required_files_for_dem <- function(dem_file, base_path) {
  cat(paste("Finding required files for:", basename(dem_file), "\n"))
  
  # Define search configurations for all required variables
  search_configs <- list(
    wexp = list(dirs = c(file.path(base_path, "WEXP")), patterns = c("WEXP", "WEXPICEFREE")),
    vdep = list(dirs = c(file.path(base_path, "VDEP")), patterns = c("VDEP", "VDEPICEFREE")),
    efah = list(dirs = c(base_path), patterns = c("EFAH", "EFAHICEFREE")),
    mpi = list(dirs = c(base_path), patterns = c("MPI", "MPIICEFREE")),
    weff = list(dirs = c(file.path(base_path, "WEFF")), patterns = c("WEFF", "WEFFICEFREE")),
    lsf = list(dirs = c(base_path), patterns = c("LSF", "LSFICEFREE")),
    asp = list(dirs = c(base_path), patterns = c("ASP", "ASPICEFREE")),
    hillsh = list(dirs = c(base_path), patterns = c("HillSh", "HILLSH", "HillShicefree")),
    mbi = list(dirs = c(base_path), patterns = c("MBI", "MBIICEFREE")),
    pdis_class = list(dirs = c(base_path), patterns = c("PDIS_CLASS", "PDIS_CLASSICEFREE")),
    rsp = list(dirs = c(base_path), patterns = c("RSP", "RSPICEFREE")),
    flac = list(dirs = c(base_path), patterns = c("FLAC", "FLACICEFREE")),
    tpw = list(dirs = c(base_path), patterns = c("TPW", "TPWICEFREE")),
    mncurv = list(dirs = c(base_path), patterns = c("MNCURV", "MNCURVICEFREE")),
    svf = list(dirs = c(base_path), patterns = c("SVF", "SVFICEFREE")),
    tas_08 = list(dirs = c(file.path(base_path, "TAS_08"), base_path), patterns = c("Deschampsia_TAS_08_Climate_Envelope", "TAS_08", "tas_08")),
    pr_01 = list(dirs = c(file.path(base_path, "PR_01"), base_path), patterns = c("Deschampsia_PR_01_Climate_Envelope", "PR_01", "pr_01"))
  )
  
  matches <- list()
  
  for (var_name in names(search_configs)) {
    config <- search_configs[[var_name]]
    
    cat(paste("  Searching for", toupper(var_name), "..."))
    
    match_file <- find_matching_file_by_extent(
      dem_file, 
      config$dirs, 
      config$patterns, 
      min_overlap = 0.7
    )
    
    if (!is.null(match_file)) {
      matches[[var_name]] <- match_file
      cat(paste(" ✓", basename(match_file), "\n"))
    } else {
      cat(" ✗ NOT FOUND\n")
    }
  }
  
  return(matches)
}

# Check strict requirements for DESCHAMPSIA - NO DUMMY DATA
check_strict_requirements_deschampsia <- function(matches) {
  # Core required variables for DESCHAMPSIA presence (excluding coordinates and DISTCOAST)
  required_presence_vars <- c("wexp", "vdep", "efah", "mpi", "weff", "lsf", "asp", "hillsh", "mbi", "pdis_class", "rsp", "flac", "tas_08", "pr_01")
  # Core required variables for abundance (excluding coordinates and DISTCOAST)
  required_abundance_vars <- c("tpw", "lsf", "mpi", "asp", "efah", "mbi", "wexp", "mncurv", "flac", "svf")
  
  # Combined requirements
  all_required <- unique(c(required_presence_vars, required_abundance_vars))
  
  missing <- all_required[!all_required %in% names(matches)]
  
  if (length(missing) > 0) {
    cat(paste("  ✗ Missing required variables:", paste(toupper(missing), collapse = ", "), "\n"))
    return(FALSE)
  }
  
  cat("  ✓ All required variables found\n")
  return(TRUE)
}

# Proper distance to coast calculation
calculate_distance_to_coast_proper <- function(raster_template, coastline_path) {
  cat("  Calculating proper distance to coast...\n")
  
  if (!file.exists(coastline_path)) {
    cat("  ✗ Coastline file not found, cannot calculate proper distances\n")
    return(NULL)
  }
  
  tryCatch({
    # Read coastline
    coastline <- st_read(coastline_path, quiet = TRUE)
    coastline <- st_transform(coastline, crs(raster_template))
    
    # Get valid pixel coordinates
    coords <- xyFromCell(raster_template, which(!is.na(values(raster_template))))
    
    if (nrow(coords) == 0) {
      cat("  ✗ No valid pixels found\n")
      return(NULL)
    }
    
    # Convert to sf points
    points <- st_as_sf(data.frame(coords), coords = c("x", "y"), crs = crs(raster_template))
    
    # Calculate minimum distance to coastline for each point
    cat(paste("  Computing distances for", nrow(points), "points...\n"))
    distances <- st_distance(points, coastline)
    min_distances <- apply(distances, 1, min)
    
    # Create distance raster
    dist_raster <- raster_template
    values(dist_raster) <- NA
    valid_cells <- which(!is.na(values(raster_template)))
    values(dist_raster)[valid_cells] <- as.numeric(min_distances)
    
    names(dist_raster) <- "DISTCOAST"
    cat("  ✓ Distance to coast calculated\n")
    return(dist_raster)
    
  }, error = function(e) {
    cat(paste("  ✗ Error calculating distance to coast:", e$message, "\n"))
    return(NULL)
  })
}

# Calculate coordinates
calculate_coordinates <- function(raster_template) {
  coords <- xyFromCell(raster_template, 1:ncell(raster_template))
  easting_raster <- raster_template
  northing_raster <- raster_template
  values(easting_raster) <- coords[, 1]
  values(northing_raster) <- coords[, 2]
  names(easting_raster) <- "Easting"
  names(northing_raster) <- "Northing"
  return(c(easting_raster, northing_raster))
}

# Convert aspect to categorical
convert_aspect_to_categorical <- function(aspect_raster) {
  cat("  Converting aspect to categorical\n")
  
  aspect_vals <- values(aspect_raster, mat = FALSE)
  aspect_cats <- rep(NA_character_, length(aspect_vals))
  
  valid_vals <- !is.na(aspect_vals)
  aspect_cats[valid_vals & (aspect_vals >= 337.5 | aspect_vals < 22.5)] <- "N"
  aspect_cats[valid_vals & aspect_vals >= 22.5 & aspect_vals < 67.5] <- "NE"
  aspect_cats[valid_vals & aspect_vals >= 67.5 & aspect_vals < 112.5] <- "E"
  aspect_cats[valid_vals & aspect_vals >= 112.5 & aspect_vals < 157.5] <- "SE"
  aspect_cats[valid_vals & aspect_vals >= 157.5 & aspect_vals < 202.5] <- "S"
  aspect_cats[valid_vals & aspect_vals >= 202.5 & aspect_vals < 247.5] <- "SW"
  aspect_cats[valid_vals & aspect_vals >= 247.5 & aspect_vals < 292.5] <- "W"
  aspect_cats[valid_vals & aspect_vals >= 292.5 & aspect_vals < 337.5] <- "NW"
  
  aspect_categorical <- aspect_raster
  values(aspect_categorical) <- aspect_cats
  
  return(aspect_categorical)
}

# Strict data cleaning - no tolerance for missing data
clean_prediction_data_strict <- function(prediction_data, predictors) {
  cat("  Strict data cleaning (no NA tolerance)...\n")
  
  # Remove rows with any NA in predictor variables
  complete_rows <- complete.cases(prediction_data[, predictors, drop = FALSE])
  n_removed <- sum(!complete_rows)
  
  if (n_removed > 0) {
    cat(paste("  Removing", n_removed, "rows with NA values\n"))
    prediction_data <- prediction_data[complete_rows, , drop = FALSE]
  }
  
  if (nrow(prediction_data) == 0) {
    cat("  ✗ No complete cases remaining\n")
    return(NULL)
  }
  
  # Handle factor variables
  if ("ASP" %in% predictors && "ASP" %in% names(prediction_data)) {
    prediction_data$ASP <- factor(as.character(prediction_data$ASP), 
                                  levels = c("FLAT", "N", "NE", "E", "SE", "S", "SW", "W", "NW"))
  }
  
  # Handle PDIS_CLASS as factor if present
  if ("PDIS_CLASS" %in% predictors && "PDIS_CLASS" %in% names(prediction_data)) {
    prediction_data$PDIS_CLASS <- as.factor(prediction_data$PDIS_CLASS)
  }
  
  cat(paste("  Clean data:", nrow(prediction_data), "rows\n"))
  return(prediction_data)
}

# Process single DEM file - STRICT MODE
process_deschampsia_tile_strict <- function(dem_file, matches, coastline_path, rf_model, qrf_model) {
  tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
  cat(paste("Processing:", tile_id, "\n"))
  
  tryCatch({
    # Load DEM as template
    altitude_raster <- rast(dem_file)
    names(altitude_raster) <- "ALTITUDE"
    
    if (sum(!is.na(values(altitude_raster))) == 0) {
      cat("  ✗ No valid altitude data\n")
      return(NULL)
    }
    
    template_raster <- altitude_raster
    all_rasters <- list(ALTITUDE = altitude_raster)
    
    # Load climate data first for masking
    cat("  Loading climate envelopes...\n")
    
    # Load TAS_08
    tas_raster <- rast(matches$tas_08)
    names(tas_raster) <- "TAS_08"
    if (!compareGeom(tas_raster, template_raster, stopOnError = FALSE)) {
      tas_raster <- resample(tas_raster, template_raster, method = "bilinear")
    }
    all_rasters$TAS_08 <- tas_raster
    
    # Load PR_01
    pr_raster <- rast(matches$pr_01)
    names(pr_raster) <- "PR_01"
    if (!compareGeom(pr_raster, template_raster, stopOnError = FALSE)) {
      pr_raster <- resample(pr_raster, template_raster, method = "bilinear")
    }
    all_rasters$PR_01 <- pr_raster
    
    # Create climate mask - areas with valid climate data
    climate_mask <- !is.na(tas_raster) & !is.na(pr_raster) & !is.na(altitude_raster)
    cat(paste("  Climate envelope covers", sum(values(climate_mask), na.rm = TRUE), "pixels\n"))
    
    # Check if there's any overlap
    if (sum(values(climate_mask), na.rm = TRUE) == 0) {
      cat("  No overlap with climate envelope - skipping\n")
      return(NULL)
    }
    
    # Apply climate mask to template
    template_raster <- mask(template_raster, climate_mask)
    all_rasters$ALTITUDE <- template_raster
    
    # Load all required variables - NO DUMMY DATA
    variable_map <- list(
      WEXP = "wexp", VDEP = "vdep", EFAH = "efah", MPI = "mpi", 
      WEFF = "weff", LSF = "lsf", ASP = "asp", HillSh = "hillsh",
      MBI = "mbi", PDIS_CLASS = "pdis_class", RSP = "rsp", FLAC = "flac",
      TPW = "tpw", MNCURV = "mncurv", SVF = "svf"
    )
    
    for (var_name in names(variable_map)) {
      match_key <- variable_map[[var_name]]
      
      if (match_key %in% names(matches)) {
        cat(paste("  Loading", var_name, "from", basename(matches[[match_key]]), "\n"))
        
        var_raster <- rast(matches[[match_key]])
        names(var_raster) <- var_name
        
        if (!compareGeom(var_raster, template_raster, stopOnError = FALSE)) {
          var_raster <- resample(var_raster, template_raster, method = "bilinear")
        }
        
        # Apply climate mask
        var_raster <- mask(var_raster, climate_mask)
        
        # Special handling for categorical variables
        if (var_name == "ASP") {
          var_raster <- convert_aspect_to_categorical(var_raster)
        } else if (var_name == "PDIS_CLASS") {
          # Keep PDIS_CLASS as numeric for now, will convert to factor in data frame
          names(var_raster) <- var_name
        }
        
        all_rasters[[var_name]] <- var_raster
      } else {
        # Check if this variable is required
        required_for_presence <- var_name %in% c("WEXP", "VDEP", "EFAH", "MPI", "WEFF", "LSF", "ASP", "HillSh", "MBI", "PDIS_CLASS", "RSP", "FLAC")
        required_for_abundance <- var_name %in% c("TPW", "LSF", "MPI", "ASP", "EFAH", "MBI", "WEXP", "MNCURV", "FLAC", "SVF")
        
        if (required_for_presence || required_for_abundance) {
          cat(paste("  ✗ Missing required variable:", var_name, "\n"))
          return(NULL)
        }
      }
    }
    
    # Calculate coordinates
    cat("  Calculating coordinates...\n")
    coord_rasters <- calculate_coordinates(template_raster)
    all_rasters$Easting <- mask(coord_rasters[[1]], climate_mask)
    all_rasters$Northing <- mask(coord_rasters[[2]], climate_mask)
    
    # Calculate proper distance to coast
    cat("  Calculating distance to coast...\n")
    dist_raster <- calculate_distance_to_coast_proper(template_raster, coastline_path)
    
    if (is.null(dist_raster)) {
      cat("  ✗ Could not calculate distance to coast\n")
      return(NULL)
    }
    
    all_rasters$DISTCOAST <- mask(dist_raster, climate_mask)
    
    # Convert to data frame
    cat("  Converting to data frame...\n")
    combined_rasters <- rast(all_rasters)
    prediction_data <- as.data.frame(combined_rasters, xy = TRUE, na.rm = TRUE)
    
    if (nrow(prediction_data) == 0) {
      cat("  ✗ No valid data after conversion\n")
      return(NULL)
    }
    
    # Check which predictors are available
    available_presence_vars <- intersect(predictors_presence, names(prediction_data))
    available_abundance_vars <- intersect(predictors_abundance, names(prediction_data))
    
    cat(paste("  Available presence predictors:", length(available_presence_vars), "of", length(predictors_presence), "\n"))
    cat(paste("  Available abundance predictors:", length(available_abundance_vars), "of", length(predictors_abundance), "\n"))
    
    # Need minimum number of predictors
    if (length(available_presence_vars) < 10) {
      cat("  ✗ Insufficient presence predictors available\n")
      return(NULL)
    }
    
    # Clean data strictly
    presence_data <- clean_prediction_data_strict(prediction_data, available_presence_vars)
    if (is.null(presence_data)) return(NULL)
    
    # Make predictions
    cat(paste("  Making predictions on", nrow(presence_data), "pixels\n"))
    
    prob_pred <- predict(rf_model, newdata = presence_data, type = "prob")[, "1"]
    binary_pred <- ifelse(prob_pred >= 0.5, 1, 0)
    
    # Abundance predictions for present pixels
    abundance_pred <- rep(NA, nrow(presence_data))
    present_idx <- which(binary_pred == 1)
    
    if (length(present_idx) > 0 && length(available_abundance_vars) >= 8) {
      cat(paste("  Calculating abundance for", length(present_idx), "present pixels\n"))
      abundance_data <- clean_prediction_data_strict(presence_data[present_idx, ], available_abundance_vars)
      
      if (!is.null(abundance_data) && nrow(abundance_data) > 0) {
        abundance_predictions <- predict(qrf_model, newdata = abundance_data, what = 0.5)
        abundance_pred[present_idx[1:length(abundance_predictions)]] <- pmax(pmin(abundance_predictions, 100), 0)
      }
    }
    
    # Create output rasters (using original template before masking for full extent)
    original_template <- rast(dem_file)
    prob_raster <- binary_raster <- abundance_raster <- original_template
    values(prob_raster) <- values(binary_raster) <- values(abundance_raster) <- NA
    
    # Map predictions back to raster
    cell_nums <- cellFromXY(original_template, presence_data[, c("x", "y")])
    valid_cells <- !is.na(cell_nums)
    
    if (sum(valid_cells) > 0) {
      values(prob_raster)[cell_nums[valid_cells]] <- prob_pred[valid_cells]
      values(binary_raster)[cell_nums[valid_cells]] <- binary_pred[valid_cells]
      values(abundance_raster)[cell_nums[valid_cells]] <- abundance_pred[valid_cells]
    }
    
    names(prob_raster) <- "Presence_Probability"
    names(binary_raster) <- "Presence_Binary"
    names(abundance_raster) <- "Abundance_Prediction"
    
    return(list(
      tile_id = tile_id,
      presence_prob = prob_raster,
      presence_binary = binary_raster,
      abundance = abundance_raster,
      climate_mask = climate_mask,
      stats = list(
        pixels_present = sum(binary_pred == 1, na.rm = TRUE),
        total_pixels = nrow(presence_data),
        mean_probability = round(mean(prob_pred, na.rm = TRUE), 3),
        climate_pixels = sum(values(climate_mask), na.rm = TRUE),
        available_presence_vars = length(available_presence_vars),
        available_abundance_vars = length(available_abundance_vars)
      )
    ))
    
  }, error = function(e) {
    cat(paste("✗ ERROR:", e$message, "\n"))
    return(NULL)
  })
}

#### LOAD MODELS ####
train_models <- function() {
  training_file <- file.path(base_dir, "hurdle_training_data_list_Deschampsia_top_predictors.rds")
  if (!file.exists(training_file)) stop("Training data file not found!")
  
  training_data_list <- readRDS(training_file)
  train_data <- training_data_list[[1]]
  
  cat("Training DESCHAMPSIA models with updated predictors...\n")
  
  # Use available predictors from training data
  available_presence <- intersect(predictors_presence, names(train_data))
  available_abundance <- intersect(predictors_abundance, names(train_data))
  
  cat(paste("Training presence model with:", length(available_presence), "variables\n"))
  cat(paste("Training abundance model with:", length(available_abundance), "variables\n"))
  
  rf_clf <- randomForest(
    x = train_data[, available_presence, drop = FALSE],
    y = as.factor(train_data$Presence),
    ntree = 300,
    importance = TRUE
  )
  
  present_data <- train_data[train_data$Presence == "1" | train_data$Presence == 1, ]
  qrf <- quantregForest(
    x = present_data[, available_abundance, drop = FALSE],
    y = as.numeric(present_data$Deschampsia),
    ntree = 300
  )
  
  return(list(rf_classifier = rf_clf, qrf_model = qrf))
}

#### MAIN PROCESSING ####
main_processing <- function() {
  cat("Loading models...\n")
  models <- train_models()
  rf_model <- models$rf_classifier
  qrf_model <- models$qrf_model
  
  # Find all DEM files
  dem_files <- find_all_dem_files(base_dir)
  
  # Coastline file
  coastline_file <- file.path(base_dir, "COAST", "add_coastline_high_res_line_v7_10.shp")
  if (!file.exists(coastline_file)) {
    stop("Coastline file not found - required for proper distance calculation")
  }
  
  successful_tiles <- c()
  failed_tiles <- c()
  
  cat(paste("\n=== PROCESSING", length(dem_files), "DEM FILES ===\n"))
  
  for (i in seq_along(dem_files)) {
    dem_file <- dem_files[i]
    tile_id <- gsub("\\.tif$", "", basename(dem_file), ignore.case = TRUE)
    
    cat(paste("\n=== PROCESSING", i, "OF", length(dem_files), ":", tile_id, "===\n"))
    
    # Find required matching files
    matches <- find_required_files_for_dem(dem_file, base_dir)
    
    # Check strict requirements
    if (!check_strict_requirements_deschampsia(matches)) {
      failed_tiles <- c(failed_tiles, tile_id)
      next
    }
    
    # Process the tile
    result <- process_deschampsia_tile_strict(dem_file, matches, coastline_file, rf_model, qrf_model)
    
    if (!is.null(result)) {
      # Save results
      tile_output_dir <- file.path(deschampsia_output_dir, paste0("tile_", result$tile_id))
      if (!dir.exists(tile_output_dir)) dir.create(tile_output_dir, recursive = TRUE)
      
      writeRaster(result$presence_prob, file.path(tile_output_dir, paste0("presence_prob_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$presence_binary, file.path(tile_output_dir, paste0("presence_binary_", result$tile_id, ".tif")), overwrite = TRUE)
      writeRaster(result$abundance, file.path(tile_output_dir, paste0("abundance_", result$tile_id, ".tif")), overwrite = TRUE)
      
      # Save climate mask for reference
      writeRaster(result$climate_mask, file.path(tile_output_dir, paste0("climate_mask_", result$tile_id, ".tif")), overwrite = TRUE)
      
      successful_tiles <- c(successful_tiles, result$tile_id)
      cat(paste("✓ SUCCESS:", result$tile_id, "- Present pixels:", result$stats$pixels_present, "\n"))
      cat(paste("  Climate pixels:", result$stats$climate_pixels, "Presence vars:", result$stats$available_presence_vars, "\n"))
    } else {
      failed_tiles <- c(failed_tiles, tile_id)
    }
    
    gc(full = TRUE)
  }
  
  # Summary
  cat(paste("\n=== FINAL SUMMARY ===\n"))
  cat(paste("Total files:", length(dem_files), "\n"))
  cat(paste("Successful:", length(successful_tiles), "\n"))
  cat(paste("Failed:", length(failed_tiles), "\n"))
  cat(paste("Success rate:", round(length(successful_tiles)/length(dem_files) * 100, 1), "%\n"))
  cat(paste("Output directory:", deschampsia_output_dir, "\n"))
  
  if (length(failed_tiles) > 0) {
    cat("\n✗ Failed tiles (missing required data):\n")
    for (tile in head(failed_tiles, 10)) {
      cat(paste("  -", tile, "\n"))
    }
    if (length(failed_tiles) > 10) {
      cat(paste("  ... and", length(failed_tiles) - 10, "more\n"))
    }
  }
  
  return(list(successful = successful_tiles, failed = failed_tiles))
}

# Execute main processing
tryCatch({
  result <- main_processing()
  cat("✓ DESCHAMPSIA processing completed successfully\n")
  cat("✓ All predictions clipped to climate envelope (TAS_08 + PR_01)\n")
  cat("✓ Proper distance-to-coast calculations applied\n")
  cat("✓ No dummy data used - only real environmental variables\n")
  cat("✓ Extended predictor set with all required variables\n")
}, error = function(e) {
  cat(paste("✗ Critical error:", e$message, "\n"))
  traceback()
  quit(status = 1)
})

cat("=== DESCHAMPSIA STRICT PROCESSING COMPLETE ===\n")




