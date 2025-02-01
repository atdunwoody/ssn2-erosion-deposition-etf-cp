################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

setwd("C:/Users/alextd/Documents/GitHub/ssn2-erosion-deposition-etf-cp")

# ETF prefixes: "ET sfm", "LM2 sfm", "LPM sfm", "MM_ET"
#               "ET lidar", "LM2 lidar", "LPM lidar", "MM_ET lidar"
# Bennett prefixes: "Bennett sfm", "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm"
#                  "Bennett lidar", "ME lidar", "MM lidar", "MW lidar", "UE lidar", "UW lidar", "UM lidar"
prefixes <- c(
  "ET lidar",
  "ET sfm",
  "Bennett sfm",
  "Bennett lidar"
  
)

# Types: "erosion", "deposition", "net"
types <- c(
  "deposition", 
  "erosion",
  "net"
)

segment_list <- c(
  20,
  10,
  5
)

# Model formula is stored in outputs folder:
# "ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt"
formula_file_name <- "ssn_formula.txt"

# If TRUE, the SSN object will be loaded from the existing path
load_ssn <- TRUE

# Bootstrapping parameters
n_bootstrap <- 1000  # Number of bootstrap samples 
set.seed(123)      # For reproducibility

################################################################################
######################### LOAD LIBRARIES ########################################
################################################################################
# Install tictoc if not already installed
if (!requireNamespace("tictoc", quietly = TRUE)) {
  install.packages("tictoc")
}
# Load Required Libraries
library(SSN2)
library(SSNbler)
library(sf)
library(dplyr)
library(purrr)
library(ggplot2)
library(lmtest)
library(spdep)
library(classInt) # For spatial weights if needed
library(knitr)    # For better output formatting (optional)
library(broom)    # For tidy model outputs
library(tictoc)
library(future)
library(furrr)
library(progressr)

warnings()

################################################################################
########################### DEFINE INPUT/OUTPUT PATHS ##########################
################################################################################

plan(multisession, workers = 10)
handlers(global = TRUE)

tic("Total Script Execution Time")

# Initialize a list to keep track of any missing CSV files
missing_csv_files <- c()

# Initialize a list to keep track of any processing failures
failures <- list()

# Loop through each prefix, segment, and type
for (prefix in prefixes) {  
  for (segment in segment_list) {
    for (type in types) {
      # Skip certain combinations based on prefix and type
      if (grepl("lidar", prefix) && type != "erosion") {
        next
      }
      
      # Wrap the entire process in tryCatch to handle errors gracefully
      tryCatch({
        
        bennett_prefixes <- c("Bennett sfm", "ME sfm", "MM sfm", "MW sfm", 
                              "UE sfm", "UW sfm", "UM sfm", "Bennett lidar", 
                              "ME lidar", "MM lidar", "MW lidar", "UE lidar", 
                              "UW lidar", "UM lidar")
        
        if (prefix %in% bennett_prefixes) {
          region <- "Bennett"
        } else {
          region <- "ETF"
        }
        
        if (prefix == "MM_ET sfm") {
          prefix_adjusted <- "MM sfm"
        } else if (prefix == "MM_ET lidar") {
          prefix_adjusted <- "MM lidar"
        } else {
          prefix_adjusted <- prefix
        }
        
        message(paste0("\n\nProcessing: | ", prefix, " | ", type, " | ", segment, "m |\n"))
        
        base_input_folder <- file.path(region, "Inputs")
        base_output_folder <- file.path(region, "Outputs")
        segment_input_folder <- file.path(base_input_folder, paste0("Segmented ", segment, "m"))
        segment_output_folder <- file.path(base_output_folder, paste0("Segmented ", segment, "m"))
        
        # Determines whether random effect of watershed is included
        if (prefix %in% c("Bennett", "ET", "Bennett sfm", "ET sfm", 
                          "Bennett lidar", "ET lidar")) {
          input_obs <- file.path(
            segment_input_folder, 
            "Combined Watersheds", 
            paste(prefix, type, "ssn points.gpkg", sep = " ")
          )
          multiple_ws <- TRUE
        } else {
          input_obs <- file.path(
            segment_input_folder, 
            "Individual Watersheds", 
            paste(prefix, type, "ssn points.gpkg", sep = " ")
          )
          multiple_ws <- FALSE
        }
        
        output_folder <- file.path(
          segment_output_folder, 
          paste0(prefix, "_", type, "_logtrans")
        )
        
        bootstrap_folder <- file.path(output_folder, "bootstrap_results")
        if (!dir.exists(bootstrap_folder)) {
          dir.create(bootstrap_folder, recursive = TRUE)
        }
        
        formula_file <- file.path(
          output_folder, 
          formula_file_name
        )
        
        if (!file.exists(formula_file)) {
          stop(paste("Formula file does not exist:", formula_file))
        }
        model_formula_str <- readLines(formula_file)
        model_formula <- as.formula(model_formula_str)
        response_var <- all.vars(model_formula)[1]
        
        message("Model formula: \n", model_formula_str)
        
        # Define the expected CSV file path following the pattern:
        # [prefix]_bootstrap_results.[type].csv
        # Example: "Bennett sfm_bootstrap_results.deposition.csv"
        expected_csv_file <- file.path(
          bootstrap_folder, 
          paste0(prefix, "_bootstrap_results.", type, ".csv")
        )
        
        # Check if the expected CSV file exists
        if (!file.exists(expected_csv_file)) {
          missing_csv_files <- c(missing_csv_files, expected_csv_file)
        }
        
        ssn_path <- file.path(
          output_folder, 
          paste0(prefix, "_", type, "_logtrans.ssn")
        )
        
        lsn_out <- file.path(output_folder, "lsn_out")
        input_streams <- file.path(
          base_input_folder, 
          "Streams", 
          "streams_100k.gpkg"
        )
        
        if (!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        # [..] Continue with your processing code here
        # For example, load SSN object, perform bootstrapping, save CSVs, etc.
        # Ensure that after processing, the expected CSV is created.
        
      }, error = function(e) {
        message("Error processing combination: ", prefix, ", ", type, ", ", segment, "m")
        message("Error message: ", e$message)
        # Optionally, add to failures list
        failures[[length(failures) + 1]] <- list(prefix = prefix, type = type, segment = segment, error = e$message)
      })
      
    }
  }
}

# After all loops, display the list of missing CSV files
if (length(missing_csv_files) > 0) {
  missing_csv_files_unique <- unique(missing_csv_files)
  message("\nFolders missing the expected bootstrap_results CSV files:")
  print(missing_csv_files_unique)
} else {
  message("\nAll expected bootstrap_results CSV files exist in the output folders.")
}

# Optionally, display any processing failures
if (length(failures) > 0) {
  message("\nSome combinations failed to process:")
  print(failures)
}

toc()  # Stop the timer and display total execution time
