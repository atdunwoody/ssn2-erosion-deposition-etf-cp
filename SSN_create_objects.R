################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET sfm", "LM2 sfm", "LPM sfm", "MM_ET"
#                "ET lidar", "LM2 lidar", "LPM lidar", "MM_ET lidar"
# Bennett prefixes: "Bennett sfm", "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm"
#                   "Bennett lidar", "ME lidar", "MM lidar", "MW lidar", "UE lidar", "UW lidar", "UM lidar"
prefixes <- c(
  "Bennett sfm", "Bennett lidar",
  "ET sfm", "ET lidar",
  "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm",
  "LM2 sfm", "LPM sfm", "MM_ET sfm", 
  "MM lidar", "UE lidar", "UW lidar", "UM lidar",
  "LM2 lidar", "LPM lidar", "MM_ET lidar"
)  # Can define single or multiple prefixes

# Types: "erosion", "deposition", "net"
types <- c(
  "erosion",
  "deposition",
  "net"
)  # Can define single or multiple (erosion, deposition, net)

segment_list <- c(
  5,
  10,
  20
)

load_ssn <- FALSE
overwrite <- TRUE

################################################################################
######################### LOAD LIBRARIES #######################################
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
# Additional Libraries for Parallelization and Progress Updates
library(future)
library(furrr)
library(progressr)

# Show warnings
warnings()

################################################################################
########################### DEFINE INPUT/OUTPUT PATHS ##########################
################################################################################

# Set up parallel plan with 10 workers
plan(multisession, workers = 10)

# Initialize global progress handlers
handlers(global = TRUE)

tic("Total Script Execution Time")

#--------------------------------------------------------------------------------
# Initialize a data frame to store duplicate coordinate info across all runs
#--------------------------------------------------------------------------------
duplicates_df <- data.frame(
  prefix  = character(),
  type    = character(),
  segment = numeric(),
  X       = numeric(),
  Y       = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each type and prefix
for (prefix in prefixes) {
  for (segment in segment_list) {
    for (type in types){
      if (('lidar' %in% prefix) && !('erosion' %in% type)) {
        next
      }
      
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
        prefix_use <- "MM sfm"
      } else if (prefix == "MM_ET lidar") {
        prefix_use <- "MM lidar"
      } else {
        prefix_use <- prefix
      }
      
      # Define the base input and output folders
      print(paste0("Processing: ", prefix_use, " - ", type, " with ", segment, "m spacing"))
      base_input_folder <- file.path(region,"Inputs")
      base_output_folder <- file.path(region,"Outputs")
      
      segment_input_folder <- file.path(base_input_folder, paste0("Segmented ", segment, "m"))
      segment_output_folder <- file.path(base_output_folder, paste0("Segmented ", segment, "m"))
      
      # Determines whether random effect of watershed is included
      if (prefix_use %in% c("Bennett", "ET", "Bennett sfm", "ET sfm", "Bennett lidar", "ET lidar")) {
        input_obs <- file.path(
          segment_input_folder, 
          "Combined Watersheds", 
          paste(prefix_use, type, "ssn points.gpkg", sep = " ")
        )
        multiple_ws <- TRUE
      } else {
        input_obs <- file.path(
          segment_input_folder, 
          "Individual Watersheds", 
          paste(prefix_use, type, "ssn points.gpkg", sep = " ")
        )
        multiple_ws <- FALSE
      }
      
      output_folder <- file.path(
        segment_output_folder, 
        paste0(prefix_use, "_", type, "_logtrans")
      )
      
      output_file <- file.path(
        output_folder, 
        paste0(prefix_use, "_ch_sfm.", type, ".norm_VIF-2_corr0.6.txt")
      )
      
      
      # Define the SSN path using the prefix and type
      ssn_path <- file.path(
        output_folder, 
        paste0(prefix_use, "_", type, "_logtrans.ssn")
      )
      
      # Define the LSN output folder
      lsn_out <- file.path(output_folder, "lsn_out")
      
      # Define the input streams path
      input_streams <- file.path(
        base_input_folder, 
        "Streams", 
        "streams_100k.gpkg"
      )
      
      # Create the output folder if it doesn't exist
      if (!dir.exists(output_folder)) {
        dir.create(output_folder, recursive = TRUE)
      }
      
      ################################################################################
      ########################### SSN2 PREPROCESSING ##################################
      ################################################################################
      
      if ((!load_ssn) && (!file.exists(ssn_path) || overwrite)) {
        # Read spatial data
        CP_streams <- st_read(input_streams)
        CP_obs <- st_read(input_obs)
        
        # Convert lines to spatial network (LSN)
        edges <- lines_to_lsn(
          streams = CP_streams,
          lsn_path = lsn_out,
          snap_tolerance = 1,
          check_topology = TRUE,
          topo_tolerance = 2,
          overwrite = TRUE,
          verbose = TRUE,
          remove_ZM = TRUE
        )
        
        # Assign observation sites to the LSN
        obs <- sites_to_lsn(
          sites = CP_obs,
          edges = edges,
          lsn_path = lsn_out,
          file_name = "obs",
          snap_tolerance = 1,
          save_local = TRUE,
          overwrite = TRUE
        )
        
        # Calculate edge distances
        edges <- updist_edges(
          edges = edges,
          lsn_path = lsn_out,
          calc_length = TRUE
        )
        
        # Assign distances to sites
        site.list <- updist_sites(
          sites = list(obs = obs),
          edges = edges,
          length_col = "Length",
          lsn_path = lsn_out
        )
        
        # Define inflow and AFV columns
        infl_col <- "flow_accum_max"
        segpi <- "flow_accum_PI"
        afv_col <- "afv_flow_accum"
        
        # Calculate AFV for edges
        edges <- afv_edges(
          edges = edges,
          infl_col = infl_col,
          segpi = segpi,
          afv_col = afv_col,
          lsn_path = lsn_out
        )
        
        # Assign AFV to sites
        site.list <- afv_sites(
          sites = site.list,
          edges = edges,
          afv_col = afv_col,
          save_local = TRUE,
          lsn_path = lsn_out
        )
        
        # Assemble SSN object
        CP_ssn <- ssn_assemble(
          edges = edges,
          lsn_path = lsn_out,
          obs_sites = site.list$obs,
          ssn_path = ssn_path,
          import = TRUE,
          overwrite = TRUE
        )
      } else {
        # Load existing SSN object if load_ssn is TRUE
        CP_ssn <- ssn_import(ssn_path, overwrite = TRUE)
        
      }
      
      # Create distance matrix
      ssn_create_distmat(CP_ssn)
      
      # -------------------- Added Code Starts Here --------------------
      # Extract obs_data and write to CSV (with and without coords)
      obs_data <- ssn_get_data(CP_ssn, name = "obs")
      
      # Convert obs_data to a data frame by dropping geometry
      obs_data_df <- st_drop_geometry(obs_data)
      
      # Get coordinate columns
      coords <- st_coordinates(obs_data)
      obs_data_with_coords <- cbind(obs_data_df, coords)
      
      # Define the output CSV file path for obs_data_with_coords
      obs_output_file <- file.path(
        output_folder, 
        paste0(prefix_use, "_obs_data.", type, ".csv")
      )
      
      # Write with coordinates
      write.csv(obs_data_with_coords, obs_output_file, row.names = FALSE)
      
      # Option 1: Using base R
      unique_points_count <- length(unique(obs_data$geometry))
      all_points_unique   <- unique_points_count == nrow(obs_data)
      
      if (all_points_unique) {
        message("Unique.")
      } else {
        message("Duplicates found.")
      }
      
      # Check if X,Y values are unique
      dup_coords <- obs_data_with_coords[duplicated(obs_data_with_coords[, c("X", "Y")]), ]
      if (nrow(dup_coords) > 0) {
        duplicates_df <- rbind(
          duplicates_df,
          data.frame(
            prefix  = rep(prefix, nrow(dup_coords)),
            type    = rep(type, nrow(dup_coords)),
            segment = rep(segment, nrow(dup_coords)),
            X       = dup_coords$X,
            Y       = dup_coords$Y,
            stringsAsFactors = FALSE
          )
        )
      }
      # -------------------- Added Code Ends Here ----------------------
      
    } # End of prefix loop
  } # End of type loop
}
################################################################################
########################### END OF SCRIPT ######################################
################################################################################

# Print information about any duplicated coordinates found across all runs
if (nrow(duplicates_df) > 0) {
  cat("\nDuplicates were found in the following prefix/type/segment combos:\n")
  print(duplicates_df)
} else {
  cat("\nNo duplicate X,Y coordinates were found.\n")
}

# End the overall timer and print it
total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", total_time$toc - total_time$tic, "seconds\n")
