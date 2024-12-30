################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET", "LM2", "LPM", "MM_ET"
# Bennett prefixes: "Bennett", "ME", "MM", "MW", "UE", "UW", "UM"
prefixes <- c("UM_test")  # Can define single or multiple prefixes

# Types: "erosion", "deposition", "net"
types <- c("erosion")  # Can define single or both erosion, deposition

# Model formula is stored in outputs folder:
#"ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt"
formula_file_name <- "ssn_formula.txt"

# SSN object already exists for each watershed/combined area, keep TRUE

# If TRUE, the SSN object will be loaded from the existing path
# ETF/Outputs/LM2_erosion_logtrans/LM2_erosion_logtrans.ssn

# If FALSE, the SSN object will be created with data from:
# Inputs/Individual Watersheds/LM2_erosion_ssn points.gpkg
# Inputs/Streams/streams_100k.gpkg
load_ssn <- TRUE

# Bootstrapping parameters
n_bootstrap <- 10  # Number of bootstrap samples (adjust as needed)
set.seed(123)         # For reproducibility

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

# Loop through each type and prefix
for (type in types) {
  for (prefix in prefixes) {
    
    bennett_prefixes <- c("Bennett", "ME", "MM", "MW", "UE", "UW", "UM", "UM_test")
    if (prefix %in% bennett_prefixes) {
      region <- "Bennett"
    } else {
      region <- "ETF"
    }
    
    if (prefix =="MM_ET") {
      prefix <- "MM"
    }
    
    # Define the base input and output folders
    print(paste0("Processing: ", prefix, " - ", type))
    base_input_folder <- file.path(region,"Inputs")
    base_output_folder <- file.path(region,"Outputs")
    
    # Determines whether random effect of watershed is included
    if (prefix %in% c("Bennett", "ET")) {
      input_obs <- file.path(
        base_input_folder, 
        "Combined Watersheds", 
        paste(prefix, type, "ssn points.gpkg", sep = " ")
      )
      multiple_ws <- TRUE
    } else {
      input_obs <- file.path(
        base_input_folder, 
        "Individual Watersheds", 
        paste(prefix, type, "ssn points.gpkg", sep = " ")
      )
      multiple_ws <- FALSE
    }
    
    output_folder <- file.path(
      base_output_folder, 
      paste0(prefix, "_", type, "_logtrans")
    )
    
    output_file <- file.path(
      output_folder, 
      paste0(prefix, "_ch_sfm.", type, ".norm_VIF-2_corr0.6.txt")
    )
    
    # SSN2 formula can be changed in the formula file
    formula_file <- file.path(
      output_folder, 
      formula_file_name
    )
    
    if (!file.exists(formula_file)) {
      stop(paste("Formula file does not exist:", formula_file))
    }
    model_formula_str <- readLines(formula_file)
    model_formula <- as.formula(model_formula_str)
    
    # Define the SSN path using the prefix and type
    ssn_path <- file.path(
      output_folder, 
      paste0(prefix, "_", type, "_logtrans.ssn")
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
    
    if (!load_ssn) {
      # Read spatial data
      CP_streams <- st_read(input_streams)
      CP_obs <- st_read(input_obs)
      
      # Uncomment if you encounter an error about Multilinestrings
      # CP_streams <- st_cast(CP_streams, "LINESTRING")
      
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
        snap_tolerance = 100,
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
    # Extract obs_data and write to CSV
    obs_data <- ssn_get_data(CP_ssn, name = "obs")
    
    # Define the output CSV file path for obs_data
    obs_output_file <- file.path(
      output_folder, 
      paste0(prefix, "_obs_data.", type, ".csv")
    )
    
    # Convert obs_data to a data frame by dropping geometry
    obs_data_df <- st_drop_geometry(obs_data)
    
    # Write the obs_data to CSV
    write.csv(obs_data_df, obs_output_file, row.names = FALSE)
    
    # Optional: If you want to include coordinates, you can do the following:
    # coords <- st_coordinates(obs_data)
    # obs_data_with_coords <- cbind(obs_data_df, coords)
    # write.csv(obs_data_with_coords, obs_output_file, row.names = FALSE)
    # -------------------- Added Code Ends Here --------------------
    
    ################################################################################
    ########################### BOOTSTRAPPING FUNCTION ##############################
    ################################################################################
    
    bootstrap_model <- function(ssn_obj, n_boot, model_formula, multiple_ws) {
      
      message("Fitting model with ", n_boot, " bootstrap samples...")
      
      # Extract the observed data
      obs_data <- ssn_get_data(ssn_obj, name = "obs")
      
      # Define the bootstrap iteration function
      bootstrap_iteration <- function(i) {
        # Resample with replacement
        boot_indices <- sample(1:nrow(obs_data), replace = TRUE)
        
        boot_data <- obs_data[boot_indices, ]
        
        # Assign a unique ID to each observation in the bootstrap sample
        boot_data <- boot_data %>% 
          mutate(
            unique_id = i * 1e6 + row_number(),  # Unique ID generator
            pid = unique_id                     # Assign unique_id to pid
          ) %>%
          st_as_sf()
        
        locID <- boot_data$locID
        
        # Update 'netgeom' to reflect the new uid while keeping the last number unchanged
        boot_data$netgeom <- mapply(function(geom, uid) {
          sub("(.*\\s)\\d+\\s(\\d+)\\)$", paste0("\\1", uid, " \\2)"), geom)
        }, boot_data$netgeom, boot_data$unique_id)
        
        
        # Remove the temporary unique_id column
        boot_data <- boot_data %>% select(-unique_id)
        
        # Reset row names to ensure they are unique
        rownames(boot_data) <- NULL  # Remove existing row names
        
        # -------------------- Added Code Starts Here --------------------
        # Convert boot_data to a data frame by dropping geometry and add bootstrap_rep
        boot_data_df <- st_drop_geometry(boot_data) %>%
          mutate(bootstrap_rep = i)
        # -------------------- Added Code Ends Here --------------------
        
        #save bootdata from iteration i to csv
        write.csv(boot_data_df, paste0("boot_data_", i, ".csv"))
        
        
        # Continue with SSN update and model fitting
        ssn_updated <- ssn_put_data(boot_data, ssn_obj, name = "obs", resize_data = TRUE)
        ssn_create_distmat(ssn_updated)
        
        response_var = all.vars(model_formula)[1]
        

        # Fit the model
        if (multiple_ws && type == "net") {
          ssn_mod_boot <- ssn_lm(
            formula = model_formula,
            ssn.object = ssn_updated,
            tailup_type = "exponential", 
            taildown_type = "none",
            euclid_type = "gaussian",
            nugget_type = "nugget",
            additive = "afv_flow_accum",
            random = ~ as.factor(ch_watershed)
          )
          model_type <- "ssn_lm"
        }
        else if (multiple_ws) {
          ssn_mod_boot <- ssn_glm(
            formula = model_formula,
            ssn.object = ssn_updated,
            family = "Gamma",
            tailup_type = "exponential",
            taildown_type = "none",
            euclid_type = "gaussian",
            nugget_type = "nugget",
            additive = "afv_flow_accum",
            random = ~ as.factor(ch_watershed)
          )
          model_type <- "ssn_glm"
        }
        else if (type == "net") {
          ssn_mod_boot <- ssn_lm(
            formula = model_formula,
            ssn.object = ssn_updated,
            tailup_type = "exponential",
            taildown_type = "none",
            euclid_type = "gaussian",
            nugget_type = "nugget",
            additive = "afv_flow_accum"
          )
          model_type <- "ssn_lm"
        } 
        else {
          ssn_mod_boot <- ssn_glm(
            formula = model_formula,
            ssn.object = ssn_updated,
            family = "Gamma",
            tailup_type = "exponential",
            taildown_type = "none",
            euclid_type = "gaussian",
            nugget_type = "nugget",
            additive = "afv_flow_accum"
          )
          model_type <- "ssn_glm"
        }
        
        # Extract tidy coefficients with confidence intervals
        tidy_mod <- tidy(ssn_mod_boot, conf.int = TRUE)
        tidy_mod$bootstrap_rep <- i  # Add bootstrap replicate identifier
        
        # Extract variance components
        varcomp_mod <- varcomp(ssn_mod_boot)
        varcomp_mod$bootstrap_rep <- i
        
        # -------------------- Added Code Starts Here --------------------
        # Return boot_data_df along with model results
        return(list(tidy_mod = tidy_mod, varcomp_mod = varcomp_mod, boot_data = boot_data_df, model_type = model_type))
        # -------------------- Added Code Ends Here --------------------
      }
      
      # Use furrr's future_map with progressr for parallel processing and progress updates
      with_progress({
        p <- progressor(along = 1:n_boot)
        
        bootstrap_results <- future_map(1:n_boot, function(i) {
          
          bootstrap_iteration(i)
        }, .options = furrr_options(seed = TRUE))
        
        p()
      })
      
      # Extract results
      results_list <- map(bootstrap_results, "tidy_mod")
      varcomp_list <- map(bootstrap_results, "varcomp_mod")
      boot_data_list <- map(bootstrap_results, "boot_data")  # Collect boot_data
      model_type <- bootstrap_results[[1]]$model_type
      # Combine all bootstrap results into a single data frame
      bootstrap_results_df <- bind_rows(results_list)
      
      # Combine all bootstrap varcomp results into a single data frame
      bootstrap_varcomp_df <- bind_rows(varcomp_list)
      
      # Combine all bootstrap_data into a single data frame
      bootstrap_data_df <- bind_rows(boot_data_list)
      
      return(list(
        bootstrap_results = bootstrap_results_df, 
        bootstrap_varcomp = bootstrap_varcomp_df,
        bootstrap_data = bootstrap_data_df,
        model_type = model_type# Include boot_data
      ))
      
    }
    
    
    ################################################################################
    ########################### BOOTSTRAPPING PROCESS ################################
    ################################################################################
    
    # Perform bootstrapping
    cat("\nStarting Bootstrapping...\n", file = output_file, append = TRUE)
    bootstrap_output <- bootstrap_model(
      ssn_obj = CP_ssn,
      n_boot = n_bootstrap,
      model_formula = model_formula,
      multiple_ws = multiple_ws
    )
    
    # Extract the individual components from the returned list
    bootstrap_results <- bootstrap_output$bootstrap_results
    bootstrap_varcomp <- bootstrap_output$bootstrap_varcomp
    bootstrap_data <- bootstrap_output$bootstrap_data  # Extract boot_data
    model_type <- bootstrap_output$model_type
    # Save bootstrap summary statistics
    bootstrap_summary <- bootstrap_results %>%
      group_by(term) %>%
      summarise(
        estimate_mean = mean(estimate, na.rm = TRUE),
        estimate_sd = sd(estimate, na.rm = TRUE),
        conf_low = quantile(estimate, 0.025, na.rm = TRUE),
        conf_high = quantile(estimate, 0.975, na.rm = TRUE),
        p_value = t.test(estimate, mu = 0)$p.value
      )
    
    # Write bootstrap summary to the output file
    cat("\nBootstrap Summary:\n", file = output_file, append = TRUE)
    capture.output(paste0("Model Type: ", model_type), file = output_file, append = TRUE)
    capture.output(print(bootstrap_summary), file = output_file, append = TRUE)
    
    # -------------------- Added Code Starts Here --------------------
    # Define the output CSV file path for bootstrap results
    bootstrap_results_file <- file.path(
      output_folder, 
      paste0(prefix, "_bootstrap_results.", type, ".csv")
    )
    
    # Write bootstrap_results to CSV
    write.csv(bootstrap_results, bootstrap_results_file, row.names = FALSE)
    
    # Define the output CSV file path for variance components
    bootstrap_varcomp_file <- file.path(
      output_folder, 
      paste0(prefix, "_varcomp_results.", type, ".csv")
    )
    
    # Write bootstrap_varcomp to CSV
    write.csv(bootstrap_varcomp, bootstrap_varcomp_file, row.names = FALSE)
    
    # Define the output CSV file path for bootstrap data
    bootstrap_data_file <- file.path(
      output_folder, 
      paste0(prefix, "_bootstrap_data.", type, ".csv")
    )
    
    # Write bootstrap_data to CSV
    write.csv(bootstrap_data, bootstrap_data_file, row.names = FALSE)
    # -------------------- Added Code Ends Here --------------------
    
    cat("\nBootstrapping Completed.\n", file = output_file, append = TRUE)
    
  } # End of prefix loop
} # End of type loop

################################################################################
########################### END OF SCRIPT ######################################
################################################################################

# End the overall timer and print it
total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", total_time$toc - total_time$tic, "seconds\n")
       