################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET", "LM2", "LPM", "MM_ET"
# Bennett prefixes: "Bennett", "ME", "MM", "MW", "UE", "UW", "UM"
prefixes <- c("ME")  # Can define single or multiple prefixes

# Types: "erosion", "deposition"
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
n_bootstrap <- 3  # Number of bootstrap samples
set.seed(123)       # For reproducibility

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
library(tictoc)   # For timing the script

#show warnings
warnings()

################################################################################
########################### DEFINE INPUT/OUTPUT PATHS ##########################
################################################################################

tic("Total Script Execution Time")

# Loop through each type and prefix
for (type in types) {
  for (prefix in prefixes) {
    
    bennett_prefixes <- c("Bennett", "ME", "MM", "MW", "UE", "UW", "UM")
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
    
    ################################################################################
    ########################### BOOTSTRAPPING FUNCTION ##############################
    ################################################################################
    
    bootstrap_model <- function(ssn_obj, n_boot, model_formula, multiple_ws) {
      
      print(paste("Fitting model with", n_boot, "bootstrap samples..."))
      
      # Extract the observed data
      obs_data <- ssn_get_data(ssn_obj, name = "obs")
      
      # Initialize a list to store results
      results_list <- vector("list", n_boot)
      
      varcomp_list <- vector("list", n_boot)
      
      for (i in 1:n_boot) {
        print(paste("Bootstrap iteration:", i))
        
        # Resample with replacement
        boot_indices <- sample(1:nrow(obs_data), replace = TRUE)
        
        boot_data <- obs_data[boot_indices, ]
        
        # Assign a unique ID to each observation in the bootstrap sample
        boot_data <- boot_data %>% 
          mutate(
            unique_id = i * 1e6 + row_number(),  # Unique ID generator
            pid = unique_id,                     # Assign unique_id to pid
            locID = unique_id                    # Assign unique_id to locID
          ) %>%
          st_as_sf()
        
        # Update 'netgeom' to reflect the new pid and locID
        boot_data$netgeom <- mapply(function(geom, uid) {
          sub("(.*\\s)(\\d+)\\s(\\d+)\\)$", paste0("\\1", uid, " ", uid, ")"), geom, uid)
        }, boot_data$netgeom, boot_data$unique_id)
        
        # Remove the temporary unique_id column
        boot_data <- boot_data %>% select(-unique_id)
        
        # Reset row names to ensure they are unique
        rownames(boot_data) <- NULL  # Remove existing row names
        
        # Continue with SSN update and model fitting
        ssn_updated <- ssn_put_data(boot_data, ssn_obj, name = "obs", resize_data = TRUE)
        ssn_create_distmat(ssn_updated)
        
        # Fit the model
        if (multiple_ws) {
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
        } else {
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
        }
        
        # Extract tidy coefficients with confidence intervals
        tidy_mod <- tidy(ssn_mod_boot, conf.int = TRUE)
        tidy_mod$bootstrap_rep <- i  # Add bootstrap replicate identifier
        
        # Store the results
        results_list[[i]] <- tidy_mod
        
        varcomp_mod <- varcomp(ssn_mod_boot)
        varcomp_mod$bootstrap_rep <- i
        
        varcomp_list[[i]] <- varcomp_mod
      }
      
      # Combine all bootstrap results into a single data frame
      bootstrap_results <- bind_rows(results_list)
      
      # Combine all bootstrap varcomp results into a single data frame
      bootstrap_varcomp <- bind_rows(varcomp_list)
      
      
      
      return(list(bootstrap_results = bootstrap_results, bootstrap_varcomp = bootstrap_varcomp))
      
    }
    
    ################################################################################
    ########################### SSN2 MODEL FITTING #################################
    ################################################################################
    
    
    # Save input parameters to the output file
    cat("Model Formula:\n", file = output_file)
    capture.output(model_formula, file = output_file, append = TRUE)
    cat("\nSSN Path:\n", file = output_file, append = TRUE)
    cat(ssn_path, file = output_file, append = TRUE)
    cat("\nLSN Output Path:\n", file = output_file, append = TRUE)
    cat(lsn_out, file = output_file, append = TRUE)
    cat("\n", file = output_file, append = TRUE)
    
    
    if (multiple_ws) {
      # SSN2 model fitting with random effect of watershed     
      ssn_mod <- ssn_glm(
        formula = model_formula,
        ssn.object = CP_ssn,
        family = "Gamma",
        tailup_type = "exponential", 
        taildown_type = "none",
        euclid_type = "gaussian",
        nugget_type = "nugget",
        additive = "afv_flow_accum",
        random = ~ as.factor(ch_watershed)
      )
    } else {
      # SSN2 model fitting for single watersheds
      ssn_mod <- ssn_glm(
        formula = model_formula,
        ssn.object = CP_ssn,
        family = "Gamma",
        tailup_type = "exponential",
        taildown_type = "none",
        euclid_type = "gaussian",
        nugget_type = "nugget",
        additive = "afv_flow_accum"
      )
    }
    
    
    # Append the test results to the file
    cat("\nsummary(ssn_mod)\n", file = output_file, append = TRUE)
    capture.output(summary(ssn_mod), file = output_file, append = TRUE)
    
    cat("\ntidy(ssn_mod, conf.int = TRUE):\n", file = output_file, append = TRUE)
    capture.output(print(tidy(ssn_mod, conf.int = TRUE), n = Inf), file = output_file, append = TRUE)
    
    cat("\nvarcomp(ssn_mod):\n", file = output_file, append = TRUE)
    capture.output(varcomp(ssn_mod), file = output_file, append = TRUE)
    
    cat("\nloocv(ssn_mod):\n", file = output_file, append = TRUE)
    capture.output(loocv(ssn_mod), file = output_file, append = TRUE)
    
    cat("\nglance(ssn_mod):\n", file = output_file, append = TRUE)
    capture.output(glance(ssn_mod), file = output_file, append = TRUE)
    
    cat("\n", file = output_file, append = TRUE)
    
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
    capture.output(print(bootstrap_summary), file = output_file, append = TRUE)
    
    # Optionally, save the full bootstrap results to a separate CSV
    bootstrap_output_file <- file.path(
      output_folder, 
      paste0(prefix, "_bootstrap_results.", type, ".csv")
    )
    write.csv(bootstrap_results, bootstrap_output_file, row.names = FALSE)
    
    
    varcomp_output_file <- file.path(
      output_folder, 
      paste0(prefix, "_varcomp_results.", type, ".csv")
    )
    write.csv(bootstrap_varcomp, varcomp_output_file, row.names = FALSE)
    
    cat("\nBootstrapping Completed.\n", file = output_file, append = TRUE)
    
  } # End of prefix loop
} # End of type loop

################################################################################
########################### END OF SCRIPT ######################################
################################################################################
# End the overall timer and print it
total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", total_time$toc - total_time$tic, "seconds\n")
