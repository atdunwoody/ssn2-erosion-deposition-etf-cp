################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET", "LM2", "LPM", "MM_ET"
# Bennett prefixes: "Bennett", "ME", "MM", "MW", "UE", "UW", "UM"
prefixes <- c("LM2")  # Can define single or multiple prefixes

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
load_ssn <- FALSE

# set seed for reproducibility
options(future.seed = TRUE)
# Add this near the beginning of your script, before any parallel processing

################################################################################
######################### LOAD LIBRARIES #######################################
################################################################################

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

# Additional Libraries for Parallelization and Progress Updates
library(future)
library(furrr)
library(progressr)
library(tictoc)

# Show warnings
warnings()

################################################################################
########################### DEFINE INPUT/OUTPUT PATHS ##########################
################################################################################

# Start the overall timer
tic("Total Script Execution Time")

# Set up parallel plan with 10 workers
plan(multisession, workers = 10)


# Initialize global progress handlers
handlers(global = TRUE)

# Ignore seed warnings


# Create a data frame of all combinations of prefixes and types
combinations <- expand.grid(prefix = prefixes, type = types, stringsAsFactors = FALSE)

# Function to process each combination
process_combination <- function(prefix, type, formula_file_name, p) {
  

  # Start timer for each prefix-type combination
  tic(paste("Processing:", prefix, "-", type))
  
  # Initialize a list to store messages or results
  result <- list()
  
  # Determine the region based on prefix
  bennett_prefixes <- c("Bennett", "ME", "MM", "MW", "UE", "UW", "UM")
  if (prefix %in% bennett_prefixes) {
    region <- "Bennett"
  } else {
    region <- "ETF"
  }
  
  if (prefix == "MM_ET") {
    prefix <- "MM"
  }
  
  # Define the base input and output folders
  message(paste0("Processing: ", prefix, " - ", type))
  base_input_folder <- file.path(region, "Inputs")
  base_output_folder <- file.path(region, "Outputs")
  
  # Determine if random effect of watershed is included
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
  
  # Start timer for SSN2 Preprocessing
  tic("SSN2 Preprocessing")
  
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
      sites = list(
        obs = obs
        ),
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
  
  print("Creating distance matrix")
  # Create distance matrix
  ssn_create_distmat(CP_ssn)
  print("FInished creating distance matrix")
  
  # End timer for SSN2 Preprocessing and log it to the output file
  preproc_time <- toc(log = TRUE, quiet = TRUE)
  cat("\nTiming for SSN2 Preprocessing:", 
      round(preproc_time$toc - preproc_time$tic, 2), 
      "seconds\n", 
      file = output_file, append = TRUE)
  
  ################################################################################
  ########################### SSN2 MODEL FITTING #################################
  ################################################################################
  
  # Start timer for SSN2 Model Fitting
  tic("SSN2 Model Fitting")
  
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
  
  # End timer for SSN2 Model Fitting and log it to the output file
  model_fit_time <- toc(log = TRUE, quiet = TRUE)
  cat("\nTiming for SSN2 Model Fitting:", 
      round(model_fit_time$toc - model_fit_time$tic, 2), 
      "seconds\n", 
      file = output_file, append = TRUE)
  
  # Save input parameters to the output file
  cat("Model Formula:\n", file = output_file)
  capture.output(model_formula, file = output_file, append = TRUE)
  cat("\nSSN Path:\n", file = output_file, append = TRUE)
  cat(ssn_path, file = output_file, append = TRUE)
  cat("\nLSN Output Path:\n", file = output_file, append = TRUE)
  cat(lsn_out, file = output_file, append = TRUE)
  cat("\n", file = output_file, append = TRUE)
  
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
  
  # Perform Leave-One-Out Cross-Validation
  loocv_results <- loocv(ssn_mod, cv_predict = TRUE, se.fit = TRUE)
  
  # Extract the response variable from the formula
  response_var <- all.vars(model_formula)[1]
  
  # Get response variable values from the SSN object
  response_values <- ssn_get_data(CP_ssn)[[response_var]]
  
  # Combine the LOOCV results and response values into a single data frame
  loocv_results_df <- data.frame(
    loocv_predictions = loocv_results$loocv_predictions,
    response_obs_values = observed_values,
    se_fit = loocv_results$se.fit
  )
  
  # Specify the output file path
  loocv_results_file <- file.path(output_folder, "LOOCV_Results.csv")
  
  # Write the LOOCV results to a CSV file
  write.csv(loocv_results_df, loocv_results_file, row.names = FALSE)
  
  ################################################################################
  ########################### SSN2 MODEL DIAGNOSTICS #############################
  ################################################################################
  
  # Extract residuals and fitted values
  residuals <- residuals(ssn_mod)         
  fitted_values <- fitted(ssn_mod)       
  
  #----------------------------#
  #   Residuals vs Fitted Plot #
  #----------------------------#
  
  # Create Residuals vs Fitted plot
  resid_fitted_plot <- ggplot(data.frame(Fitted = fitted_values, Residuals = residuals), 
                              aes(x = Fitted, y = Residuals)) +
    geom_point(color = "blue") +
    geom_hline(yintercept = 0, color = "red") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "green") +
    labs(title = "Residuals vs Fitted Values",
         x = "Fitted Values",
         y = "Residuals") +
    theme_minimal()
  
  # Save the plot to the output folder
  ggsave(filename = file.path(output_folder, "Residuals_vs_Fitted.png"), 
         plot = resid_fitted_plot, width = 8, height = 6)
  
  
  #----------------------------#
  #         Q-Q Plot           #
  #----------------------------#
  
  # ggplot2 Q-Q Plot
  qq_plot_gg <- ggplot(data.frame(Residuals = residuals), aes(sample = Residuals)) +
    stat_qq(color = "blue") +
    stat_qq_line(color = "red") +
    labs(title = "Q-Q Plot of Residuals") +
    theme_minimal()
  
  # Save the ggplot2 Q-Q plot
  ggsave(filename = file.path(output_folder, "QQ_Plot_Residuals_ggplot2.png"), 
         plot = qq_plot_gg, width = 8, height = 6)
  
  #----------------------------#
  #       Histogram Plot       #
  #----------------------------#
  
  hist_gg <- ggplot(data.frame(Residuals = residuals), aes(x = Residuals)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, 
                   fill = "lightblue", color = "black") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(residuals), sd = sd(residuals)),
                  color = "red", linewidth = 1) +
    labs(title = "Histogram of Residuals with Normal Curve",
         x = "Residuals",
         y = "Density") +
    theme_minimal()
  
  # Save the ggplot2 Histogram
  ggsave(filename = file.path(output_folder, "Histogram_Residuals_ggplot2.png"), 
         plot = hist_gg, width = 8, height = 6)
  
  #----------------------------#
  #      Statistical Tests     #
  #----------------------------#
  
  # Define the path to the output statistics file
  stats_test_file <- file.path(output_folder, "Statistical_Tests.txt")
  
  # Initialize the output file by overwriting any existing content
  if (file.exists(stats_test_file)) file.remove(stats_test_file)
  
  #-------------------------------#
  # Shapiro-Wilk Test for Normality
  #-------------------------------#
  
  # Perform Shapiro-Wilk Test on residuals
  shapiro_test <- shapiro.test(residuals)
  
  # Append the test name to the file
  cat("Shapiro-Wilk Test for Normality:\n", file = stats_test_file)
  
  # Append the test results to the file
  capture.output(shapiro_test, file = stats_test_file, append = TRUE)
  
  # Optionally, add a newline for readability
  cat("\n", file = stats_test_file, append = TRUE)
  
  #------------------------------------------#
  # Moran's I Test for Spatial Autocorrelation
  #------------------------------------------#
  
  # Extract spatial coordinates from the SSN object
  ssn_data <- ssn_get_data(ssn_mod)
  
  # Check if 'geometry' column exists
  if (!"geometry" %in% names(ssn_data)) {
    stop("The SSN object does not contain a 'geometry' column.")
  }
  
  # Extract coordinates
  coords <- st_coordinates(ssn_data$geometry)
  
  # Create a spatial weights matrix using k-nearest neighbors (k = 12)
  knn <- knearneigh(coords, k = 12)
  nb <- knn2nb(knn)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Compute Moran's I Test
  moran_test <- moran.test(residuals, lw, zero.policy = TRUE)
  
  # Append the test name to the file
  cat("Moran's I Test for Spatial Autocorrelation:\n", file = stats_test_file, append = TRUE)
  
  # Append the test results to the file
  capture.output(moran_test, file = stats_test_file, append = TRUE)
  
  # Optionally, add a newline for readability
  cat("\n", file = stats_test_file, append = TRUE)
  
  #-----------------------------------------#
  # Summary Statement Based on Test Results
  #-----------------------------------------#
  
  # Define the significance level
  significance_level <- 0.05
  
  # Extract p-values from the test results
  shapiro_p <- shapiro_test$p.value
  moran_p <- moran_test$p.value
  
  # Initialize a variable to hold the summary statement
  summary_statement <- ""
  
  # Determine pass/fail for Shapiro-Wilk Test
  if (shapiro_p > significance_level) {
    shapiro_result <- "passed (Residuals are normally distributed)."
  } else {
    shapiro_result <- paste0("failed (p-value = ", signif(shapiro_p, 4), ").")
  }
  
  # Determine pass/fail for Moran's I Test
  if (moran_p > significance_level) {
    moran_result <- "passed (No significant spatial autocorrelation detected)."
  } else {
    moran_result <- paste0("failed (p-value = ", signif(moran_p, 4), ").")
  }
  
  # Construct the summary statement
  summary_statement <- "Summary of Statistical Tests:\n"
  summary_statement <- paste0(summary_statement, "- Shapiro-Wilk Test for Normality: ", shapiro_result, "\n")
  summary_statement <- paste0(summary_statement, "- Moran's I Test for Spatial Autocorrelation: ", moran_result, "\n")
  
  # Append the summary statements to the file
  cat(summary_statement, file = stats_test_file, append = TRUE)
  
  # Add a final newline for readability
  cat("\n", file = stats_test_file, append = TRUE)
  
  print(paste0("Finished SSN test and statistical evaluation for: ", 
               prefix, 
               " - ", 
               type)
  )
  p()
  
}

################################################################################
########################### Parallel Processing #################################
################################################################################

# Initialize progress handlers
with_progress({
  p <- progressor(along = 1:nrow(combinations))
  
  # Use furrr's future_pmap to apply the process_combination function in parallel
  results <- furrr::future_pmap(
    combinations,
    function(prefix, type) {
      process_combination(prefix, type, formula_file_name, p)
    },
    .options = furrr_options(seed = TRUE)
  )
})

################################################################################
########################### END OF SCRIPT ######################################
################################################################################

# End the overall timer and print it
total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", 
    round(total_time$toc - total_time$tic, 2), 
    "seconds\n")
