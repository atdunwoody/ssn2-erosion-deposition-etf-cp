################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

prefixes <- c(
  "ETF sfm", "ETF lidar",
  # "LM2 sfm", "LPM sfm", "MM_ET sfm",
  # "LM2 lidar", "LPM lidar", "MM_ET lidar",
  "CPF sfm" , "CPF lidar"
  # "ME sfm", "UE sfm", "UW sfm", "UM sfm",
  # "MM lidar", "UE lidar", "UW lidar", "UM lidar"
)
types <- c(
  "net change",
  "deposition",
  "erosion"
           )
segments <- c(
  #20, 
  10 
  # 5
  )  

param_grid <- expand.grid(prefix = prefixes, type = types, segment = segments, stringsAsFactors = FALSE)

corr <- 0.7

# Model formula is stored in outputs folder:
# "ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt"
formula_file_name <- "ssn_formula.txt"

# SSN object already exists for each watershed/combined area, keep TRUE
# If TRUE, the SSN object will be loaded from the existing path:
# ETF/Outputs/LM2_erosion_logtrans/LM2_erosion_logtrans.ssn
# If FALSE, the SSN object will be created with data from:
# Inputs/Individual Watersheds/LM2_erosion_ssn points.gpkg
# Inputs/Streams/streams_100k.gpkg
load_ssn <- FALSE

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
library(progressr)
library(tictoc)

# Show warnings
warnings()


################################################################################
########################### DEFINE INPUT/OUTPUT PATHS ##########################
################################################################################

# Start the overall timer
tic("Total Script Execution Time")

run_ssn_analysis <- function(prefix, type, segment) {
  tryCatch({
    # Early exit conditions based on prefix and type
    if ((grepl("lidar", prefix) && type != "erosion") || prefix %in% c("ME lidar", "MW lidar")) {
      message("Skipping processing for prefix: ", prefix, " and type: ", type)
      return(NULL)
    }
    
    # Determine region and prefix_use
    if (prefix %in% c("CPF sfm", "ME sfm", "MM sfm", "MW sfm", 
                      "UE sfm", "UW sfm", "UM sfm", "CPF lidar", 
                      "ME lidar", "MM lidar", "MW lidar", "UE lidar", 
                      "UW lidar", "UM lidar")) {
      region <- "CPF"
    } else {
      region <- "ETF"
    }
    prefix_use <- switch(prefix,
                         "MM_ET sfm" = "MM sfm",
                         "MM_ET lidar" = "MM lidar",
                         prefix)
    
    message("Processing: ", prefix_use, " - ", type, " with ", segment, "m spacing")
    
    # Define input/output folders based on region, prefix, type, and segment.
    base_input_folder <- file.path(region, "Inputs")
    base_output_folder <- file.path(region, "Outputs")
    segment_input_folder <- file.path(base_input_folder, paste0("Segmented ", segment, "m"))
    segment_output_folder <- file.path(base_output_folder, paste0("Segmented ", segment, "m"))
    
    if (prefix_use %in% c("CPF", "ETF", "CPF sfm", "ETF sfm", "CPF lidar", "ETF lidar")) {
      input_obs <- file.path(segment_input_folder, 
                             "Combined Watersheds", 
                             paste(prefix_use, type, "ssn points.gpkg", sep = " "))
      multiple_ws <- TRUE
    } else {
      input_obs <- file.path(segment_input_folder, 
                             "Individual Watersheds", 
                             paste(prefix_use, type, "ssn points.gpkg", sep = " "))
      multiple_ws <- FALSE
    }
    
    output_folder <- file.path(segment_output_folder, paste0(prefix_use, "_", type, "_logtrans"))
    message("Output folder: ", output_folder)
    
    formula_file <- file.path(output_folder, "ssn_formula.txt")
    if (!file.exists(formula_file)) {
      stop("Formula file does not exist: ", formula_file)
    }
    
    model_formula_str <- readLines(formula_file)
    model_formula <- as.formula(model_formula_str)
    response_var <- all.vars(model_formula)[1]
    
    # Define other parameters
    ssn_path <- file.path(output_folder, paste0(prefix_use, "_", type, "_logtrans.ssn"))
    lsn_out <- file.path(output_folder, "lsn_out")
    input_streams <- file.path(base_input_folder, "Streams", "streams_100k.gpkg")
    output_file <- file.path(output_folder, paste0(response_var, " ", segment, "m_VIF-3_corr", corr, ".txt"))
    
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    
    ################################################################################
    ########################### SSN2 PREPROCESSING ##################################
    ################################################################################
    
    if ((!load_ssn)) {
      # Read spatial data
      CP_streams <- st_read(input_streams)
      CP_obs <- st_read(input_obs)
      print(paste0("CP_obs path: ", input_obs))
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
        snap_tolerance = 10,
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
    
    # Fit the model
    if (multiple_ws && type == "net change") {
      ssn_mod <- ssn_lm(
        formula = model_formula,
        ssn.object = CP_ssn,
        tailup_type = "exponential",
        taildown_type = "none",
        estmethod = "reml",
        euclid_type = "spherical",
        nugget_type = "nugget",
        additive = "afv_flow_accum",
        random = ~ as.factor(ch_watershed)
      )
      model_type <- "ssn_lm"
    }
    else if (multiple_ws) {
      ssn_mod <- ssn_glm(
        formula = model_formula,
        ssn.object = CP_ssn,
        family = "Gamma",
        tailup_type = "exponential",
        taildown_type = "none",
        estmethod = "reml",
        euclid_type = "spherical",
        nugget_type = "nugget",
        additive = "afv_flow_accum",
        random = ~ as.factor(ch_watershed)
      )
      model_type <- "ssn_glm"
    }
    else if (type == "net change") {
      ssn_mod <- ssn_lm(
        formula = model_formula,
        ssn.object = CP_ssn,
        tailup_type = "exponential",
        taildown_type = "none",
        estmethod = "ml",
        euclid_type = "spherical",
        nugget_type = "nugget",
        additive = "afv_flow_accum"
      )
      model_type <- "ssn_lm"
    } else {
      ssn_mod <- ssn_glm(
        formula = model_formula,
        ssn.object = CP_ssn,
        family = "Gamma",
        tailup_type = "exponential",
        taildown_type = "none",
        estmethod = "ml",
        euclid_type = "spherical",
        nugget_type = "nugget",
        additive = "afv_flow_accum"
      )
      model_type <- "ssn_glm"
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
    
    # Add response variable observations to the LOOCV results
    loocv_results$obs <- ssn_get_data(ssn_mod)[[response_var]]
    
    # Specify the output file path
    loocv_results_file <- file.path(output_folder, "loocv_results.csv")
    
    # Write the LOOCV results to a CSV file
    write.csv(loocv_results, loocv_results_file, row.names = FALSE)
    
    tidy_results <- tidy(ssn_mod, conf.int = TRUE)
    tidy_csv <- file.path(output_folder, "tidy_results.csv")
    write.csv(tidy_results, tidy_csv, row.names = FALSE)
    
    varcomp_results <- varcomp(ssn_mod)
    varcomp_csv <- file.path(output_folder, "varcomp_results.csv")
    write.csv(varcomp_results, varcomp_csv, row.names = FALSE)
    
    ################################################################################
    ########################### SSN2 MODEL DIAGNOSTICS #############################
    ################################################################################
    
    # Extract residuals and fitted values
    residuals <- residuals(ssn_mod)         
    fitted_values <- fitted(ssn_mod)       
    
    #----------------------------#
    #   Residuals vs Fitted Plot #
    #----------------------------#
    
    resid_fitted_plot <- ggplot(data.frame(Fitted = fitted_values, Residuals = residuals), 
                                aes(x = Fitted, y = Residuals)) +
      geom_point(color = "blue") +
      geom_hline(yintercept = 0, color = "red") +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "green") +
      labs(title = "Residuals vs Fitted Values",
           x = "Fitted Values",
           y = "Residuals") +
      theme_minimal()
    
    ggsave(filename = file.path(output_folder, "Residuals_vs_Fitted.png"), 
           plot = resid_fitted_plot, width = 8, height = 6)
    
    #----------------------------#
    #         Q-Q Plot           #
    #----------------------------#
    
    qq_plot_gg <- ggplot(data.frame(Residuals = residuals), aes(sample = Residuals)) +
      stat_qq(color = "blue") +
      stat_qq_line(color = "red") +
      labs(title = "Q-Q Plot of Residuals") +
      theme_minimal()
    
    ggsave(filename = file.path(output_folder, "QQ_Plot_Residuals.png"), 
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
    
    ggsave(filename = file.path(output_folder, "Histogram_Residuals.png"), 
           plot = hist_gg, width = 8, height = 6)
    
    #----------------------------#
    #      Statistical Tests     #
    #----------------------------#
    
    stats_test_file <- file.path(output_folder, "Statistical_Tests.txt")
    if (file.exists(stats_test_file)) file.remove(stats_test_file)
    
    # Shapiro-Wilk Test for Normality
    shapiro_test <- shapiro.test(residuals)
    cat("Shapiro-Wilk Test for Normality:\n", file = stats_test_file)
    capture.output(shapiro_test, file = stats_test_file, append = TRUE)
    cat("\n", file = stats_test_file, append = TRUE)
    
    # Moran's I Test for Spatial Autocorrelation
    ssn_data <- ssn_get_data(ssn_mod)
    if (!"geometry" %in% names(ssn_data)) {
      stop("The SSN object does not contain a 'geometry' column.")
    }
    
    coords <- st_coordinates(ssn_data$geometry)
    knn <- knearneigh(coords, k = 12)
    nb <- knn2nb(knn)
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(residuals, lw, zero.policy = TRUE)
    
    cat("Moran's I Test for Spatial Autocorrelation:\n", file = stats_test_file, append = TRUE)
    capture.output(moran_test, file = stats_test_file, append = TRUE)
    cat("\n", file = stats_test_file, append = TRUE)
    
    significance_level <- 0.05
    shapiro_p <- shapiro_test$p.value
    moran_p <- moran_test$p.value
    
    if (shapiro_p > significance_level) {
      shapiro_result <- "passed (Residuals are normally distributed)."
    } else {
      shapiro_result <- paste0("failed (p-value = ", signif(shapiro_p, 4), ").")
    }
    
    if (moran_p > significance_level) {
      moran_result <- "passed (No significant spatial autocorrelation detected)."
    } else {
      moran_result <- paste0("failed (p-value = ", signif(moran_p, 4), ").")
    }
    
    summary_statement <- "Summary of Statistical Tests:\n"
    summary_statement <- paste0(summary_statement, "- Shapiro-Wilk Test for Normality: ", shapiro_result, "\n")
    summary_statement <- paste0(summary_statement, "- Moran's I Test for Spatial Autocorrelation: ", moran_result, "\n")
    
    cat(summary_statement, file = stats_test_file, append = TRUE)
    cat("\n", file = stats_test_file, append = TRUE)
    
    
    # Calculate mean and standard deviation of residuals and response variable
    mean_residuals <- mean(residuals, na.rm = TRUE)
    sd_residuals <- sd(residuals, na.rm = TRUE)
    
    mean_response <- mean(abs(ssn_data[[response_var]]), na.rm = TRUE)
    sd_response <- sd(ssn_data[[response_var]], na.rm = TRUE)
    RMSPE <- loocv_results$stats$RMSPE
    bias <- loocv_results$stats$bias
    RAV <- loocv_results$stats$RAV
    glance_results <- glance(ssn_mod)
    AIC <- glance_results$AIC
    
    response_name <- paste0(response_var)
    # Replace . with _ in the response variable name for file naming
    response_name <- gsub("\\.", "_", response_name)
    model_evaluation_file <- file.path(output_folder, 
                                       paste0(response_name, 
                                              "_corr_", corr, "_evaluation.csv"))
    num_parameters <- length(all.vars(model_formula)) - 1
    
    # num_observations is equal to the number of rows in the CP obs data frame
    num_obs <- as.numeric(nrow(ssn_data))
    
    model_evaluation_data <- data.frame(
      Mean_Residuals = mean_residuals,
      SD_Residuals = sd_residuals,
      Mean_Response = mean_response,
      SD_Response = sd_response,
      RMSPE = RMSPE,
      Bias = bias,
      RAV = RAV,
      AIC = AIC,
      num_parameters = num_parameters,
      num_obs = num_obs
      
    )
    
    write.csv(model_evaluation_data, model_evaluation_file, row.names = FALSE)
    
    message("Finished processing for: ", prefix, " | ", type, " | ", segment)
    
    # Optionally, return some summary information
    return(list(prefix = prefix, type = type, segment = segment, error = NULL))
    
  }, error = function(e) {
    error_msg <- paste("Error for", prefix, "-", type, "-", segment, ":", e$message)
    message(error_msg)
    
    # Return error info instead of NULL so you can record it later
    return(list(prefix = prefix, type = type, segment = segment, error = e$message))
  })
}

################################################################################
########################### PARALLEL PROCESSING ################################
################################################################################

library(future)
library(furrr)
library(progressr)

# Set up the parallel plan (adjust according to your system)
options(future.seed = TRUE)
plan(multisession, workers = 10)

# Optionally, set up progress reporting
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = 1:nrow(param_grid))
  results <- future_pmap(
    param_grid, 
    function(prefix, type, segment) {
    p()  # update progress
    run_ssn_analysis(prefix, type, segment)
  }, 
  .options = furrr_options(seed = TRUE)
  )
})

# Filter out the runs that returned an error
error_runs <- purrr::keep(results, ~ !is.null(.$error))

if (length(error_runs) > 0) {
  cat("The following model runs encountered errors:\n")
  for (err in error_runs) {
    cat("Prefix:", err$prefix, "| Type:", err$type, "| Segment:", err$segment, "\n")
    cat("Error Message:", err$error, "\n\n")
  }
} else {
  cat("No errors encountered.\n")
}


