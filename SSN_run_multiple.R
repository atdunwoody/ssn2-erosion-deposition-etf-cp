################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET sfm", "LM2 sfm", "LPM sfm", "MM_ET"
#                "ET lidar", "LM2 lidar", "LPM lidar", "MM_ET lidar"
# Bennett prefixes: "Bennett sfm", "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm"
#                   "Bennett lidar", "ME lidar", "MM lidar", "MW lidar", "UE lidar", "UW lidar", "UM lidar"
prefixes <- c(
  "Bennett sfm"
  # "ET sfm"
  # "LM2 sfm"
)  # Can define single or multiple prefixes

# Types: "erosion", "deposition", "net"
types <- c(
  # "deposition",
  "erosion"
  #, "net"
)  # Can define single or both erosion, deposition

segment_list <- c(
  # 5, 
  # 10, 
  20
)


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
for (segment in segment_list) {
  for (type in types) {
    for (prefix in prefixes) {
      
      bennett_prefixes <- c("Bennett sfm", "ME sfm", "MM sfm", "MW sfm", 
                            "UE sfm", "UW sfm", "UM sfm", "Bennett lidar", 
                            "ME lidar", "MM lidar", "MW lidar", "UE lidar", 
                            "UW lidar", "UM lidar")
      if (prefix %in% bennett_prefixes) {
        region <- "Bennett"
      } else {
        region <- "ETF"
      }
      
      if (prefix =="MM_ET sfm") {
        prefix <- "MM sfm"
      }
      else if (prefix == "MM_ET lidar") {
        prefix <- "MM lidar"
      }
      
      # Define the base input and output folders
      print(paste0("Processing: ", prefix, " - ", type, " with ", segment, "m spacing"))
      base_input_folder <- file.path(region,"Inputs")
      base_output_folder <- file.path(region,"Outputs")
      
      segment_input_folder <- file.path(base_input_folder, paste0("Segmented ", segment, "m"))
      segment_output_folder <- file.path(base_output_folder, paste0("Segmented ", segment, "m"))
      
      # Determines whether random effect of watershed is included
      if (prefix %in% c("Bennett", "ET", "Bennett sfm", "ET sfm", "Bennett lidar", "ET lidar")) {
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
########################### SSN2 PREPROESSING ##################################
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
    
    
    # Fit the model
    if (multiple_ws && type == "net") {
      ssn_mod <- ssn_lm(
        formula = model_formula,
        ssn.object = CP_ssn,
        tailup_type = "exponential", 
        taildown_type = "none",
        euclid_type = "gaussian",
        nugget_type = "nugget",
        additive = "afv_flow_accum",
        random = ~ as.factor(ch_watershed)
      )
      model_type <- "ssn_lm"
    } else if (multiple_ws) {
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
      model_type <- "ssn_glm"
    } else if (type == "net") {
      ssn_mod <- ssn_lm(
        formula = model_formula,
        ssn.object = CP_ssn,
        tailup_type = "exponential",
        taildown_type = "none",
        euclid_type = "gaussian",
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
        euclid_type = "gaussian",
        nugget_type = "nugget",
        additive = "afv_flow_accum"
      )
      model_type <- "ssn_glm"
    }
    

    cat("\nModel Type:\n", file = output_file, append = TRUE)
    cat(model_type, file = output_file, append = TRUE)
    
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
      loocv_predictions = loocv_results$cv_predict,
      obs_values = response_values,
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
    
    # Create a spatial weights matrix using k-nearest neighbors (k = 100)
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

  }
  }
}
