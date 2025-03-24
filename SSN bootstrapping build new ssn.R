################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# Print working directory
getwd()

# ETF prefixes: "ET sfm", "LM2 sfm", "LPM sfm", "MM_ET"
#               "ET lidar", "LM2 lidar", "LPM lidar", "MM_ET lidar"
# CPF prefixes: "CPF sfm", "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm"
#              "CPF lidar", "ME lidar", "MM lidar", "MW lidar", "UE lidar", "UW lidar", "UM lidar"
prefixes <- c(
  # "ET lidar",
  # "ET sfm",
  # "CPF sfm",
  "CPF lidar"
)

# Types: "erosion", "deposition", "net"
types <- c(
  # "deposition", 
  # "net",
  "erosion"
)

segment_list <- c(
  20
  # 10
  #,5
)

# Model formula is stored in outputs folder:
# "ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt"
formula_file_name <- "ssn_formula.txt"

# If TRUE, the SSN object will be loaded from the existing path
load_ssn <- TRUE

# Bootstrapping parameters
n_bootstrap <- 100  # Number of bootstrap samples 
set.seed(123)      # For reproducibility

################################################################################
######################### LOAD LIBRARIES ########################################
################################################################################
if (!requireNamespace("tictoc", quietly = TRUE)) {
  install.packages("tictoc")
}
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

# Initialize a list to keep track of any failures
failures <- list()

# Loop through each segment, type, and prefix
for (prefix in prefixes) {   
  for (segment in segment_list) {
    for (type in types) {  
      if (grepl("lidar", prefix) && type != "erosion") {
        next
      }
      # Wrap the entire process in tryCatch to capture high-level errors
      tryCatch({
        
        CPF_prefixes <- c("CPF sfm", "CPF lidar")
        region <- if (prefix %in% CPF_prefixes) "CPF" else "ETF"
        
        message(paste0("\n\nProcessing: | ", prefix, " | ", type, " | ", segment, "m |\n"))
        
        base_input_folder <- file.path(region, "Inputs")
        base_output_folder <- file.path(region, "Outputs")
        segment_input_folder <- file.path(base_input_folder, paste0("Segmented ", segment, "m"))
        segment_output_folder <- file.path(base_output_folder, paste0("Segmented ", segment, "m"))
        
        output_folder <- file.path(segment_output_folder, paste0(prefix, "_", type, "_logtrans"))
        bootstrap_folder <- file.path(output_folder, "bootstrap_results")
        if (!dir.exists(bootstrap_folder)) {
          dir.create(bootstrap_folder, recursive = TRUE)
        }
        
        formula_file <- file.path(output_folder, formula_file_name)
        if (!file.exists(formula_file)) {
          stop(paste("Formula file does not exist:", formula_file))
        }
        model_formula_str <- readLines(formula_file)
        model_formula <- as.formula(model_formula_str)
        response_var <- all.vars(model_formula)[1]
        message("Model formula: \n", model_formula_str)
        
        output_file <- file.path(
          bootstrap_folder, 
          paste0(response_var, " ", segment, "m_VIF-2_corr0.6_bootstrap.txt")
        )
        
        input_streams <- file.path(base_input_folder, "Streams", "streams_100k.gpkg")
        if (!dir.exists(output_folder)) {
          dir.create(output_folder, recursive = TRUE)
        }
        
        bootstrap_results_file <- file.path(bootstrap_folder, paste0(prefix, "_bootstrap_results.", type, ".csv"))
        if (file.exists(bootstrap_results_file)) {
          message(paste0("Bootstrap results file already exists: ", bootstrap_results_file))
          next
        }
        
        ################################################################################
        ########################### BOOTSTRAPPING FUNCTION ##############################
        ################################################################################
        bootstrap_model <- function(input_obs, n_boot, model_formula, multiple_ws) {
          
          message(paste0("\nFitting model for | ", prefix, " | ", type, " | ", segment, "m | with ", n_boot, " bootstrap samples..."))
          
          obs_data <- st_read(input_obs)
          message("DEBUG: Read obs_data with ", nrow(obs_data), " rows.")
          
          # We'll allow up to 1 attempt per bootstrap iteration
          max_attempts <- 10
          
          # A helper function that attempts one bootstrap iteration and returns NULL if there is an error. 
          bootstrap_iteration_once <- function(i) {
            tryCatch({
              # New bootstrap sample for each iteration
              boot_indices <- sample(1:nrow(obs_data), replace = TRUE)
              boot_data <- obs_data[boot_indices, ]
              message("DEBUG: Bootstrap iteration ", i, " - boot_data has ", nrow(boot_data), " rows.")
              
              CP_streams <- st_read(input_streams)
              message("DEBUG: Read CP_streams with ", nrow(CP_streams), " features.")
              
              ssn_path <- file.path(output_folder, "ssn_boot", paste0(prefix, "_", type, "_logtrans_boot_", i, ".ssn"))
              lsn_out <- file.path(output_folder, "lsn_out_boot", paste0("lsn_out_", i))
              
              # Create folders if they don't exist
              if (!dir.exists(dirname(ssn_path))) {
                dir.create(dirname(ssn_path), recursive = TRUE)
              }
              if (!dir.exists(dirname(lsn_out))) {
                dir.create(dirname(lsn_out), recursive = TRUE)
              }
              
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
              
              obs <- sites_to_lsn(
                sites = boot_data,
                edges = edges,
                lsn_path = lsn_out,
                file_name = "obs",
                snap_tolerance = 10,
                save_local = TRUE,
                overwrite = TRUE
              )
              
              edges <- updist_edges(edges = edges, lsn_path = lsn_out, calc_length = TRUE)
              
              site.list <- updist_sites(
                sites = list(obs = obs), 
                edges = edges,
                length_col = "Length", 
                lsn_path = lsn_out
              )
              
              infl_col <- "flow_accum_max"
              segpi <- "flow_accum_PI"
              afv_col <- "afv_flow_accum"
              
              edges <- afv_edges(
                edges = edges, 
                infl_col = infl_col, 
                segpi = segpi,
                afv_col = afv_col, 
                lsn_path = lsn_out
              )
              site.list <- afv_sites(
                sites = site.list, 
                edges = edges, 
                afv_col = afv_col,
                save_local = TRUE, 
                lsn_path = lsn_out
              )
              
              boot_ssn <- ssn_assemble(
                edges = edges,
                lsn_path = lsn_out,
                obs_sites = site.list$obs,
                ssn_path = ssn_path,
                import = TRUE,
                overwrite = TRUE
              )
              message("DEBUG: Assembled SSN object.")
              
              ssn_create_distmat(boot_ssn)
              
              # DEBUG: Check the SSN object structure (if possible)
              message("DEBUG: SSN object class: ", class(boot_ssn))
              
              # Fit the model
              if (type == "net") {
                ssn_mod_boot <- ssn_lm(
                  formula = model_formula,
                  ssn.object = boot_ssn,
                  tailup_type = "exponential",
                  taildown_type = "none",
                  estmethod = "ml",
                  euclid_type = "spherical",
                  nugget_type = "nugget",
                  additive = "afv_flow_accum",
                  random = ~ as.factor(ch_watershed)
                )
                model_type <- "ssn_lm"
              } else {
                ssn_mod_boot <- ssn_glm(
                  formula = model_formula,
                  ssn.object = boot_ssn,
                  family = "Gamma",
                  tailup_type = "exponential",
                  taildown_type = "none",
                  estmethod = "ml",
                  euclid_type = "spherical",
                  nugget_type = "nugget",
                  additive = "afv_flow_accum",
                  random = ~ as.factor(ch_watershed)
                )
                model_type <- "ssn_glm"
              }
              message("DEBUG: Model fitting complete. Class of model: ", class(ssn_mod_boot))
              
              # Extract outputs and add iteration ID
              tidy_mod <- tidy(ssn_mod_boot, conf.int = TRUE)
              tidy_mod$bootstrap_rep <- i
              message("DEBUG: Tidy model output columns: ", paste(names(tidy_mod), collapse=", "))
              
              varcomp_mod <- varcomp(ssn_mod_boot)
              varcomp_mod$bootstrap_rep <- i
              
              loocv_mod <- loocv(ssn_mod_boot, cv_predict = TRUE, se.fit = TRUE)
              loocv_mod$bootstrap_rep <- i
              
              glance_mod <- glance(ssn_mod_boot)
              glance_mod$bootstrap_rep <- i
              
              residuals <- residuals(ssn_mod_boot)
              fitted_values <- fitted(ssn_mod_boot)
              
              res_fit_df <- data.frame(
                residuals = residuals,
                fitted_values = fitted_values,
                bootstrap_rep = i
              )
              
              # Return all outputs; note: boot_data_df is replaced by boot_data here
              return(list(
                tidy_mod = tidy_mod,
                varcomp_mod = varcomp_mod,
                loocv_mod = loocv_mod,
                glance_mod = glance_mod,
                res_fit_df = res_fit_df,
                boot_data = boot_data,  # returning boot_data for troubleshooting
                model_type = model_type
              ))
              
            }, error = function(e) {
              message("Error in iteration ", i, ": ", e$message)
              return(NULL)
            })
          }
          
          # Core iteration function: try up to max_attempts
          bootstrap_iteration <- function(i) {
            for (attempt_num in seq_len(max_attempts)) {
              if (attempt_num > 1){
                message("\nBootstrap iteration ", i, ", attempt ", attempt_num)
              }
              result <- bootstrap_iteration_once(i)
              high_covariance <- FALSE
              
              if (!is.null(result)) {
                # Identify bootstrap repetitions with extreme variance proportions
                selected_bootstrap_reps <- result$varcomp_mod$proportion > 0.9999
                varcomp_mod <- result$varcomp_mod
                varcomp_filtered <- varcomp_mod %>%
                  filter(!varcomp %in% c("taildown_de", "euclid_de", "nugget", "tailup_de"))
                selected_covariates <- varcomp_filtered$proportion == 0
                if (any(selected_bootstrap_reps) || any(selected_covariates)) {
                  high_covariance <- TRUE
                  message("DEBUG: High covariance detected in iteration ", i)
                }
              }
              
              if (!is.null(result) && (high_covariance == FALSE)) {
                return(result)
              } else {
                message("Retrying iteration ", i, "...")
              }
            }
            message("All attempts failed for iteration ", i)
            return(result)
          }
          
          # Run all bootstrap iterations using future_map
          with_progress({
            p <- progressor(along = 1:n_boot)
            bootstrap_results <- future_map(
              1:n_boot,
              function(i) {
                p(sprintf("Bootstrap iteration %d", i))
                bootstrap_iteration(i)
              },
              .options = furrr_options(seed = TRUE)
            )
          })
          
          # Extract outputs from each iteration
          results_list <- map(bootstrap_results, "tidy_mod")
          varcomp_list <- map(bootstrap_results, "varcomp_mod")
          loocv_list <- map(bootstrap_results, "loocv_mod")
          glance_list <- map(bootstrap_results, "glance_mod")
          res_fit_df <- map(bootstrap_results, "res_fit_df")
          boot_data_list <- map(bootstrap_results, "boot_data")
          model_type <- bootstrap_results[[1]]$model_type
          
          # Combine into data frames
          bootstrap_results_df <- bind_rows(results_list)
          bootstrap_varcomp_df <- bind_rows(varcomp_list)
          bootstrap_loocv <- bind_rows(loocv_list)
          bootstrap_glance <- bind_rows(glance_list)
          bootstrap_res_fit_df <- bind_rows(res_fit_df)
          bootstrap_data_df <- bind_rows(boot_data_list)
          
          # DEBUG: Check that the tidy model output contains the "term" column
          if (!"term" %in% names(bootstrap_results_df)) {
            stop("DEBUG: 'term' column not found in bootstrap results. Available columns: ", 
                 paste(names(bootstrap_results_df), collapse=", "))
          }
          
          return(list(
            bootstrap_results = bootstrap_results_df,
            bootstrap_varcomp = bootstrap_varcomp_df,
            bootstrap_loocv = bootstrap_loocv,
            bootstrap_glance = bootstrap_glance,
            bootstrap_res_fit = bootstrap_res_fit_df,
            bootstrap_data = bootstrap_data_df,
            model_type = model_type
          ))
        }
        
        # Remove the output file if it exists
        if (file.exists(output_file)) {
          file.remove(output_file)
        }
        
        cat("\nStarting Bootstrapping...\n", file = output_file, append = TRUE)
        
        bootstrap_output <- bootstrap_model(
          input_obs = input_obs,
          n_boot = n_bootstrap,
          model_formula = model_formula,
          multiple_ws = multiple_ws
        )
        
        bootstrap_results <- bootstrap_output$bootstrap_results
        bootstrap_varcomp <- bootstrap_output$bootstrap_varcomp
        bootstrap_loocv <- bootstrap_output$bootstrap_loocv
        bootstrap_glance <- bootstrap_output$bootstrap_glance
        bootstrap_res_fit <- bootstrap_output$bootstrap_res_fit
        bootstrap_data <- bootstrap_output$bootstrap_data
        model_type <- bootstrap_output$model_type
        
        # DEBUG: Check structure of bootstrap_results before grouping
        message("DEBUG: Names of bootstrap_results columns: ", paste(names(bootstrap_results), collapse=", "))
        
        bootstrap_summary <- bootstrap_results %>%
          group_by(term) %>%
          summarise(
            estimate_mean = mean(estimate, na.rm = TRUE),
            estimate_sd   = sd(estimate, na.rm = TRUE),
            conf_low      = quantile(estimate, 0.025, na.rm = TRUE),
            conf_high     = quantile(estimate, 0.975, na.rm = TRUE),
            p_value       = t.test(estimate, mu = 0)$p.value
          )
        
        cat("\nBootstrap Summary:\n", file = output_file, append = TRUE)
        capture.output(paste0("Model Type: ", model_type), file = output_file, append = TRUE)
        capture.output(print(bootstrap_summary), file = output_file, append = TRUE)
        
        # Write outputs to CSV files
        bootstrap_results_file <- file.path(bootstrap_folder, paste0(prefix, "_bootstrap_results.", type, ".csv"))
        write.csv(bootstrap_results, bootstrap_results_file, row.names = FALSE)
        
        bootstrap_varcomp_file <- file.path(bootstrap_folder, paste0(prefix, "_varcomp_results.", type, ".csv"))
        write.csv(bootstrap_varcomp, bootstrap_varcomp_file, row.names = FALSE)
        
        bootstrap_loocv_file <- file.path(bootstrap_folder, paste0(prefix, "_loocv_results.", type, ".csv"))
        write.csv(bootstrap_loocv, bootstrap_loocv_file, row.names = FALSE)
        
        bootstrap_glance_file <- file.path(bootstrap_folder, paste0(prefix, "_glance_results.", type, ".csv"))
        write.csv(bootstrap_glance, bootstrap_glance_file, row.names = FALSE)
        
        bootstrap_res_fit_file <- file.path(bootstrap_folder, paste0(prefix, "_residuals_fitted.", type, ".csv"))
        write.csv(bootstrap_res_fit, bootstrap_res_fit_file, row.names = FALSE)
        
        bootstrap_data_file <- file.path(bootstrap_folder, paste0(prefix, "_bootstrap_data.", type, ".csv"))
        write.csv(bootstrap_data, bootstrap_data_file, row.names = FALSE)
        
        #----------------------------#
        #   Residuals vs Fitted Plot #
        #----------------------------#
        fitted_values <- bootstrap_res_fit$fitted_values
        residuals <- bootstrap_res_fit$residuals
        
        resid_fitted_plot <- ggplot(data.frame(Fitted = fitted_values, Residuals = residuals), 
                                    aes(x = Fitted, y = Residuals)) +
          geom_point(color = "blue") +
          geom_hline(yintercept = 0, color = "red") +
          geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "green") +
          labs(title = "Residuals vs Fitted Values",
               x = "Fitted Values",
               y = "Residuals") +
          theme_minimal()
        
        ggsave(filename = file.path(bootstrap_folder, "Residuals_vs_Fitted.png"), 
               plot = resid_fitted_plot, width = 8, height = 6)
        
        #----------------------------#
        #         Q-Q Plot           #
        #----------------------------#
        qq_plot_gg <- ggplot(data.frame(Residuals = residuals), aes(sample = Residuals)) +
          stat_qq(color = "blue") +
          stat_qq_line(color = "red") +
          labs(title = "Q-Q Plot of Residuals") +
          theme_minimal()
        
        ggsave(filename = file.path(bootstrap_folder, "QQ_Plot_Residuals.png"), 
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
        
        ggsave(filename = file.path(bootstrap_folder, "Histogram_Residuals.png"), 
               plot = hist_gg, width = 8, height = 6)
        
        cat("\nBootstrapping Completed.\n", file = output_file, append = TRUE)
        
      }, error = function(e) {
        message(paste0(
          "Error for prefix: ", prefix, 
          " | type: ", type, 
          " | segment: ", segment, 
          "\nMessage: ", e$message
        ))
        failures[[length(failures) + 1]] <<- list(
          prefix = prefix, 
          type = type, 
          segment = segment, 
          error_message = e$message
        )
      })
      
    }
  }
}

total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", total_time$toc - total_time$tic, "seconds\n")

if (length(failures) > 0) {
  message("\nThe following parameter combinations failed:")
  print(failures)
} else {
  message("\nNo parameter combinations failed.")
}
