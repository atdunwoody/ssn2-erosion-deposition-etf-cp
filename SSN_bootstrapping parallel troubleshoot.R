################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET sfm", "LM2 sfm", "LPM sfm", "MM_ET"
#               "ET lidar", "LM2 lidar", "LPM lidar", "MM_ET lidar"
# Bennett prefixes: "Bennett sfm", "ME sfm", "MM sfm", "MW sfm", "UE sfm", "UW sfm", "UM sfm"
#                  "Bennett lidar", "ME lidar", "MM lidar", "MW lidar", "UE lidar", "UW lidar", "UM lidar"
prefixes <- c(
  "ET sfm",
  "ET lidar",
  "Bennett lidar",
  "Bennett sfm"
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
n_bootstrap <- 10  # Number of bootstrap samples 
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

# Initialize a list to keep track of any failures
failures <- list()

# Loop through each segment, type, and prefix
for (prefix in prefixes) {  
  for (segment in segment_list) {
    for (type in types) {
    
      # Wrap the entire process in tryCatch
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
        
        if (prefix =="MM_ET sfm") {
          prefix <- "MM sfm"
        } else if (prefix == "MM_ET lidar") {
          prefix <- "MM lidar"
        }
        
        print(paste0("Processing: ", prefix, " - ", type, " with ", segment, "m spacing"))
        
        base_input_folder <- file.path(region,"Inputs")
        base_output_folder <- file.path(region,"Outputs")
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
        
        output_file <- file.path(
          output_folder, 
          paste0(prefix, "_ch_sfm.", type, ".norm_VIF-2_corr0.6.txt")
        )
        
        formula_file <- file.path(
          output_folder, 
          formula_file_name
        )
        
        if (!file.exists(formula_file)) {
          stop(paste("Formula file does not exist:", formula_file))
        }
        model_formula_str <- readLines(formula_file)
        model_formula <- as.formula(model_formula_str)
        
        print(model_formula_str)
        
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
        
        ################################################################################
        ########################### SSN2 PREPROCESSING ##################################
        ################################################################################
        if (!load_ssn) {
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
          
          obs <- sites_to_lsn(
            sites = CP_obs,
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
          
          CP_ssn <- ssn_assemble(
            edges = edges,
            lsn_path = lsn_out,
            obs_sites = site.list$obs,
            ssn_path = ssn_path,
            import = TRUE,
            overwrite = TRUE
          )
        } else {
          CP_ssn <- ssn_import(ssn_path, overwrite = TRUE)
        }
        
        ssn_create_distmat(CP_ssn)
        
        obs_data <- ssn_get_data(CP_ssn, name = "obs")
        obs_output_file <- file.path(output_folder, paste0(prefix, "_obs_data.", type, ".csv"))
        obs_data_df <- st_drop_geometry(obs_data)
        write.csv(obs_data_df, obs_output_file, row.names = FALSE)
        
        ################################################################################
        ########################### DEBUGGING HELPERS ###################################
        ################################################################################
        near_zero_var_cols <- function(df, freq_cut = 95/5, unique_cut = 10) {
          n <- nrow(df)
          nzv_vars <- c()
          for (col in colnames(df)) {
            vals <- df[[col]]
            if (is.numeric(vals)) {
              tbl <- table(vals)
              if (length(tbl) == 1) {
                nzv_vars <- c(nzv_vars, col)
              } else {
                freq_ratio <- max(tbl) / sort(tbl, decreasing = TRUE)[2]
                pct_unique <- 100 * (length(unique(vals)) / n)
                if (freq_ratio > freq_cut && pct_unique < unique_cut) {
                  nzv_vars <- c(nzv_vars, col)
                }
              }
            }
          }
          nzv_vars
        }
        
        message("\n--- TROUBLESHOOTING CHECKS ---")
        response_var <- all.vars(model_formula)[1]
        if (response_var %in% colnames(obs_data_df)) {
          zero_neg_count <- sum(obs_data_df[[response_var]] <= 0, na.rm = TRUE)
          message("Count of zero or negative ", response_var, ": ", zero_neg_count)
        }
        
        pred_vars <- all.vars(model_formula)[-1]
        for (v in pred_vars) {
          if (v %in% colnames(obs_data_df)) {
            na_count <- sum(is.na(obs_data_df[[v]]))
            inf_count <- sum(is.infinite(obs_data_df[[v]]))
            variance <- var(obs_data_df[[v]], na.rm = TRUE)
            message(v, " | NA count: ", na_count, " | Inf count: ", inf_count, "| Variance: ", variance)
          }
        }
        
        nzv_cols <- near_zero_var_cols(obs_data_df)
        if (length(nzv_cols) > 0) {
          message("Near-zero variance columns identified: ", paste(nzv_cols, collapse = ", "))
        }
        
        if ("locID" %in% names(obs_data_df)) {
          dup_locID <- sum(duplicated(obs_data_df$locID))
          message("Number of duplicate locID: ", dup_locID)
        }
        
        # Uncomment to get quartile ranges of data
        # numeric_cols <- obs_data_df %>% select(where(is.numeric))
        # message("--- Data summary (numeric columns) ---")
        # print(summary(numeric_cols))
        

          ################################################################################
          ########################### BOOTSTRAPPING FUNCTION ##############################
          ################################################################################
          bootstrap_model <- function(ssn_obj, n_boot, model_formula, multiple_ws) {
            
            message("Fitting model with ", n_boot, " bootstrap samples...")
            
            obs_data <- ssn_get_data(ssn_obj, name = "obs")
            
            # We'll allow up to 5 attempts for each bootstrap iteration
            max_attempts <- 5
            
            # A helper function that attempts to run one bootstrap iteration
            # and returns NULL if there is an error. 
            bootstrap_iteration_once <- function(i) {
              tryCatch({
                # New bootstrap sample each attempt
                boot_indices <- sample(1:nrow(obs_data), replace = TRUE)
                boot_data <- obs_data[boot_indices, ]
                
                boot_data <- boot_data %>%
                  mutate(
                    unique_id = i * 1e6 + row_number(),
                    pid = unique_id
                  ) %>%
                  st_as_sf()
                
                # Update netgeom with the newly generated unique_id
                boot_data$netgeom <- mapply(
                  function(geom, uid) {
                    sub("(.*\\s)\\d+\\s(\\d+)\\)$",
                        paste0("\\1", uid, " \\2)"), 
                        geom)
                  },
                  boot_data$netgeom,
                  boot_data$unique_id
                )
                
                boot_data <- boot_data %>% select(-unique_id)
                rownames(boot_data) <- NULL
                
                boot_data_df <- st_drop_geometry(boot_data) %>%
                  mutate(bootstrap_rep = i)
                
                # Optionally write boot_data_df to CSV
                if (!dir.exists(file.path(output_folder, "boot_data"))) {
                  dir.create(file.path(output_folder, "boot_data"), recursive = TRUE)
                }
                write.csv(
                  boot_data_df,
                  file.path(output_folder, "boot_data", paste0("boot_data_", i, ".csv")),
                  row.names = FALSE
                )
                
                # Update SSN object with these bootstrapped data
                ssn_updated <- ssn_put_data(boot_data, ssn_obj, name = "obs", resize_data = TRUE)
                ssn_create_distmat(ssn_updated)
                
                message("Bootstrap iteration ", i, " data size: ", nrow(boot_data_df))
                
                # Fit the model
                if (multiple_ws && type == "net") {
                  ssn_mod_boot <- ssn_lm(
                    formula = model_formula,
                    ssn.object = ssn_updated,
                    tailup_type = "exponential",
                    taildown_type = "none",
                    estmethod = "ml",
                    euclid_type = "spherical",
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
                    estmethod = "ml",
                    euclid_type = "spherical",
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
                    estmethod = "ml",
                    euclid_type = "spherical",
                    nugget_type = "nugget",
                    additive = "afv_flow_accum"
                  )
                  model_type <- "ssn_lm"
                } else {
                  ssn_mod_boot <- ssn_glm(
                    formula = model_formula,
                    ssn.object = ssn_updated,
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
                
                # On success, return the results
                tidy_mod <- tidy(ssn_mod_boot, conf.int = TRUE)
                tidy_mod$bootstrap_rep <- i
                
                varcomp_mod <- varcomp(ssn_mod_boot)
                varcomp_mod$bootstrap_rep <- i
                
                list(
                  tidy_mod = tidy_mod,
                  varcomp_mod = varcomp_mod,
                  boot_data = boot_data_df,
                  model_type = model_type
                )
                
              }, error = function(e) {
                # Return NULL so we know the attempt failed
                message("Error in iteration ", i, ": ", e$message)
                NULL
              })
            }
            
            # Now we define the core iteration function that attempts up to max_attempts
            # times until it succeeds or we run out of attempts.
            bootstrap_iteration <- function(i, p) {
              for (attempt_num in seq_len(max_attempts)) {
                message("\nBootstrap iteration ", i, ", attempt ", attempt_num)
                result <- bootstrap_iteration_once(i)
                high_covariance <- FALSE
                
                
                if (!is.null(result)) {
                  # Identify bootstrap repetitions with proportion > 0.99
                  selected_bootstrap_reps <- result$varcomp_mod$proportion > 0.99
                  varcomp_mod <- result$varcomp_mod
                  
                  # Filter out unwanted varcomp entries
                  varcomp_filtered <- varcomp_mod %>%
                    filter(!varcomp %in% c("taildown_de", "euclid_de", "nugget", "tailup_de"))
                  
                  selected_covariates <- varcomp_filtered$proportion == 0
                  

                  if (any(selected_bootstrap_reps) || any(selected_covariates)) {
                    high_covariance <- TRUE
                  }
                }
                
                if (!is.null(result) && (high_covariance == FALSE)) {
                  # If it succeeds, return the result immediately
                  return(result)
                } else {
                  message("Retrying iteration ", i, "...")
                }
              }
              # If we've reached here, all attempts have failed
              stop("Exceeded max attempts (", max_attempts, 
                   ") for bootstrap iteration ", i)
              p(sprintf("Bootstrap iteration %d", i))
            }
            
            # Use future_map or a loop to run all bootstrap iterations
            with_progress({
              p <- progressor(along = 1:n_boot)
              bootstrap_results <- future_map(
                1:n_boot,
                function(i, p) {
                  bootstrap_iteration(i, p())
                },
                .options = furrr_options(seed = TRUE)
              )
            })
            
            # Extract model results
            results_list <- map(bootstrap_results, "tidy_mod")
            varcomp_list <- map(bootstrap_results, "varcomp_mod")
            boot_data_list <- map(bootstrap_results, "boot_data")
            model_type <- bootstrap_results[[1]]$model_type
            
            # Combine into data frames
            bootstrap_results_df <- bind_rows(results_list)
            bootstrap_varcomp_df <- bind_rows(varcomp_list)
            bootstrap_data_df <- bind_rows(boot_data_list)
            
            list(
              bootstrap_results = bootstrap_results_df,
              bootstrap_varcomp = bootstrap_varcomp_df,
              bootstrap_data = bootstrap_data_df,
              model_type = model_type
            )
          }
          
        
        cat("\nStarting Bootstrapping...\n", file = output_file, append = TRUE)
        
        bootstrap_output <- bootstrap_model(
          ssn_obj = CP_ssn,
          n_boot = n_bootstrap,
          model_formula = model_formula,
          multiple_ws = multiple_ws
        )
        
        bootstrap_results <- bootstrap_output$bootstrap_results
        bootstrap_varcomp <- bootstrap_output$bootstrap_varcomp
        bootstrap_data <- bootstrap_output$bootstrap_data
        model_type <- bootstrap_output$model_type
        
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
        
        bootstrap_results_file <- file.path(output_folder, paste0(prefix, "_bootstrap_results.", type, ".csv"))
        write.csv(bootstrap_results, bootstrap_results_file, row.names = FALSE)
        
        bootstrap_varcomp_file <- file.path(output_folder, paste0(prefix, "_varcomp_results.", type, ".csv"))
        write.csv(bootstrap_varcomp, bootstrap_varcomp_file, row.names = FALSE)
        
        bootstrap_data_file <- file.path(output_folder, paste0(prefix, "_bootstrap_data.", type, ".csv"))
        write.csv(bootstrap_data, bootstrap_data_file, row.names = FALSE)
        
        cat("\nBootstrapping Completed.\n", file = output_file, append = TRUE)
        
      }, error = function(e) {
        message(
          paste0(
            "Error for prefix: ", prefix, 
            " | type: ", type, 
            " | segment: ", segment, 
            "\nMessage: ", e$message
          )
        )
        # Record the failed combination
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

# Print failures at the end if any
if (length(failures) > 0) {
  message("\nThe following parameter combinations failed:")
  print(failures)
} else {
  message("\nNo parameter combinations failed.")
}



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

# Load necessary library
library(dplyr)

# Calculate variance components
varcomp_mod <- varcomp(ssn_mod)

# Print the original variance components
print(varcomp_mod$varcomp)

# Remove rows where varcomp column equals "taildown_de" using base R
varcomp_filtered <- varcomp_mod$varcomp[varcomp_mod$varcomp != "taildown_de"]

# Print the filtered variance components
print(varcomp_filtered)

# Print the filtered variance components
print(varcomp_filtered)


# Print the filtered variance components
print(varcomp_filtered)


