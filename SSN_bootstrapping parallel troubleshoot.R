################################################################################
########################### USER DEFINED VARIABLES #############################
################################################################################

# ETF prefixes: "ET", "LM2", "LPM", "MM_ET"
# Bennett prefixes: "Bennett", "ME", "MM", "MW", "UE", "UW", "UM"
prefixes <- c("ET_test")  # Can define single or multiple prefixes

# Types: "erosion", "deposition", "net"
types <- c("erosion")  # Can define single or both erosion, deposition

# Model formula is stored in outputs folder:
#"ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt"
formula_file_name <- "ssn_formula.txt"

# SSN object already exists for each watershed/combined area, keep TRUE
load_ssn <- TRUE

# Bootstrapping parameters
n_bootstrap <- 10  # Number of bootstrap samples (adjust as needed)
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
library(classInt)
library(knitr)
library(broom)
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

for (type in types) {
  for (prefix in prefixes) {
    
    bennett_prefixes <- c("Bennett", "ME", "MM", "MW", "UE", "UW", "UM", "UM_test", "Bennett_test")
    if (prefix %in% bennett_prefixes) {
      region <- "Bennett"
    } else {
      region <- "ETF"
    }
    
    if (prefix =="MM_ET") {
      prefix <- "MM"
    }
    
    print(paste0("Processing: ", prefix, " - ", type))
    base_input_folder <- file.path(region,"Inputs")
    base_output_folder <- file.path(region,"Outputs")
    
    # Determines whether random effect of watershed is included
    if (prefix %in% c("Bennett", "ET", "ET_test")) {
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
    
    output_folder <- file.path(base_output_folder, paste0(prefix, "_", type, "_logtrans"))
    output_file <- file.path(output_folder, paste0(prefix, "_ch_sfm.", type, ".norm_VIF-2_corr0.6.txt"))
    formula_file <- file.path(output_folder, formula_file_name)
    ssn_path <- file.path(output_folder, paste0(prefix, "_", type, "_logtrans.ssn"))
    lsn_out <- file.path(output_folder, "lsn_out")
    input_streams <- file.path(base_input_folder, "Streams", "streams_100k.gpkg")
    
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    if (!file.exists(formula_file)) {
      stop(paste("Formula file does not exist:", formula_file))
    }
    
    model_formula_str <- readLines(formula_file)
    model_formula <- as.formula(model_formula_str)
    
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
      site.list <- updist_sites(sites = list(obs = obs), edges = edges,
                                length_col = "Length", lsn_path = lsn_out)
      
      infl_col <- "flow_accum_max"
      segpi <- "flow_accum_PI"
      afv_col <- "afv_flow_accum"
      
      edges <- afv_edges(edges = edges, infl_col = infl_col, segpi = segpi,
                         afv_col = afv_col, lsn_path = lsn_out)
      site.list <- afv_sites(sites = site.list, edges = edges, afv_col = afv_col,
                             save_local = TRUE, lsn_path = lsn_out)
      
      CP_ssn <- ssn_assemble(
        edges = edges,
        lsn_path = lsn_out,
        obs_sites = site.list$obs,
        ssn_path = ssn_path,
        import = TRUE,
        overwrite = TRUE
      )
    } else {
      # Load existing SSN
      CP_ssn <- ssn_import(ssn_path, overwrite = TRUE)
    }
    
    ssn_create_distmat(CP_ssn)
    
    # Extract obs_data and write to CSV for inspection
    obs_data <- ssn_get_data(CP_ssn, name = "obs")
    obs_output_file <- file.path(output_folder, paste0(prefix, "_obs_data.", type, ".csv"))
    obs_data_df <- st_drop_geometry(obs_data)
    write.csv(obs_data_df, obs_output_file, row.names = FALSE)
    
    ################################################################################
    ########################### DEBUGGING HELPERS ###################################
    ################################################################################
    # Quick function to check for near-zero variance
    near_zero_var_cols <- function(df, freq_cut = 95/5, unique_cut = 10) {
      # This is a simple approach: ratio freqCut, etc. You can also use caret::nearZeroVar.
      n <- nrow(df)
      nzv_vars <- c()
      for (col in colnames(df)) {
        vals <- df[[col]]
        # Only numeric
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
    
    # Run some checks on obs_data
    message("\n--- TROUBLESHOOTING CHECKS ---")
    
    # 1. Check if response variable has zero/negative values (for Gamma)
    response_var <- all.vars(model_formula)[1]
    if (response_var %in% colnames(obs_data_df)) {
      zero_neg_count <- sum(obs_data_df[[response_var]] <= 0, na.rm = TRUE)
      message("Count of zero or negative ", response_var, ": ", zero_neg_count)
    }
    
    # 2. Check for missing or infinite values in predictor columns
    pred_vars <- all.vars(model_formula)[-1]
    for (v in pred_vars) {
      if (v %in% colnames(obs_data_df)) {
        na_count <- sum(is.na(obs_data_df[[v]]))
        inf_count <- sum(is.infinite(obs_data_df[[v]]))
        message("NA count for ", v, ": ", na_count, " | Inf count: ", inf_count)
      }
    }
    
    # 3. Check near-zero variance columns
    nzv_cols <- near_zero_var_cols(obs_data_df)
    if (length(nzv_cols) > 0) {
      message("Near-zero variance columns identified: ", paste(nzv_cols, collapse = ", "))
    }
    
    # 4. Check for duplicates in locID or netgeom
    if ("locID" %in% names(obs_data_df)) {
      dup_locID <- sum(duplicated(obs_data_df$locID))
      message("Number of duplicate locID: ", dup_locID)
    }
    if ("netgeom" %in% names(obs_data_df)) {
      dup_netgeom <- sum(duplicated(obs_data_df$netgeom))
      message("Number of duplicate netgeom: ", dup_netgeom)
    }
    
    # 5. Print a summary of the numeric columns for quick reference
    numeric_cols <- obs_data_df %>% select(where(is.numeric))
    message("--- Data summary (numeric columns) ---")
    print(summary(numeric_cols))
    
    ################################################################################
    ########################### BOOTSTRAPPING FUNCTION ##############################
    ################################################################################
    bootstrap_model <- function(ssn_obj, n_boot, model_formula, multiple_ws) {
      
      message("Fitting model with ", n_boot, " bootstrap samples...")
      
      obs_data <- ssn_get_data(ssn_obj, name = "obs")
      
      bootstrap_iteration <- function(i) {
        boot_indices <- sample(1:nrow(obs_data), replace = TRUE)
        boot_data <- obs_data[boot_indices, ]
        
        # Assign a unique ID to each observation in the bootstrap sample
        boot_data <- boot_data %>%
          mutate(
            unique_id = i * 1e6 + row_number(),
            pid = unique_id
          ) %>%
          st_as_sf()
        
        # Update 'netgeom' to reflect the new uid while keeping the tail number
        boot_data$netgeom <- mapply(
          function(geom, uid) {
            sub("(.*\\s)\\d+\\s(\\d+)\\)$", paste0("\\1", uid, " \\2)"), geom)
          }, 
          boot_data$netgeom, 
          boot_data$unique_id
        )
        
        boot_data <- boot_data %>% select(-unique_id)
        rownames(boot_data) <- NULL
        
        # Write out the boot_data (without geometry) to CSV
        boot_data_df <- st_drop_geometry(boot_data) %>%
          mutate(bootstrap_rep = i)
        if (!dir.exists(file.path(output_folder, "boot_data"))) {
          dir.create(file.path(output_folder, "boot_data"), recursive = TRUE)
        } 
        
        write.csv(boot_data_df, file.path(output_folder, "boot_data", paste0("boot_data_", i, ".csv")))
        
        ssn_updated <- ssn_put_data(boot_data, ssn_obj, name = "obs", resize_data = TRUE)
        ssn_create_distmat(ssn_updated)
        
        # Print shape of data to check if we have enough rows
        message("Bootstrap iteration ", i, " data size: ", nrow(boot_data_df))
        
        # Fit model
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
        } else if (multiple_ws) {
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
        } else if (type == "net") {
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
          model_type <- "ssn_glm"
        }
        
        # Return results
        tidy_mod <- tidy(ssn_mod_boot, conf.int = TRUE)
        tidy_mod$bootstrap_rep <- i
        
        varcomp_mod <- varcomp(ssn_mod_boot)
        varcomp_mod$bootstrap_rep <- i
        
        list(tidy_mod = tidy_mod, varcomp_mod = varcomp_mod, boot_data = boot_data_df, model_type = model_type)
      }
      
      with_progress({
        p <- progressor(along = 1:n_boot)
        bootstrap_results <- future_map(1:n_boot, function(i) {
          p(sprintf("Bootstrap iteration %d", i))
          bootstrap_iteration(i)
        }, .options = furrr_options(seed = TRUE))
      })
      
      results_list <- map(bootstrap_results, "tidy_mod")
      varcomp_list <- map(bootstrap_results, "varcomp_mod")
      boot_data_list <- map(bootstrap_results, "boot_data")
      model_type <- bootstrap_results[[1]]$model_type
      
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
    
    ################################################################################
    ########################### BOOTSTRAPPING PROCESS ##############################
    ################################################################################
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
    
    # Write to CSV files
    bootstrap_results_file <- file.path(output_folder, paste0(prefix, "_bootstrap_results.", type, ".csv"))
    write.csv(bootstrap_results, bootstrap_results_file, row.names = FALSE)
    
    bootstrap_varcomp_file <- file.path(output_folder, paste0(prefix, "_varcomp_results.", type, ".csv"))
    write.csv(bootstrap_varcomp, bootstrap_varcomp_file, row.names = FALSE)
    
    bootstrap_data_file <- file.path(output_folder, paste0(prefix, "_bootstrap_data.", type, ".csv"))
    write.csv(bootstrap_data, bootstrap_data_file, row.names = FALSE)
    
    cat("\nBootstrapping Completed.\n", file = output_file, append = TRUE)
    
  }
}

total_time <- toc(log = TRUE, quiet = TRUE)
cat("\nTotal Script Execution Time:", total_time$toc - total_time$tic, "seconds\n")
