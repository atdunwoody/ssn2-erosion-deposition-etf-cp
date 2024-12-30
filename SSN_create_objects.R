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

  }
}
    