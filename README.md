# ssn2-erosion-deposition-etf-cp

## Overview
This script is used to process multiple SSN2 (Spatial Stream Network) models for different watershed areas in the East-Troublesome (ET) and Cameron Peak (Bennett) burn scars. The script processes multiple watershed regions (prefixes) with specified response variables (erosion or deposition). It includes SSN2 object creation/loading, model fitting using generalized linear models (GLM), the use of random effects for combined watershed models and model diagnostics.

## User Inputs into SSN_run_multiple.R
The user needs to define the following parameters at the start of the script:

### 1. Watershed Prefixes
```r
prefixes <- c("LM2", "LPM")
```
- These are the prefixes for the watershed regions to be processed.
- Bennett Options:
```r
prefixes <- c(Bennett", "ME", "MM", "MW", "UE", "UW", "UM")
```
- ET options
```r
prefixes <- c("ET", "LM2", "LPM", "MM_ET")
```

### 2. Erosion/Deposition variable setting
```r
types <- c("erosion", "deposition")
```
- Specifies the response variable to be fitted for each watershed.
- Options: Can be a single type ("erosion" or "deposition") or both.

### 3. Formula file 
```r
formula_file_name <- "ssn_formula.txt"
```
- Ensure that ssn_formula.txt exists in the output folder for each prefix and type specified above.
- Example Formula:
```r
  ch_sfm.erosion.norm ~ ch_curvature.median + ch_valley_width + 
    ch_change.in.slope.over.width + ch_channel.width.over.valley.width+     hs_eastness.median + hs_curvature.median + hs_ndvi.max + 
    hs_dnbr.median + hs_drainage_density
```
### 4. Load SSN
```r
load_ssn <- TRUE
```
- If TRUE, the script will load existing SSN objects.
- If FALSE, it will create new SSN objects using data in Inputs folder

## Outputs folder

### 1. Model Summary and Results

- **`<prefix>_ch_sfm.<type>.norm_VIF-2_corr0.6.txt`**
- e.g. `LM2_ch_sfm.erosion.norm_VIF-2_corr0.6.txt`

  This text file contains information about the fitted statistical model, including:

  - **Model Formula:** The formula used for the model.
  - **SSN Path:** Path to the Spatial Stream Network (SSN) object file.
  - **LSN Output Path:** Path to the Local Stream Network (LSN) output directory.
  - **Model Summary:** Detailed summary of the fitted model (`summary(ssn_mod)`).
  - **Tidy Model Outputs:** Cleaned and organized model coefficients with confidence intervals (`tidy(ssn_mod, conf.int = TRUE)`).
  - **Variance Components:** Variance components of the model (`varcomp(ssn_mod)`).
  - **Leave-One-Out Cross-Validation (LOOCV):** Leave one out cross validation results (`loocv(ssn_mod)`).
  - **Glance:** Model diagnostics and goodness-of-fit metrics (`glance(ssn_mod)`).

### 2. SSN Object

- **`<prefix>_<type>_logtrans.ssn`**
- e.g. `LM2_erosion_logtrans.ssn`

  This file stores the assembled Spatial Stream Network (SSN2) object for the specific prefix and type. It contains all spatial and attribute data necessary for modeling.

### 3. Local Stream Network (LSN) Outputs

- **`lsn_out/`**

  Directory containing Local Stream Network (LSN) processing outputs, including:

  - Edge distances
  - AFV (Area Flow Variable) calculations
  - Other intermediary spatial network data

### 4. Diagnostic Plots

- **`Residuals_vs_Fitted.png`**

  Scatter plot of residuals versus fitted values to assess the homoscedasticity and linearity of the model.

- **`QQ_Plot_Residuals_ggplot2.png`**

  Q-Q plot of residuals to evaluate the normality assumption of the model residuals.

- **`Histogram_Residuals_ggplot2.png`**
 
  Histogram of residuals overlaid with a normal distribution curve to visualize the distribution of residuals.

### 5. Statistical Tests

- **`Statistical_Tests.txt`**
 
  This text file includes the results of statistical tests conducted on the model residuals:

  - **Shapiro-Wilk Test for Normality:** Assesses whether the residuals are normally distributed.
  - **Moran's I Test for Spatial Autocorrelation:** Evaluates the presence of spatial autocorrelation in the residuals.
  - **Summary Statement:** Provides a concise interpretation of the test results, indicating whether the model passed or failed each test based on a significance level of 0.05.

## Additional Notes

- **Random Effects:**  
  If multiple watersheds are involved (`multiple_ws = TRUE`), the model includes a random effect for the watershed factor (`random = ~ as.factor(ch_watershed)`).

- **Model Family and Link Function:**  
  The models are fitted using the Gamma family with an exponential tail-up autocorrelation structure and guassian euclidean autocorrelation structure. This was selected based on a sensitivty test with all iterations of autocorrelation structures for tailup, taildown, and euclidean variance types, using AIC and RMSPE to evaluate goodness of model fit.  

- **LSN Snap Tolerance:**  
  The script uses a snap tolerance of 1 meter (if stream input data is projected) for stream snapping during LSN creation. This can be adjusted in the script if input points are further away from the stream network.

## Troubleshooting

- **Missing Formula File:**  
  If the script stops with an error about a missing formula file (e.g., `ssn_formula.txt`), ensure that the formula file exists in the specified output folder  (e.g.,`ETF/Outputs/LM2_erosion_logtrans/ssn_formula.txt`).

### Contact
For questions or issues, please contact [Alex Thornton-Dunwoody](mailto:alextd@colostate.edu) 
