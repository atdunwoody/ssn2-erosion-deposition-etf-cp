library(dplyr)
library(ggplot2)
library(readr)

plot_bootstrap_results <- function(bootstrap_csv_path, 
                                   output_file = NULL,
                                   plot_width  = 8, 
                                   plot_height = 4) {
  
  # Read the bootstrap results
  df <- read_csv(bootstrap_csv_path)
  
  # Compute the mean per term (to reorder the terms by mean)
  df_summary <- df %>%
    group_by(term) %>%
    summarize(mean_est = mean(estimate), .groups = "drop")
  
  # Reorder 'term' based on ascending mean estimate
  df <- df %>%
    mutate(term = factor(
      term, 
      levels = df_summary$term[order(df_summary$mean_est)]
    ))
  
  # Create the plot using a standard boxplot
  # whis = c(0.05, 0.95) means whiskers extend to 5th and 95th percentile
  # outlier.shape = NA hides any points beyond the whiskers
  p <- ggplot(df, aes(x = term, y = estimate)) +
    geom_boxplot(
      whis          = c(0.05, 0.95),
      outlier.shape = NA,
      fill          = "gray70",
      alpha         = 0.8
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Bootstrap Summary Plot",
      x     = "Term",
      y     = "Estimate"
    ) +
    # Flip coordinates to make box plots horizontal
    coord_flip()
  
  # Print the plot to the current device
  print(p)
  
  # Optionally save to file with specified width/height
  if (!is.null(output_file)) {
    ggsave(filename = output_file, plot = p, 
           width = plot_width, height = plot_height)
  }
  
  # Return the plot object (in case the user wants to modify or save it differently)
  invisible(p)
}

results_path <- "C:\\Users\\alextd\\Documents\\GitHub\\ssn2-erosion-deposition-etf-cp\\ETF\\Outputs\\Segmented 20m\\ET lidar_erosion_logtrans\\bootstrap_results\\ET lidar_bootstrap_results.erosion.csv"

plot_bootstrap_results(results_path, 
#                       output_file = "bootstrap_summary_plot.png",
                        plot_width  = 20,  # customize as needed
                       plot_height = 6)


