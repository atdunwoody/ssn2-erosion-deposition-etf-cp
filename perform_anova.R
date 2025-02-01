# Install and load necessary packages

library(dplyr)
library(ggplot2)
library(car)
library(emmeans)


# Display all column names
print(names(df))


csv_path = "C:\\Users\\alextd\\Documents\\GitHub\\ssn2-erosion-deposition-etf-cp\\Bootstrap Summary\\aggregated_bootstrap_results.csv"

# Read the data
df <- read.csv(csv_path)

# Convert relevant columns to factors
df$Region <- as.factor(df$Region)
df$`Segmentation Interval` <- as.factor(df$`segmentation.Interval`)
df$`Dependent Variable` <- as.factor(df$`Dependent.Variable`)

# Define the second factor (e.g., Fire)
# Replace this with your actual second factor if different
df$Fire <- sample(c("Yes", "No"), size = nrow(df), replace = TRUE)
df$Fire <- as.factor(df$Fire)

# Perform two-way ANOVA
anova_model <- aov(estimate ~ Region * Fire, data = df)
summary(anova_model)

# Check ANOVA assumptions
par(mfrow = c(2, 2))
plot(anova_model)

# Perform post hoc Tukey tests using emmeans
emm <- emmeans(anova_model, ~ Region * Fire)
pairwise_comparisons <- contrast(emm, method = "pairwise", adjust = "tukey")
summary(pairwise_comparisons)

# Alternatively, using TukeyHSD
tukey_results <- TukeyHSD(anova_model, which = "Region:Fire")
print(tukey_results)

# Visualize interaction effects
ggplot(df, aes(x = Region, y = estimate, fill = Fire)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Two-Way ANOVA: Region and Fire Effects",
       y = "Estimate",
       x = "Region") +
  theme_minimal()
