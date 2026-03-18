# Set seed for reproducibility
set.seed(123)

# Number of p-values to simulate
n <- 10000  # Larger number of p-values for better coverage

# Simulate p-values from a chi-squared distribution
p_values <- pchisq(rchisq(n, df=1), df=1)

# Parameters for inflation and mixing
inflation_factor <- 50  # Increase this factor to significantly inflate p-values
mixing_proportion <- 0.1  # Reduce this proportion to decrease the influence of the uniform distribution

# Apply a significant inflation factor
inflated_p_values <- p_values * inflation_factor

# Ensure the p-values remain between 0 and 1
inflated_p_values <- pmin(inflated_p_values, 1)

# Introduce variability by combining with a uniform distribution
uniform_p_values <- runif(n, min = 0, max = 1)
mixed_p_values <- (inflated_p_values * (1 - mixing_proportion) + uniform_p_values * mixing_proportion)

# Ensure the p-values remain between 0 and 1
mixed_p_values <- pmin(mixed_p_values, 1)

# Create Q-Q plot function
qqplot <- function(p_values, main = "Q-Q plot of Simulated p-values", xlim = c(0, 5), ylim = c(0, 8), pch = 19, col = "blue4", cex = 1, las = 1) {
  # Compute the theoretical quantiles
  theo_quantiles <- -log10(ppoints(length(p_values)))
  
  # Compute the observed quantiles
  obs_quantiles <- -log10(sort(p_values))
  
  # Plot the Q-Q plot
  plot(theo_quantiles, obs_quantiles, main = main, xlim = xlim, ylim = ylim, pch = pch, col = col, cex = cex, las = las, xlab = "Theoretical -log10(p)", ylab = "Observed -log10(p)")
  
  # Add reference line
  abline(a = 0, b = 1, col = "red")
}

# Generate the Q-Q plot with adjusted inflation factor and mixing proportion
qqplot(mixed_p_values)