# Install and load necessary packages
install.packages("rrBLUP")
install.packages("corrplot")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("hglm")

library(rrBLUP)
library(corrplot)
library(ggplot2)
library(dplyr)
library(hglm)

# Comb size vs diet
#comb_diet_count <- combined %>% count(feed, breed)

anova_weight_feed <- aov(comb_mass.x ~ feed, data = combined)
summary(anova_weight_feed)

#print(comb_diet_count)
ggplot(combined, aes(x = feed, y = comb_mass.x, color = breed)) +
  geom_jitter(width = 0.4, alpha = 1) +  # Jitter to avoid overlap
  stat_summary(fun = mean, geom = "point", size = 5, shape = 20, color = "black") +
  labs(title = "Feed vs Comb Weight", x = "feed", y = "comb weight (gm)") + 
  theme_minimal()  # Minimalist theme for better visualization

# Comb size vs housing (cage.pen)
#comb_housing_count <- combined %>% count(cage.pen, breed)

anova_comb_housing <- aov(comb_mass.x ~ cage.pen, data = combined)
summary(anova_comb_housing)

#print(comb_housing_count)
ggplot(combined, aes(x = cage.pen, y = comb_mass.x, color = breed)) +
  geom_jitter(width = 0.4, alpha = 1) +  
  stat_summary(fun = mean, geom = "point", size = 5, shape = 20, color = "black") +
  labs(title = "Housing vs Comb Weight", x = "housing", y = "comb weight (gm)") + 
  theme_minimal()  

# Comb size vs crossbreeds
#comb_breed_count <- combined %>% count(breed)

comb_breed_lm <- lm(comb_mass.x ~ breed, data = combined)
anova(comb_breed_lm)
summary(comb_breed_lm)

#print(comb_breed_count)
ggplot(combined, aes(x = breed, y = comb_mass.x, color = breed)) +
  geom_jitter(width = 0.4, alpha = 1) +  
  stat_summary(fun = mean, geom = "point", size = 5, shape = 20, color = "black") +
  labs(title = "Crossbreeds vs Comb weight", x = "breed", y = "Comb Weight (gm)") + 
  theme_minimal()  

# Comb size vs body mass
#comb_body_mass_count <- combined %>% count(weight, breed)

correlation_comb_weight <- cor.test(combined$comb_mass.x, combined$weight)
print(correlation_comb_weight)

#print(comb_body_mass_count)
# Calculate means
mean_comb_g <- mean(combined$comb_g)
mean_weight <- mean(combined$weight)
# Create scatter plot
ggplot(combined, aes(x = comb_mass.x, y = weight, color = breed)) +
  geom_point() +
  # Add mean points
  geom_point(data = data.frame(comb_g = mean_comb_g, weight = mean_weight),
             aes(x = comb_g, y = weight),
             color = "black", size = 3, shape = 17) +
  # Add mean text labels
  geom_text(data = data.frame(comb_g = mean_comb_g, weight = mean_weight),
            aes(x = comb_g, y = weight, label = paste("Mean: (", round(mean_comb_g, 2), ",", round(mean_weight, 2), ")")),
            vjust = 9, hjust = -3.75, color = "red", size = 3) +
  # Add linear regression line
  geom_smooth(method = "lm", se = FALSE, aes(colour = breed)) +
  labs(title = "Comb weight vs Body mass",x = "Comb weight (gm)", y = "Body Mass (kg)") +
  theme_minimal()

# Comb size vs load_N (bone breaking strength)
# comb_loadN_count <- combined %>% count(load_N, breed)
#print(comb_loadN_count)
# Calculate means
mean_comb_g <- mean(combined$comb_g, na.rm = TRUE)
mean_load_N <- mean(combined$load_N, na.rm = TRUE)

# Create scatter plot
ggplot(combined, aes(x = comb_g, y = load_N, color = breed)) +
  geom_jitter(width = 0.4, alpha = 1) +  
  # Add mean points
  geom_point(data = data.frame(comb_g = mean_comb_g, load_N = mean_load_N),
             aes(x = comb_g, y = load_N),
             color = "black", size = 3, shape = 17) +
  # Add mean text labels
  geom_text(data = data.frame(comb_g = mean_comb_g, load_N = mean_load_N),
            aes(x = comb_g, y = load_N, label = paste("Mean: (", round(mean_comb_g, 2), ",", round(mean_load_N, 2), ")")),
            vjust = 9, hjust = -3, color = "red", size = 3) +
  # Add linear regression lines for each breed
  geom_smooth(method = "lm", se = FALSE, aes(color = breed)) +
  labs(title = "Comb weight vs Bone breaking strength", x = "Comb Weight (gm)", y = "Bone breaking strength (N)") + 
  theme_minimal() n


# Plot histogram of chicken weight
hist({pheno}P$weight, col = "black", xlab = 'Weight', ylab = 'Frequency', border = 'white', breaks = 10, main = 'Weight Histogram')

# Plot histogram of chicken comb weight
hist(comb_g$comb_g, col = "blue", xlab = 'Comb Weight', ylab = 'Frequency', border = 'white', breaks = 10, main = 'Comb Weight Histogram')

qqnorm(cleaned_data)
qqline(cleaned_data, col = "red")  # Adds a reference line

# Boxplots
boxplot(pheno$weight, xlab = 'Boxplot', col = 'red', ylab = 'weight', ylim = c(1.6, 2.7), main = 'Boxplot of Weight')
boxplot(comb_g$comb_g, xlab = 'Boxplot', col = 'pink', ylab = 'comb_weight(gm)', ylim = c(15.00, 30.00), main = 'Boxplot of Weight')

# Plot histograms

ggplot(Pheno_file, aes(x = breed, fill = cage.pen)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  facet_wrap(~ cage.pen) +
  labs(title = "Phenotype Distribution in Different Housing Systems",
       x = "Breed",
       y = "Count") +
  theme_minimal()

ggplot(Pheno_file, aes(x = comb_mass, fill = cage.pen)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  facet_wrap(~ cage.pen) +
  labs(title = "Phenotype Distribution in Different Housing Systems",
       x = "combmass",
       y = "Frequency") +
  theme_minimal()

# Load PCA results
pca <- read.table("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/PLINK/PCA/genotype_qc_pca.eigenvec")
combine_pca <- inner_join(pheno, pca, by =c("animal_id" = "V1" ))

# Create the second plot
plot_pca <- ggplot(combine_pca, aes(x = V3, y = V4, colour = breed)) + 
  geom_point(width = 0.4, alpha = 1) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()
print(plot_pca)





