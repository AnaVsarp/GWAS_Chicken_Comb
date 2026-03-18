install.packages("rrBLUP")
install.packages("corrplot")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("hglm")
install.packages("qqman")

library(rrBLUP)
library(corrplot)
library(ggplot2)
library(dplyr)
library(hglm)
library(qqman)
library(gridExtra)

#Source function file
source("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/R_files/GWAS_Functions.R")

#Load Input Files
Geno_file <- read.table("A:/Books/SEM 4/GWAS_1/genotypes_qc.raw", header = TRUE)
Pheno_file <- read.csv("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/GWAS_files/Pheno_file.csv", header = TRUE)
pheno_covariate <-read.csv("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/GWAS_files/pheno_covariate.csv", header = TRUE)

#Select the required columns to merge that are present in Pheno file
#pheno_covariate <- dplyr::select(Pheno_file, feed, cage.pen, breed, weight)

#Select the required columns to merge
Pheno_IID <- Pheno_file %>% dplyr::select(IID,breed,comb_mass)

#Merge the Genotype file with the phenotypefile to unify all the required data
data_file <- merge(Geno_file, Pheno_IID, by = "IID")

#Seperate the Genotype data based on the breed 
Bovans_Geno_data <- filter(data_file, breed == "Bovans")
LSL_Geno_data <- filter(data_file, breed == "LSL")

#Remove phenotypes after merger to retain only the Genotype data
Bovans_Geno_data <- Bovans_Geno_data %>% dplyr::select(-breed, -comb_mass)
LSL_Geno_data <- LSL_Geno_data %>% dplyr::select(-breed,-comb_mass)

#Seperate the Phenotype data based on the breed
Bovans_Pheno_data <- subset(Pheno_file, breed == "Bovans")
LSL_Pheno_data <- subset(Pheno_file, breed =="LSL")

#Seperate the covariate data based on breeds
pheno_covariate_Bovans <- subset(pheno_covariate, breed == "Bovans")
pheno_covariate_LSL <- subset(pheno_covariate, breed == "LSL")

#Individual Breed GWAS analysis
# Prune SNPs
geno_pruned_Bovans <- prune_snp_matrix(Bovans_Geno_data)
geno_pruned_LSL <- prune_snp_matrix(LSL_Geno_data)

# Construct GRM
G_B <- construct_grm(geno_pruned_Bovans)
G_L <- construct_grm(geno_pruned_LSL)

# Decomposition of the GRM to create an incidence matrix for the additive genetic effect
Z_grm_B <- decompose_grm(G_B, FALSE)
Z_grm_L <- decompose_grm(G_L, FALSE)

# Add an intercept column (a column of ones) to the covariate matrix
X_Bovans <-model.matrix(~ 1 + weight + feed + cage.pen , pheno_covariate_Bovans)
X_LSL <-model.matrix(~ 1 + weight + feed + cage.pen , pheno_covariate_LSL)

# Number of columns for each random effect; here only one
RandC_B <- ncol(Z_grm_B)
RandC_L <- ncol(Z_grm_L)

# Run GWAS for individual breeds
gwas_B <- run_gwasA(pheno = Bovans_Pheno_data$comb_mass, X = X_Bovans, Z = Z_grm_B, RandC = RandC_B, geno_pruned_Bovans)
gwas_L <- run_gwasA(pheno = LSL_Pheno_data$comb_mass, X = X_LSL, Z = Z_grm_L, RandC = RandC_L, geno_pruned_LSL)

# Placeholder for marker positions
marker_pos <- read.table("A:/Books/SEM 4/Thesis/GWAS/Script/Work_data/Input_files/genotypes.map")
colnames(marker_pos) <- c("chromosome", "marker", "distance", "position")

#Removing extra characters to make markerID uniform in all files
gwas_B$marker_id <- sub("_[ATGC]$", replacement = "", gwas_B$marker_id)
gwas_L$marker_id <- sub("_[ATGC]$", replacement = "", gwas_L$marker_id)

#Join the output to the markerID file
gwas_pos_B <- inner_join(gwas_B, marker_pos, by = c("marker_id" = "marker"))
gwas_pos_L <- inner_join(gwas_L, marker_pos, by = c("marker_id" = "marker"))

# Look at significant positions in the two scans
sig_gwas_pos_B <- filter(gwas_pos_B, p < 1e-4)
sig_gwas_pos_L <- filter(gwas_pos_L, p < 1e-4)

#To highlight sig_SNP in Man. plot
sig_pos_B <- sig_gwas_pos_B$marker_id
sig_pos_L <- sig_gwas_pos_L$marker_id

# Prepare the data for qqman
gwas_results_B <- data.frame(SNP = gwas_pos_B$marker_id, CHR = gwas_pos_B$chromosome, BP = gwas_pos_B$position, P = gwas_pos_B$p)
gwas_results_L <- data.frame(SNP = gwas_pos_L$marker_id, CHR = gwas_pos_L$chromosome, BP = gwas_pos_L$position, P = gwas_pos_L$p)

gwas_results_B$CHR <- as.numeric(gwas_results_B$CHR)
gwas_results_L$CHR <- as.numeric(gwas_results_L$CHR)

gwas_results_B <- gwas_results_B %>%filter(!is.na(CHR) & !is.na(BP) & !is.na(P))
gwas_results_L <- gwas_results_L %>%filter(!is.na(CHR) & !is.na(BP) & !is.na(P))

# Generate Manhattan plot
manhattan(gwas_results_B,main = "Manhattan Plot: Bovans", highlight = sig_pos_B, col = c("blue4", "orange3"), annotatePval = 0.0001, genomewideline = -log10(5e-8), suggestiveline = -log10(1e-4),chrlabs = c(1:22, "X", "Y"))
manhattan(gwas_results_L,main = "Manhattan Plot: LSL", highlight = sig_pos_L, annotatePval = 0.0001, genomewideline = -log10(5e-8), suggestiveline = -log10(1e-4),chrlabs = c(1:22, "X", "Y"))

# Create Manhattan plot for gwasA_pos (Bovans)
manhattan_gwas_B <- ggplot(gwas_pos_B, aes(x = as.numeric(factor(chromosome)), y = -log10(p))) +
  geom_point(aes(color = as.factor(chromosome)), size = 1.5) +  # Color by chromosome
  labs(title = "Manhattan Plot: Bovans", x = "Chromosome", y = "-log10(p)") +
  geom_hline(yintercept = -log10(1e-4), linetype = "dotted", color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot for Bovans
print(manhattan_gwas_B)

# Create Manhattan plot for gwasA_pos (LSL)
manhattan_gwas_L <- ggplot(gwas_pos_L, aes(x = as.numeric(factor(chromosome)), y = -log10(p))) +
  geom_point(aes(color = as.factor(chromosome)), size = 1.5) +  # Color by chromosome
  labs(title = "Manhattan Plot: LSL", x = "Chromosome", y = "-log10(p)") +
  geom_hline(yintercept = -log10(1e-4), linetype = "dotted", color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot for LSL
print(manhattan_gwas_L)

# Display both plots side by side
grid.arrange(manhattan_gwas_B, manhattan_gwas_L, ncol = 2)

#Create qqplot
qq(gwas_results_B$P, main = "Q-Q plot of GWAS: Bovans p-values", xlim = c(0, 6), ylim = c(0,8), pch = 19, col = "blue4", cex = 1, las = 1)
qq(gwas_results_L$P, main = "Q-Q plot of GWAS: LSL p-values", xlim = c(0, 6), ylim = c(0,8), pch = 19, col = "blue4", cex = 1, las = 1)


# Step 1: Calculate chi-squared statistics
gwas_results_B$chi_sq <- qchisq(gwas_results_B$P, df=1, lower.tail=FALSE)
gwas_results_L$chi_sq <- qchisq(gwas_results_L$P, df=1, lower.tail=FALSE)

# Step 2: Calculate the observed and expected median of chi-squared statistics
observed_median_B <- median(gwas_results_B$chi_sq)
expected_median_B <- qchisq(0.5, df=1) # This is 0.456 for chi-squared with 1 degree of freedom

observed_median_L <- median(gwas_results_L$chi_sq)
expected_median_L <- qchisq(0.5, df=1) # This is 0.456 for chi-squared with 1 degree of freedom

# Step 3: Calculate lambda genomic inflation factor
lambda_B <- observed_median_B/ expected_median_B
lambda_L <- observed_median_L/ expected_median_L

# Print the lambda value
print(paste("Genomic inflation factor (lambda):", lambda_B))
print(paste("Genomic inflation factor (lambda):", lambda_L))
