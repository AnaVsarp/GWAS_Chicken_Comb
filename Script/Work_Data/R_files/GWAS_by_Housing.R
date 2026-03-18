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
Pheno_IID <- Pheno_file %>% dplyr::select(IID,cage.pen,comb_mass)

#Merge the Genotype file with the phenotypefile to unify all the required data
data_file <- merge(Geno_file, Pheno_IID, by = "IID")

#Seperate the Genotype data based on the breed 
Pen_Geno_data <- filter(data_file, cage.pen == "PEN")
Cage_Geno_data <- filter(data_file, cage.pen == "CAGE")

#Remove phenotypes after merger to retain only the Genotype data
Pen_Geno_data <- Pen_Geno_data %>% dplyr::select(-cage.pen, -comb_mass)
Cage_Geno_data <- Cage_Geno_data %>% dplyr::select(-cage.pen,-comb_mass)

#Seperate the Phenotype data based on the breed
Pen_Pheno_data <- subset(Pheno_file, cage.pen == "PEN")
Cage_Pheno_data <- subset(Pheno_file, cage.pen =="CAGE")

#Seperate the covariate data based on breeds
pheno_covariate_Pen <- subset(pheno_covariate, cage.pen == "PEN")
pheno_covariate_Cage <- subset(pheno_covariate, cage.pen == "CAGE")

#Individual Breed GWAS analysis
# Prune SNPs
geno_pruned_Pen <- prune_snp_matrix(Pen_Geno_data)
geno_pruned_Cage <- prune_snp_matrix(Cage_Geno_data)

# Construct GRM
G_P <- construct_grm(geno_pruned_Pen)
G_C <- construct_grm(geno_pruned_Cage)

# Decomposition of the GRM to create an incidence matrix for the additive genetic effect
Z_grm_P <- decompose_grm(G_P, FALSE)
Z_grm_C <- decompose_grm(G_C, FALSE)

# Add an intercept column (a column of ones) to the covariate matrix
X_Pen <-model.matrix(~ 1 + weight + feed + breed, pheno_covariate_Pen)
Z_group_P <- model.matrix(~ factor(group), pheno_covariate_Pen)

X_Cage <-model.matrix(~ 1 + weight + feed + breed, pheno_covariate_Cage)
Z_group_C <- model.matrix(~ factor(group), pheno_covariate_Cage)

# Z_grm_P and Z_group_P are matrices or data frames
num_cols_Z_grm_P <- ncol(Z_grm_P)
num_cols_Z_group_P <- ncol(Z_group_P)

num_cols_Z_grm_C <- ncol(Z_grm_C)
num_cols_Z_group_C <- ncol(Z_group_C)

# Number of columns for each random effect; here only one
RandC_P <- c(num_cols_Z_grm_P, num_cols_Z_group_P)
RandC_C <- c(num_cols_Z_grm_C, num_cols_Z_group_C)

# Run GWAS for individual breeds
gwas_P <- run_gwasA(pheno = Pen_Pheno_data$comb_mass, X = X_Pen, Z = cbind(Z_grm_P, Z_group_P), RandC = RandC_P, geno_pruned_Pen)
gwas_C <- run_gwasA(pheno = Cage_Pheno_data$comb_mass, X = X_Cage, Z = cbind(Z_grm_C,Z_group_C), RandC = RandC_C, geno_pruned_Cage)

# Placeholder for marker positions
marker_pos <- read.table("A:/Books/SEM 4/Thesis/GWAS/Script/Work_data/Input_files/genotypes.map")
colnames(marker_pos) <- c("chromosome", "marker", "distance", "position")

#Removing extra characters to make markerID uniform in all files
gwas_P$marker_id <- sub("_[ATGC]$", replacement = "", gwas_P$marker_id)
gwas_C$marker_id <- sub("_[ATGC]$", replacement = "", gwas_C$marker_id)

#Join the output to the markerID file
gwas_pos_P <- inner_join(gwas_P, marker_pos, by = c("marker_id" = "marker"))
gwas_pos_C <- inner_join(gwas_C, marker_pos, by = c("marker_id" = "marker"))

# Look at significant positions in the two scans
sig_gwas_pos_P <- filter(gwas_pos_P, p < 1e-4)
sig_gwas_pos_C <- filter(gwas_pos_C, p < 1e-4)

#To highlight sig_SNP in Man. plot
sig_pos_P <- sig_gwas_pos_P$marker_id
sig_pos_C <- sig_gwas_pos_C$marker_id

# Prepare the data for qqman
gwas_results_P <- data.frame(SNP = gwas_pos_P$marker_id, CHR = gwas_pos_P$chromosome, BP = gwas_pos_P$position, P = gwas_pos_P$p)
gwas_results_C <- data.frame(SNP = gwas_pos_C$marker_id, CHR = gwas_pos_C$chromosome, BP = gwas_pos_C$position, P = gwas_pos_C$p)

gwas_results_P$CHR <- as.numeric(gwas_results_P$CHR)
gwas_results_C$CHR <- as.numeric(gwas_results_C$CHR)

gwas_results_P <- gwas_results_P %>%filter(!is.na(CHR) & !is.na(BP) & !is.na(P))
gwas_results_C <- gwas_results_C %>%filter(!is.na(CHR) & !is.na(BP) & !is.na(P))

# Generate Manhattan plot
manhattan(gwas_results_P,main = "Manhattan Plot: Pen", highlight = sig_pos_P, col = c("blue4", "orange3"), annotatePval = 0.0001, genomewideline = -log10(5e-8), suggestiveline = -log10(1e-4),chrlabs = c(1:22, "X", "Y"))
manhattan(gwas_results_C,main = "Manhattan Plot: Cage", highlight = sig_pos_C, annotatePval = 0.0001, genomewideline = -log10(5e-8), suggestiveline = -log10(1e-4),chrlabs = c(1:22, "X", "Y"))

# Create Manhattan plot for gwasA_pos (Pen)
manhattan_gwas_P <- ggplot(gwas_pos_P, aes(x = as.numeric(factor(chromosome)), y = -log10(p))) +
  geom_point(aes(color = as.factor(chromosome)), size = 1.5) +  # Color by chromosome
  labs(title = "Manhattan Plot: Pen", x = "Chromosome", y = "-log10(p)") +
  geom_hline(yintercept = -log10(1e-4), linetype = "dotted", color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot for Pen
print(manhattan_gwas_P)

# Create Manhattan plot for gwasA_pos (Cage)
manhattan_gwas_C <- ggplot(gwas_pos_C, aes(x = as.numeric(factor(chromosome)), y = -log10(p))) +
  geom_point(aes(color = as.factor(chromosome)), size = 1.5) +  # Color by chromosome
  labs(title = "Manhattan Plot: Cage", x = "Chromosome", y = "-log10(p)") +
  geom_hline(yintercept = -log10(1e-4), linetype = "dotted", color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot for Cage
print(manhattan_gwas_C)

# Display both plots side by side
grid.arrange(manhattan_gwas_P, manhattan_gwas_C, ncol = 2)

#Create qqplot
qq(gwas_results_P$P, main = "Q-Q plot of GWAS: Pen p-values", xlim = c(0, 6), ylim = c(0,8), pch = 19, col = "blue4", cex = 1, las = 1)
qq(gwas_results_C$P, main = "Q-Q plot of GWAS: Cage p-values", xlim = c(0, 6), ylim = c(0,8), pch = 19, col = "blue4", cex = 1, las = 1)


# Step 1: Calculate chi-squared statistics
gwas_results_P$chi_sq <- qchisq(gwas_results_P$P, df=1, lower.tail=FALSE)
gwas_results_C$chi_sq <- qchisq(gwas_results_C$P, df=1, lower.tail=FALSE)

# Step 2: Calculate the observed and expected median of chi-squared statistics
observed_median_P <- median(gwas_results_P$chi_sq)
expected_median_P <- qchisq(0.5, df=1) # This is 0.456 for chi-squared with 1 degree of freedom

observed_median_C <- median(gwas_results_C$chi_sq)
expected_median_C <- qchisq(0.5, df=1) # This is 0.456 for chi-squared with 1 degree of freedom

# Step 3: Calculate lambda genomic inflation factor
lambda_P <- observed_median_P/ expected_median_P
lambda_C <- observed_median_C/ expected_median_C

# Print the lambda value
print(paste("Genomic inflation factor (lambda):", lambda_P))
print(paste("Genomic inflation factor (lambda):", lambda_C))



