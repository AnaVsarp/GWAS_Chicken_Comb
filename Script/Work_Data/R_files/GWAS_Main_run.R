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

# Load the raw genotype file form PLINK and othe phenotype files
Geno_1 <- read.table("A:/Books/SEM 4/GWAS_1/genotypes_qc.raw", header = TRUE)
pheno <- read.csv("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/GWAS_files/Pheno_file.csv", header = TRUE)
comb_g <- read.csv("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/Input_files/comb_g.csv",header = TRUE)

#Source function file
source("A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/R_files/GWAS_Functions.R")

#Filter the samples whose genotype data is available
pheno_filter <- dplyr::filter(pheno, IID %in% Geno_1$IID)

#Check if the data is oriented in a same order in all the data files
identical(Geno_1$IID, pheno$IID)
#order(Geno_1$IID)
#order(pheno$animal_id)
Geno_2<- Geno_1[order(Geno_1$IID),]
pheno_2<- pheno_filter[order(pheno_filter$IID),]
identical(Geno_2$IID, pheno_2$IID)

#Attach comb data to phenotype and genotype files to remove all the non combmass data samples
comb <- dplyr::select(comb_g, IID)
Pheno_file <- dplyr::inner_join(pheno_2, comb, by = c("IID" = "IID"))
Geno_Main <- inner_join(Geno_2, comb, by = c("IID"= "IID"))
#Geno_join <- subset(Geno_join, select = -c(comb_mass))
identical(Geno_Main$IID, Pheno_file$IID)

#Form covariates file
pheno_covariate <- dplyr::select(Pheno_file, IID, feed, cage.pen, breed, weight, comb_mass, group)

# Prune SNPs
geno_pruned <- prune_snp_matrix(Geno_Main)

# Construct GRM
G <- construct_grm(geno_pruned)

# Decomposition of the GRM to create an incidence matrix for the additive genetic effect
Z_grm <- decompose_grm(G, FALSE)

# Add an intercept column (a column of ones) to the covariate matrix
X <-model.matrix(~ 1 + cage.pen + breed + weight + feed , pheno_covariate)

# Number of columns for each random effect; here only one
RandC <- ncol(Z_grm)

# Run GWAS without dominance term
gwas <- run_gwasA(pheno = Pheno_file$comb_mass, X = X, Z = Z_grm, RandC = RandC, geno_pruned)

write.csv(gwas, "A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/GWAS_Files/gwas.csv", row.names = FALSE )

# Placeholder for marker positions
marker_pos <- read.table("A:/Books/SEM 4/Thesis/GWAS/Script/Work_data/Input_files/genotypes.map")
colnames(marker_pos) <- c("chromosome", "marker", "distance", "position")

gwas$marker_id <- sub("_[ATGC]$", replacement = "", gwas$marker_id)

# Add positions to results
gwas_pos <- inner_join(gwas, marker_pos, by = c("marker_id" = "marker"))

# Look at significant positions in the two scans
sig_gwas_pos <- filter(gwas_pos, p < 1e-5)

#To highlight sig_SNP in Man. plot
sig_pos <-sig_gwas_pos$marker_id

# Prepare the data for qqman
gwas_results <- data.frame(SNP = gwas_pos$marker_id, CHR = gwas_pos$chromosome, BP = gwas_pos$position, P = gwas_pos$p)
gwas_results$CHR <- as.numeric(gwas_results$CHR)
gwas_results <- gwas_results %>%filter(!is.na(CHR) & !is.na(BP) & !is.na(P))

write.csv( gwas_results, "A:/Books/SEM 4/Thesis/GWAS/Script/Work_Data/GWAS_Files/gwas_results.csv", row.names = FALSE )

# Calculate FDR adjusted p-values
gwas_results <- data.frame(SNP = gwas_pos$marker_id, CHR = gwas_pos$chromosome, BP = gwas_pos$position, P = gwas_pos$p)
gwas_results$CHR <- as.numeric(gwas_results$CHR)
gwas_results <- gwas_results %>% filter(!is.na(CHR) & !is.na(BP) & !is.na(P))

# Generate Manhattan plot
manhattan(gwas_results,main = "Manhattan Plot", highlight = sig_pos, col = c("blue4", "orange3"), annotatePval = 0.0001, genomewideline = -log10(5e-8), suggestiveline = -log10(1e-5),chrlabs = c(1:22, "X", "Y"))


# Create Manhattan plot for gwas_pos by chromosome
manhattan_gwas <- ggplot(gwas_pos, aes(x = as.numeric(factor(chromosome)), y = -log10(p))) +
  geom_point(aes(color = as.factor(chromosome)), size = 1.5) +  # Color by chromosome
  labs(title = "Manhattan Plot", x = "Chromosome", y = "-log10(p)") +
  geom_hline(yintercept = -log10(1e-5), linetype = "dotted", color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(manhattan_gwas)

# Combine all plots into one layout
#combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)

# Print the combined plot
#print(combined_plot)

#Create qqplot
qq(gwas_results$P, main = "Q-Q plot of GWAS p-values", xlim = c(0, 6), ylim = c(0,8), pch = 19, col = "blue4", cex = 1, las = 1)

# Step 1: Calculate chi-squared statistics
gwas_results$chi_sq <- qchisq(gwas_results$P, df=1, lower.tail=FALSE)

# Step 2: Calculate the observed and expected median of chi-squared statistics
observed_median <- median(gwas_results$chi_sq)
expected_median <- qchisq(0.5, df=1) # This is 0.456 for chi-squared with 1 degree of freedom

# Step 3: Calculate lambda genomic inflation factor
lambda <- observed_median / expected_median

# Print the lambda value
print(paste("Genomic inflation factor (lambda):", lambda))