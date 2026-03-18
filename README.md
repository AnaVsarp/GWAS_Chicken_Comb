# GWAS Chicken Comb Analysis Project

# Description

The study was aimed to identify genetic variants associated with comb traits in commercial laying hens using Genome-Wide Association Studies (GWAS). Comb traits are important phenotypic indicators linked to thermoregulation, reproductive health, and overall fitness in poultry.

This project involves processing genotype data, performing quality control, and applying statistical association models to detect significant SNPs. Post-GWAS analyses are to be conducted to interpret the functional relevance of identified variants and understand their contribution to phenotypic variation.

# Analysis to Perform:
## 1. Genotype Data Preprocessing:
  > Genotype data was formatted and prepared for analysis by applying filtering criteria to remove SNPs and samples with high missingness, low minor allele frequency (MAF), to retain high quality data.

## 2. Quality Control (QC):
 > Quality control was performed to ensure data reliability by filtering low-quality samples and variants, identifying outliers, and assessing population structure to minimize confounding effects in the analysis

3. Genome-Wide Association Study (GWAS):
 > Association testing was performed between SNPs and comb traits using statistical models. Linear mixed models was implemented to account for population stratification. Significantly associated genome loci were identified.

4. Post-GWAS Analysis:

Fine-mapping of significant regions to pinpoint candidate variants.
Annotation of SNPs to identify nearby genes and regulatory elements.
Functional interpretation of genetic variants influencing comb traits.

5. Data Integration and Interpretation:

Integration of GWAS results with biological knowledge and genomic annotations.
Identification of candidate genes and pathways associated with the phenotype.
Interpretation of results in the context of poultry genetics and breeding.

6. Visualization and Results:

Generation of Manhattan plots and QQ plots for GWAS results.
Visualization of significant loci and genomic regions.
Summary statistics and graphical representation of findings.

7. Reproducibility and Workflow Management:

Structured pipeline for reproducible analysis.
Use of scripts and documented workflows for transparency.
Ensuring consistency in data processing and analysis steps.

⚠️ Data Availability:

Large genotype files (e.g., .ped, .raw) are excluded from this repository due to file size limitations.
Only scripts, workflows, and processed summary results are included.
Raw data can be made available upon reasonable request.

Author:

Prasanna Venkat
Bioinformatics | Genomics | GWAS Analysis

📄 License:

This project is intended for academic and research purposes only.

