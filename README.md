# GWAS Chicken Comb Analysis Project

# Description

The study was aimed to identify genetic variants associated with comb traits in commercial laying hens using Genome-Wide Association Studies (GWAS). Comb traits are important phenotypic indicators linked to thermoregulation, reproductive health, and overall fitness in poultry.

This project involves processing genotype data, performing quality control, and applying statistical association models to detect significant SNPs. Post-GWAS analyses are to be conducted to interpret the functional relevance of identified variants and understand their contribution to phenotypic variation.

# Analyses Performed
## 1. Genotype Data Preprocessing
  > Genotype data was formatted and prepared for analysis by applying filtering criteria to remove SNPs and samples with high missingness, low minor allele frequency (MAF), to retain high quality data.

## 2. Quality Control (QC)
 > Quality control was performed to ensure data reliability by filtering low-quality samples and variants, identifying outliers, and assessing population structure to minimize confounding effects in the analysis

## 3. Genome-Wide Association Study (GWAS)
 > Association testing was performed between SNPs and comb traits using statistical models. Linear mixed models was implemented to account for population stratification. Significantly associated genome loci were identified.

## 4. Data Integration and Interpretation
 > GWAS results was integrated with genome annotations using biological evidences. Candidate genes and pathways were identified that are associated with the genes. The results were interpreted in the context of poultry genetics and breeding.

## 5. Visualization and Results
 > Generation of Manhattan plots and QQ plots for GWAS results to visualise significant loci and genomic regions. Summary statistics were also plotted.

## 6. Reproducibility and Workflow Management
 > Workflow were documented in a structured scripts for reproducible analysis and for transparency, ensuring consistency in data processing and analysis steps.

 > Large genotype files (e.g., .ped, .raw) are excluded from this repository due to file size limitations. Only scripts, workflows, and processed summary results are included. Raw data can be made available upon reasonable request.

## Directory Tree
 <img width="383" height="286" alt="image" src="https://github.com/user-attachments/assets/af3093d3-32e1-44a3-b616-66fe0190d3e6" />

## Author
 ### Prasanna Venkataramanan

## License

This project is intended for academic and research purposes only.

## Acknowledgements
  1. I would like to thank Uppsala University for providing the academic environment and resources to support this research.
  2. I would like to thank my project supervisor [Dr. Martin Johnsson](https://github.com/mrtnj) for guidance and support throughout this project.

