# Pain_RheumatoidArthritic_Depression_Fusion
Codes for the paper “Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.” 
E-mail: qianwang_bnu@mail.bnu.edu.cn
This repository provides R code for statistical modeling, survival analysis, mediation analysis, and Mendelian randomization (MR) supporting the findings of the study investigating the multidimensional links among rheumatoid arthritis (RA), depression and multi site chronic pain (MCP).
Code

0_anova_linear_analysis.R
One-way ANOVA and polynomial trend modeling assessing group differences in MCP-related multi-omic signatures across healthy, depression, RA, and comorbid groups.

1_cox.R
Cox proportional hazards modeling evaluating bidirectional associations between baseline RA and incident depression (and vice versa), including Schoenfeld residual tests, age stratification, VIF evaluation and logistic regression sensitivity analyses.

2_cox.R
Cox survival models testing whether single MCP-derived omic signatures prospectively predict risk of depression or RA, with full covariate adjustment and diagnostic checks.

3_mediation_loadings.R
Mediation analysis evaluating whether MCP-related blood signatures mediate the bidirectional RA ↔ depression relationship (bootstrap = 5,000 iterations, ACME/ADE/total effect estimation + proportion mediated).

4_protein_MR_analysis.R
Two-sample Mendelian randomization for RA → protein traits (and bidirectionally for other omics). Implements IVW (primary), MR-Egger, weighted/median estimators, heterogeneity tests, MR-PRESSO outlier correction, and automated looping over 119 proteins.

5_MR_plot_forest.R
Generation of publication-ready forest plots summarizing MR causal estimates (OR + 95% CI), formatted for inclusion in main or supplementary figures.


🔒 Data Availability Statement

All analyses were conducted using UK Biobank genetic, phenotypic and biomarker data, which require institutional approval and formal data access agreements.

Due to data privacy restrictions, individual-level data cannot be shared.
To support reproducibility, this repository includes:

✔ Complete analysis scripts
✔ Code for preprocessing, modeling, MR, mediation, and plotting
✔ Variable requirements and data structure documentation inside each script

Researchers with approved access to UK Biobank can fully reproduce the pipeline and results.
