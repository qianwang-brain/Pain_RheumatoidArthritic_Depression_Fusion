###############################################################################
## Script: 4_protein_MR_analysis.R

## Purpose:
##   Causal inference using Mendelian Randomization (MR).
##   This script conducts bidirectional two-sample MR analyses to examine
##   whether genetically predicted MCP-related multi-omic signatures exert
##   causal effects on RA or depression. The section shown here corresponds to
##   the direction: RA → protein traits, iterating over 119 plasma proteins.
##
## Methods overview:
##   – Summary statistics for RA (exposure) and plasma protein GWAS (outcomes)
##     are imported from independent datasets as described in Supplementary
##     Table 11.
##   – Exposure and outcome datasets are harmonised to ensure allele alignment,
##     allowing automatic strand flips when required.
##   – MR is performed using:
##        * Inverse-variance weighted (IVW, primary)
##        * Weighted median
##        * Simple median
##        * MR-Egger regression
##        * Simple mode / weighted mode
##     Odds ratios (ORs) and confidence intervals are extracted for all methods.
##
##   – For each of the 119 proteins, the script:
##        * Reads GWAS summary statistics,
##        * Harmonises exposure–outcome SNPs,
##        * Performs MR and heterogeneity tests,
##        * Optionally runs MR-PRESSO,
##        * Saves all MR results into protein-specific CSV files,
##        * Updates a summary table recording final IVW p-values and any
##          associated error or pleiotropy flags.
## Notes:
##   – All analyses follow the same framework as described in the Methods:
##       (1) SNP extraction at genome-wide significance,
##       (2) LD clumping,
##       (3) harmonisation,
##       (4) IVW plus complementary MR estimators,
##       (5) heterogeneity & pleiotropy assessment,
##       (6) MR-PRESSO correction when applicable.
##   – This script corresponds to the RA → protein direction; parallel scripts
##     are used for depression → proteins and multi-omic layers (metabolites,
##     blood cell traits, biochemistry markers).
# Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###############################################################################




rm(list = ls())
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(MRPRESSO)
protein_names <- fread("E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/protein_name.csv",
                       header = TRUE)
names(protein_names) <- "protein"

results_summary <- data.table(
  protein = character(),
  pval_method3 = numeric(),
  error_type = character()
)

exposure <- fread("E:/Gene_Dataset/RA_MR/finngen_R11_M13_RHEUMA/finngen_R11_M13_RHEUMA",
                  header = TRUE)
exposure <- subset(exposure, pval < 5e-8)

write.csv(
  exposure,
  'C:/Users/qianw/AppData/Local/R/win-library/4.4/TwoSampleMR/exposure.csv'
)
exposure <- system.file('exposure.csv', package = "TwoSampleMR")

exposure_exp_dat <- read_exposure_data(
  filename = exposure,
  sep = ",",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval"
)

exposure_dat_clumped <- clump_data(
  exposure_exp_dat,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)
rm(exposure, exposure_exp_dat)


for (current_protein in protein_names$protein) {
  tryCatch({
    outcome_file <- sprintf("merged_%s_with_rsid.regenie", current_protein)
    outcome_path <- file.path("E:/Gene_Dataset/UKB_protein_GWAS/protein_merge",
                              outcome_file)
    outcome <- fread(outcome_path, header = TRUE)
    outcome$P <- 10 ^ (-outcome$LOG10P)
    c <- merge(exposure_dat_clumped,
               outcome,
               by.x = "SNP",
               by.y = "ID")
  
    write.csv(c,
              'C:/Users/qianw/AppData/Local/R/win-library/4.4/TwoSampleMR/outcome.csv')
    rm(c)
    outcome <- system.file('outcome.csv', package = "TwoSampleMR")
    outcome_dat <- read_outcome_data(
      snps = exposure_dat_clumped$SNP,
      filename = outcome,
      sep = ",",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "P"
    )
    
    dat <- harmonise_data(exposure_dat = exposure_dat_clumped,
                          outcome_dat = outcome_dat)
    result =mr(dat, method_list = c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe","mr_egger_regression_bootstrap",
                                    "mr_simple_median","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
    or <- generate_odds_ratios(result)
    output_path <- sprintf(
      "E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/results_protein2/%s_RA.csv",
      current_protein
    )
    write.csv(or, output_path, row.names = FALSE)
    
    
    hetegeneity<-mr_heterogeneity(dat)
    q_pval <- hetegeneity %>%
      filter(method == "Inverse variance weighted") %>%
      pull(Q_pval) %>%
      signif(3)
    
    if (length(q_pval) == 1 && !is.na(q_pval)) {
      ivw_method <- ifelse(q_pval < 0.05,
                           "Inverse variance weighted (multiplicative random effects)",
                           "Inverse variance weighted (fixed effects)")
      
      ivw_pval <- or %>%
        filter(method == ivw_method) %>%
        pull(pval)
      
      if (length(ivw_pval) == 0) {
        ivw_pval <- NA
        error_reason <- paste0("OR 表中未找到 ", ivw_method)
      } else {
        error_reason <- ""
      }
      
    } else {
      ivw_pval <- NA
      error_reason <- "异质性检验 Q 值缺失或无效"
    }
    
    pleitropy<-mr_pleiotropy_test(dat)
    
    if (!is.na(ivw_pval) && ivw_pval < 0.05 && pleitropy$pval < 0.05) {
      mr_presso_results <- mr_presso(BetaOutcome = "beta.outcome",
                                     BetaExposure = "beta.exposure",
                                     SdOutcome = "se.outcome",
                                     SdExposure = "se.exposure",
                                     OUTLIERtest = TRUE,
                                     DISTORTIONtest = TRUE,
                                     data = dat,
                                     NbDistribution = 5000,
                                     SignifThreshold = 0.05)

      outlier_indices <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      outlier_snps <- dat$SNP[outlier_indices]
    
      dat_clean <- subset(dat, !(SNP %in% outlier_snps))
    
      result_clean <- mr(dat_clean, method_list = c("mr_egger_regression","mr_ivw_mre","mr_ivw_fe","mr_egger_regression_bootstrap",
                                                    "mr_simple_median","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
      or_clean <- generate_odds_ratios(result_clean)
      output_path <- sprintf(
        "E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/results_protein2/upgrade_according to Preitropy/%s_RA.csv",
        current_protein
      )
      write.csv(or_clean, output_path, row.names = FALSE)
      
      ivw_pval <- or_clean %>%
        filter(method == ivw_method) %>%
        pull(pval)
      
    } else {
      message("未触发 MR-PRESSO（因 IVW 或 pleiotropy 不显著）")
    }
    
    # 更新汇总表
    results_summary <- rbindlist(list(
      results_summary,
      list(
        protein = current_protein,
        pval_method3 = ifelse(length(ivw_pval) == 1, ivw_pval, NA),
        error_type = error_reason
      )
    ))
  }, error = function(e) {
    results_summary <<- rbindlist(list(
      results_summary,
      list(
        protein = current_protein,
        pval_method3 = NA,
        error_type = e$message
      )
    ))
    message(sprintf("[错误] %s: %s", current_protein, e$message))
    
  })
}

output_summary_path <- "E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/results_protein2/summary_pvals_RA_protein.csv"

write.csv(
  results_summary,
  file = output_summary_path,
  row.names = FALSE,
  fileEncoding = "GBK" 
)


