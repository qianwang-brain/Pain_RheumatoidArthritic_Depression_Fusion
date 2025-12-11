###########################################################################
# Script: 0_anova_linear_analysis.R
#
# Purpose:
#   Linear gradients of MCP-related multi-omic signatures across disease burden.
#   To formally assess how the identified multi-omic signatures differ across
#   groups, the script performs a one-way ANOVA followed by polynomial trend
#   analysis.
#
# Description:
#   – Loads group-stratified multi-omic results from a CSV file.
#   – Reorders the group factor into an ordered structure:
#        healthy < Depression < RA < Comorbid.
#   – Performs a one-way ANOVA to evaluate overall group differences.
#   – Applies polynomial contrasts (linear, quadratic, cubic) to quantify
#        systematic trends across disease burden levels.
#   – Extracts linear/quadratic/cubic trend statistics from the model
#        using summary.lm for clearer interpretation of gradient effects.
#
# Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###########################################################################


library(readr)
library(dplyr)

df <- read_csv("C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/3_MCCAR2_results_2/results_plot/Linear_gradient/Blood1_ttest_all.csv")

df$group <- factor(
  trimws(df$group),
  levels = c("healthy", "Depression", "RA" ,"Comorbid"),
  ordered = TRUE
)

anova_res <- aov(a1_h ~ group, data = df)
summary(anova_res)


## 3) polynomial contrasts-------------------------------
contrasts(df$group) <- contr.poly(nlevels(df$group))
anova_trend <- aov(a1_h ~ group, data = df)

## (1) ANOVA
summary(anova_trend)
## (2)
summary.lm(anova_trend)

