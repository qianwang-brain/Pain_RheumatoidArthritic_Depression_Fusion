###############################################################################
## Script: 1_Cox.R
##
## Purpose:
##   Associations with the incidence of depression and RA.
##   This script implements Cox proportional hazards models to quantify the
##   bidirectional association between depression and RA, and to test whether
##   MCP-related blood signatures are prospectively associated with incident
##   depression or RA.
##
## Methods overview:
##   – Fit Cox proportional hazards models using follow-up time as the
##     underlying timescale:
##        * RA_label   → incident depression
##        * dep_label  → incident RA
##        * blood MCP-related signatures (e.g. L1) → incident outcomes.
##   – Stratify the baseline hazard by age strata (e.g. <65 vs ≥65) to relax
##     violations of the proportional hazards assumption.
##   – Adjust all models for age, sex, body mass index (BMI), alcohol
##     consumption, smoking status, ethnicity, household income, university
##     education, deprivation index, and metabolic syndrome.
##   – Extract hazard ratios (HRs), 95% confidence intervals, and p-values for
##     the main exposures (RA_label, dep_label, L1) into compact result tables
##     (result_RA, result_dep, result_L1) for downstream meta-analysis and FDR
##     multiple-testing correction (as described in Supplementary Methods 2).
##   – Evaluate the proportional hazards assumption using scaled and
##     rank-based Schoenfeld residual tests (cox.zph); no consistent or global
##     violations are expected.
##   – Assess multicollinearity among covariates using variance inflation
##     factors (VIF), verifying that inflation is low (VIF < ~1.8).
##   – Run logistic regression models with the same covariate set as a
##     robustness / sensitivity analysis against the time-to-event Cox models.
##
## Sections in this script:
##   RA → DEP      : incident depression as outcome, baseline RA_label as exposure.
##   DEP → RA      : incident RA as outcome, baseline dep_label as exposure.
##   Blood signatures:
##                   incident outcomes in relation to MCP-related blood
##                   signature scores analysed separately (e.g. L1).
# Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###############################################################################


##--------------------------------RA----DEP
library(mediation)
library(survival)
library(survminer)

datas <- read.csv('C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation1/1_RA_label_cox/data_cox.csv')

# 将年龄按临床意义或数据分布分层（示例：<65岁 vs ≥65岁）
datas$age_strata <- ifelse(datas$age < 65, "Younger", "Older")
datas$age_strata <- factor(datas$age_strata, levels = c("Younger", "Older"))


cox_model <- coxph(Surv(time, status) ~ RA_label+strata(age_strata)+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
summary(cox_model)

# 验证Cox模型的比例风险假设
test.ph <- cox.zph(cox_model)
print(test.ph)  # 若全局检验p>0.05，满足假设


test.ph <- cox.zph(cox_model,"rank")
print(test.ph)  # 若全局检验p>0.05，满足假设


# 提取特定变量的HR结果
hr_data <- summary(cox_model)$coefficients["RA_label", ]
confint_data <- summary(cox_model)$conf.int["RA_label", ]

# 创建结果数据框
result_RA <- data.frame(
  variable = "RA Status",
  HR = exp(coef(cox_model)["RA_label"]),
  lower = confint_data[3],
  upper = confint_data[4],
  pvalue = hr_data[5]
)

all_vars <- ls(all.names = TRUE)  # 包括隐藏变量（以`.`开头的变量）
rm_vars <- setdiff(all_vars, "result_RA")
rm(list = rm_vars)


## evaluate vif of predictiors
a=lm(status~RA_label+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
vif=car::vif(a)#aGSIF (the last column, named adjusted generalized standard error inflation factor) values above sqrt(2.5) may be of concern



##--------------------------------DEP----RA
library(mediation)
library(survival)
library(survminer)

datas <- read.csv('C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation2/dep_label_cox/data_cox.csv')

# 将年龄按临床意义或数据分布分层（示例：<65岁 vs ≥65岁）
datas$age_strata <- ifelse(datas$age < 65, "Younger", "Older")
datas$age_strata <- factor(datas$age_strata, levels = c("Younger", "Older"))


cox_model <- coxph(Surv(time, status) ~ dep_label+strata(age_strata)+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
summary(cox_model)

# 验证Cox模型的比例风险假设
test.ph <- cox.zph(cox_model)
print(test.ph)  # 若全局检验p>0.05，满足假设


test.ph <- cox.zph(cox_model,"rank")
print(test.ph)  # 若全局检验p>0.05，满足假设


# 提取特定变量的HR结果
hr_data <- summary(cox_model)$coefficients["dep_label", ]
confint_data <- summary(cox_model)$conf.int["dep_label", ]

# 创建结果数据框
result_dep <- data.frame(
  variable = "dep Status",
  HR = exp(coef(cox_model)["dep_label"]),
  lower = confint_data[3],
  upper = confint_data[4],
  pvalue = hr_data[5]
)

all_vars <- ls(all.names = TRUE)  # 包括隐藏变量（以`.`开头的变量）
rm_vars <- setdiff(all_vars, "result_dep")
rm(list = rm_vars)

## evaluate vif of predictiors
a=lm(status~dep_label+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
vif=car::vif(a)#aGSIF (the last column, named adjusted generalized standard error inflation factor) values above sqrt(2.5) may be of concern
max(vif)
##glm
logit_model <- glm(status ~ dep_label+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome, family = binomial(link = "logit"), data = datas)
summary(logit_model)




##Blood signatures separately
library(mediation)
library(survival)
library(survminer)

datas <- read.csv('C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation1/data_mediation1_blood.csv')
# 将年龄按临床意义或数据分布分层（示例：<65岁 vs ≥65岁）
datas$age_strata <- ifelse(datas$age < 65, "Younger", "Older")
datas$age_strata <- factor(datas$age_strata, levels = c("Younger", "Older"))


cox_model <- coxph(Surv(time, status) ~ L1+strata(age_strata)+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
summary(cox_model)

# 验证Cox模型的比例风险假设
test.ph <- cox.zph(cox_model)
print(test.ph)  # 若全局检验p>0.05，满足假设


test.ph <- cox.zph(cox_model,"rank")
print(test.ph)  # 若全局检验p>0.05，满足假设

# 提取特定变量的HR结果
hr_data <- summary(cox_model)$coefficients["L1", ]
confint_data <- summary(cox_model)$conf.int["L1", ]

# 创建结果数据框
result_L1 <- data.frame(
  variable = "L1",
  HR = exp(coef(cox_model)["L1"]),
  lower = confint_data[3],
  upper = confint_data[4],
  pvalue = hr_data[5]
)


all_vars <- ls(all.names = TRUE)  # 包括隐藏变量（以`.`开头的变量）
rm_vars <- setdiff(all_vars, c("result_RA","result_L1"))
rm(list = rm_vars)

## evaluate vif of predictiors
a=lm(status~L1+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
vif=car::vif(a)#aGSIF (the last column, named adjusted generalized standard error inflation factor) values above sqrt(2.5) may be of concern















