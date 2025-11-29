###############################################################################
## Script: 3_mediation_loadings.R
##
## Purpose:
##   Mediation analyses.
##   This script evaluates whether MCP-related blood signatures mediate the
##   bidirectional association between RA and depression, focusing here on
##   the pathway:
##       prevalent baseline RA status → blood signature weights1 → incident depression.
##
## Methods overview:
##   – Load the mediation dataset for the RA → depression direction.
##   – Specify:
##        * Outcome model (c): parametric survival model (survreg) with
##          follow-up time and incident depression status as the outcome,
##          RA_label as the exposure, and weights1 as a putative mediator,
##          adjusted for age, sex, BMI, alcohol consumption, smoking,
##          ethnicity, household income, university education, deprivation
##          index, and metabolic syndrome.
##        * Mediator model (b): linear regression of weights1 on RA_label and the
##          same covariate set.
##   – Apply the `mediate()` function (mediation R package) with 5,000
##     bootstrap iterations (sims = 5000, boot = TRUE) to estimate:
##        * total effect (tau),
##        * average direct effect (ADE, d.avg),
##        * average indirect / mediated effect (ACME, z.avg),
##        * proportion mediated (n.avg),
##       along with 95% confidence intervals and p-values.
##   – Store key mediation estimates in a compact matrix and export them
##     as a CSV file for downstream aggregation and multiple testing
##     correction using FDR, consistent with the main analysis plan.
##   – Derive percentage of effect mediated (pm) and its confidence
##     interval (lower / upper) for reporting in tables and figures.
##
# Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###############################################################################

library(mediation)
library(survival)
library(survminer)

datas <- read.csv('C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation1/data_mediation1_blood.csv')

c <- survreg(Surv(time, status) ~ RA_label+L1+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
summary(c)


b<-lm(L1~RA_label+age+gender+bmi+alcohol+smooking+ethnicity+income+university_education+deprivation+metabolic_syndrome,data = datas)
summary(b)



# bootstrapping test, it's better to use the bootstrapping test, if the treat group was not specificied, the default is '1'
contcont <- mediate(b, c, sims=5000,boot=T, treat="RA_label", mediator="L1",outcome='time')

# 输出中介模型结果
print(contcont)

a<-matrix(nrow=1,ncol=20)
a[1,1]=contcont$tau.coef
a[1,2]=contcont$tau.ci[1]
a[1,3]=contcont$tau.ci[2]
a[1,4]=contcont$tau.p
a[1,5]=0
a[1,6]=contcont$d.avg
a[1,7]=contcont$d.avg.ci[1]
a[1,8]=contcont$d.avg.ci[2]
a[1,9]=contcont$d.avg.p
a[1,10]=0
a[1,11]=contcont$z.avg
a[1,12]=contcont$z.avg.ci[1]
a[1,13]=contcont$z.avg.ci[2]
a[1,14]=contcont$z.avg.p
a[1,15]=0
a[1,16]=contcont$n.avg
a[1,17]=contcont$n.avg.ci[1]
a[1,18]=contcont$n.avg.ci[2]
a[1,19]=contcont$n.avg.p
a[1,20]=nrow(datas)

print(a[1,16])
write.csv(a,'C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation1/results_mediation1_L1_allcov.csv')

#
result<-read.csv('C:/Program Files/MATLAB/R2020b/bin/ukb_RA/pain_guided_blood/4_Mediation/mediation1/results_mediation1_L1_allcov.csv',header=T)
result$pm=100*result$V16
result$lower=100*result$V17
result$upper=100*result$V18

