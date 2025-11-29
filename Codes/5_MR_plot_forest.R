###############################################################################
## Script: 5_MR_plot_forest.R
##
## Purpose:
##   To visualise Mendelian randomization (MR) results for RA → protein traits
##   as a publication-ready forest plot, summarising odds ratios (ORs) and
##   95% confidence intervals across all tested MCP-related protein outcomes.
##
# Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###############################################################################

library(data.table)
library(forestplot)
library(dplyr)

# 读取数据（示例路径需替换为实际路径）
data <- fread("E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/results_protein2/pval_matched_rows_Final.csv", 
              header = TRUE)

# 数据格式化
combined_data <- data %>%
  select(exposure,outcome, or, or_lci95, or_uci95, pval) %>%
  mutate(
    `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lci95, or_uci95),
    p_value = formatC(pval, format = "e", digits = 2)
  )

# 添加表头（网页1方法）
header <- tibble(
  exposure = "exposure",
  outcome="outcome",
  `OR (95% CI)` = "OR (95% CI)",
  p_value = "P Value"
)
combined_data <- bind_rows(header, combined_data)

cairo_pdf("E:/Project_R/R_code/UKB_project_RA/RA_MR_analysis/forestplot_ALL.pdf", 
          width = 8,       # 宽度12英寸（A4纸宽度）
          height = 11,       # 高度8英寸
          family = "Arial") # 设置英文字体

# 绘制森林图（网页2的绘图参数优化）
forestplot(
  labeltext = as.matrix(combined_data[, c("exposure","outcome","OR (95% CI)", "p_value")]),
  mean = c(NA, combined_data$or[-1]),
  lower = c(NA, combined_data$or_lci95[-1]),
  upper = c(NA, combined_data$or_uci95[-1]),
  graph.pos = 4,  # 调整图形位置到第4列[6](@ref)
  boxsize = 0.35,
  zero = 1,
  xlog = FALSE,
  clip = c(floor(min(combined_data$or_lci95, na.rm = TRUE)), 
           ceiling(max(combined_data$or_uci95, na.rm = TRUE))),
  col = fpColors(box = "#3c5a66", line = "black"),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.8),
    ticks = gpar(cex = 0.6),
    xlab = gpar(cex = 0.8)
  ),
  xlab = "OR (95% CI)",
  lwd.xaxis = 1.0,
  xticks = seq(1, 1.35, by = 0.05),  
)

dev.off() 




