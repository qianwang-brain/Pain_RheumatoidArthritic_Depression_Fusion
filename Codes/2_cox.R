###############################################################################
## Script: 2_Cox.R
## Purpose:
##   Visualization of Cox proportional hazards results.
##   This section generates a forest-style plot showing hazard ratios (HR) and
##   their 95% confidence intervals for the main exposure variables used in the
##   Cox models (e.g., RA_label, dep_label, or MCP-related blood signatures).

### Reference:
#   Multi-Site Chronic Pain Reveals Shared Neuro-Immune-Metabolic Alterations Underlying Rheumatoid Arthritis and Depression.
###############################################################################

##----------------Plot
library(ggplot2)
library(dplyr)
library(forcats)

data <- data %>%
  mutate(variable = fct_rev(fct_inorder(variable)))

x_limits <- c(floor(min(data$lower)*10)/10, ceiling(max(data$upper)*10)/10)

ggplot(data, aes(x = HR, y = variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper, color = variable),
    height = 0, 
    linewidth = 1.2, 
    lineend = "butt"
  ) +
  geom_point(
    aes(fill = variable), 
    size = 3.5, 
    shape = 23, 
    color = "white"
  ) +
  scale_x_continuous(
    breaks = seq(x_limits[1], x_limits[2], 0.1),
    limits = x_limits,
    expand = expansion(mult = 0.02),
    name = "Hazard Ratio (95% CI)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12), 
    legend.position = "none"
  ) +
  scale_color_hue(l = 40) +
  scale_fill_hue(l = 40)

ggsave("E:/Project_R/R_code/UKB_project_RA/Cox_results/before_mediation1/cox1.pdf", 
       width = 9, 
       height = 0.7 * nrow(data),
       device = "pdf", 
       units = "in",
       dpi = 300)





















