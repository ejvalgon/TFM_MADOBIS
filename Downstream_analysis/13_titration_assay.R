# ============================================================
# SPIKE-IN LINEAR REGRESSION PLOT
# ============================================================
# This script reads the spike-in read counts obtained from idxstats,
# extracts the known spike-in percentage from the sample name,
# fits a linear regression model, and plots the relationship between
# spike-in percentage and spike-in reads using ggplot2.
#
# The goal is to check whether the number of spike-in reads increases
# linearly with the amount of spike-in chromatin added to each sample.
# ============================================================


# ===============================
# LOAD REQUIRED LIBRARIES
# ===============================

library(tidyverse)


# ===============================
# READ SPIKE-IN COUNTS TABLE
# ===============================
# The input table should contain one row per sample.
# Expected columns:
#   sample       -> sample name, containing the SPIKE percentage
#   spike_reads  -> number of reads mapped to the spike-in genome

scale_factor_df <- read_tsv("../Script_Limpios/Titration_assay/scale_factors.tsv")


# ===============================
# EXTRACT SPIKE-IN PERCENTAGE
# ===============================
# The spike-in percentage is encoded in the sample name after the word "SPIKE".
#
# Examples:
#   SPIKE1_S6     -> 1
#   SPIKE5_1_S4   -> 5
#   SPIKE10_1_S2  -> 10
#   SPIKE20_S1    -> 20
#
# The regular expression extracts only the numeric value immediately
# following "SPIKE".

df_plot <- scale_factor_df %>%
  mutate(
    spike_percentage = str_extract(sample, "(?<=SPIKE)\\d+") %>% as.numeric()
  ) %>%
  arrange(spike_percentage)


# ===============================
# FIT LINEAR REGRESSION MODEL
# ===============================
# The model tests whether spike-in reads scale linearly with the known
# spike-in percentage.
#
# Model:
#   spike_reads = intercept + slope * spike_percentage

model <- lm(spike_reads ~ spike_percentage, data = df_plot)


# ===============================
# EXTRACT MODEL STATISTICS
# ===============================
# Extract the slope, intercept, and R-squared value so they can be shown
# directly on the plot.

intercept <- coef(model)[1]
slope <- coef(model)[2]
r2 <- summary(model)$r.squared

label_text <- paste0(
  "y = ", round(slope, 2), "x + ", round(intercept, 2),
  "\nR² = ", round(r2, 4)
)


# ===============================
# PLOT LINEAR REGRESSION
# ===============================
# geom_point() shows the observed spike-in read counts.
# geom_smooth(method = "lm") adds the fitted linear regression line.
# The confidence interval is disabled with se = FALSE to keep the plot clean.

p_spike_regression <- ggplot(df_plot, aes(x = spike_percentage, y = spike_reads)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = label_text,
    hjust = 1.1,
    vjust = -0.8,
    size = 5
  ) +
  labs(
    x = "Spike-in chromatin percentage",
    y = "Spike-in reads"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 20)
  ) +
  theme_classic(base_size = 20)

p_spike_regression


# ===============================
# SAVE PLOT
# ===============================

dir.create("../Script_Limpios/Titration_assay/plots", showWarnings = FALSE)

ggsave(
  filename = "../Script_Limpios/Titration_assay/plots/spikein_linear_regression_ggplot.png",
  plot = p_spike_regression,
  width = 6,
  height = 5,
  dpi = 2400
)

ggsave(
  filename = "../Script_Limpios/Titration_assay/plots/spikein_linear_regression_ggplot.jpg",
  plot = p_spike_regression,
  width = 6,
  height = 5,
  dpi = 2400
)