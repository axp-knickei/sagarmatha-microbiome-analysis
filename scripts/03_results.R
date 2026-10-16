# ============================================
# 03_results.R — Extract results, tables, figures
# ============================================
# Purpose: Extract coefficient tables and a residual diagnostic from the legacy
# single fitted-model object.
# Inputs: results/model_fitted.rds written by scripts/02_run_glmm.R.
# Outputs: Fixed/random-effect CSVs and residual diagnostic PNG/SVG files under
# the generic results subdirectories.
# Dependencies: lme4, ggplot2, here, and svglite.
# Usage: Rscript scripts/03_results.R
# Workflow: Final stage of the legacy 01 -> 02 -> 03 pipeline.

library(lme4)
library(ggplot2)
library(here)

# --- Load model ---
model <- readRDS(here("results", "model_fitted.rds"))

# --- Fixed effects table ---
fixef_table <- as.data.frame(summary(model)$coefficients)
fixef_table$Term <- rownames(fixef_table)
fixef_table <- fixef_table[, c("Term", setdiff(names(fixef_table), "Term"))]

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(fixef_table, here("results", "tables", "fixed_effects.csv"), row.names = FALSE)

# --- Random effects ---
ranef_df <- as.data.frame(ranef(model, condVar = TRUE))
write.csv(ranef_df, here("results", "tables", "random_effects.csv"), row.names = FALSE)

# --- Diagnostic plot: residuals vs fitted ---
diagnostics_df <- data.frame(
  .fitted = fitted(model),
  .residuals = residuals(model, type = "pearson")
)

p <- ggplot(diagnostics_df, aes(x = .fitted, y = .residuals)) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Fitted",
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_minimal()

dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("results", "figures", "residuals_vs_fitted.png"), p, width = 7, height = 5)
ggsave(
  here("results", "figures", "residuals_vs_fitted.svg"),
  p,
  width = 7,
  height = 5,
  device = svglite::svglite
)

cat("Results saved to results/tables/ and results/figures/\n")
