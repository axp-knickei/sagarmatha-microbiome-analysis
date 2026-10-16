# ============================================
# 02_run_glmm.R - Fit Configured Mixed Model
# ============================================
# Legacy single-model entrypoint. With family = "gaussian" this runs an LMM
# through lme4::lmer(); non-Gaussian families use lme4::glmer().
#
# Inputs: config/config.yml and its configured processed RDS file.
# Outputs: results/model_fitted.rds and results/logs/session_info.txt.
# Dependencies: lme4, dplyr, here, sessioninfo, and yaml.
# Usage: Rscript scripts/02_run_glmm.R
# Workflow: Legacy stage 02 after 01_data_prep.R; not the primary per-MAG fit.

library(lme4)
library(dplyr)
library(here)
library(sessioninfo)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

config <- yaml::read_yaml(here("config", "config.yml"))

# --- Load processed data ---
df <- readRDS(here(config$data$processed))

# --- Fit model ---
formula <- as.formula(config$model$formula)
family <- tolower(config$model$family)
optimizer <- config$model$optimizer
glm_families <- c("binomial", "poisson", "gamma", "inverse.gaussian")

cat("Fitting model:", deparse(formula), "\n")
cat("Family:", family, "\n")
cat("Optimizer:", optimizer, "\n\n")

if (family == "gaussian") {
  cat("Model type: Gaussian LMM via lme4::lmer()\n\n")
  model <- lmer(
    formula = formula,
    data = df,
    control = lmerControl(optimizer = optimizer)
  )
} else {
  if (!family %in% glm_families) {
    stop(
      "Unsupported model family: ", family,
      ". Supported families are: gaussian, ",
      paste(glm_families, collapse = ", "),
      call. = FALSE
    )
  }

  family_function <- switch(
    family,
    binomial = binomial,
    poisson = poisson,
    gamma = Gamma,
    inverse.gaussian = inverse.gaussian
  )
  cat("Model type: non-Gaussian GLMM via lme4::glmer()\n\n")
  model <- glmer(
    formula = formula,
    data = df,
    family = family_function(),
    control = glmerControl(optimizer = optimizer)
  )
}

# --- Save model ---
results_dir <- here("results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(model, file.path(results_dir, "model_fitted.rds"))

# --- Print summary ---
print(summary(model), correlation = FALSE)

# --- Save session info ---
log_dir <- here("results", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
sink(here(log_dir, "session_info.txt"))
print(session_info())
sink()

cat("\nModel saved to: results/model_fitted.rds\n")
cat("Session info saved to:", here(log_dir, "session_info.txt"), "\n")
