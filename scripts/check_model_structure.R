# ============================================
# check_model_structure.R - Singularity Check
# ============================================
# Inspects the fitted model structure without changing the configured model.

library(lme4)
library(here)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

model_path <- here("results", "model_fitted.rds")
report_path <- here("results", "logs", "singularity_check.txt")

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

sink(report_path, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Model Singularity Check\n")
cat("=======================\n\n")
cat("Model file:", model_path, "\n")
cat("Report file:", report_path, "\n\n")

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path, call. = FALSE)
}

model <- readRDS(model_path)

singular <- lme4::isSingular(model)
varcorr <- lme4::VarCorr(model)
varcorr_df <- as.data.frame(varcorr)

cat("lme4::isSingular(model):", singular, "\n\n")

cat("lme4::VarCorr(model):\n")
print(varcorr)
cat("\n")

cat("Variance components table:\n")
print(varcorr_df)
cat("\n")

subject_intercept <- subset(
  varcorr_df,
  grp == "Subject" & var1 == "(Intercept)" & is.na(var2)
)

if (nrow(subject_intercept) > 0) {
  subject_variance <- subject_intercept$vcov[[1]]
  cat("Subject random-intercept variance:", subject_variance, "\n")
  cat("Subject random-intercept standard deviation:", subject_intercept$sdcor[[1]], "\n")

  if (isTRUE(all.equal(subject_variance, 0))) {
    cat("Interpretation: Subject random-intercept variance is exactly zero.\n")
  } else if (subject_variance < 1e-8) {
    cat("Interpretation: Subject random-intercept variance is near zero.\n")
  } else {
    cat("Interpretation: Subject random-intercept variance is not near zero.\n")
  }
} else {
  cat("Subject random-intercept variance was not found in VarCorr(model).\n")
}

cat("\nCandidate comparison models to consider later, without switching automatically:\n")
cat("1. Y ~ days_from_baseline_subject\n")
cat("2. Y ~ days_from_baseline_subject + (1 | Subject)\n")

cat("\nRecommendation:\n")
cat("Do not treat this singular model as final. Review whether the random intercept is needed, whether modeling should be per MAG, whether Y needs transformation, and whether a simpler fixed-effect model should be compared before HPC execution.\n")

cat("\nSaved singularity check to results/logs/singularity_check.txt\n")
