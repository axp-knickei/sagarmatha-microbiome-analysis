#!/usr/bin/env Rscript
# Purpose: Visualize sample coverage, LOSO predictions, baseline comparisons,
# sensitivity analyses, and feature eligibility from a completed validation run.
# Inputs: Audit, prediction, comparison, and feature CSVs in RF_VALIDATION_RUN_DIR.
# Outputs: Matching PNG, PDF, and SVG figures in run-specific figure directories.
# Dependencies: dplyr, tidyr, ggplot2, and R/common.R.
# Usage: Normally invoked by run_random_forest_subject_validation.R.
# Workflow: Stage 10; creates figures without refitting models.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
source("analysis/random_forest_subject_validation/R/common.R")
cfg <- read_validation_config(); run_dir <- Sys.getenv("RF_VALIDATION_RUN_DIR", ""); if (!nzchar(run_dir)) stop("Set RF_VALIDATION_RUN_DIR.")
save_three <- function(p, stem, width = 7, height = 5) {
  ggsave(file.path(run_dir, "figures/png", paste0(stem, ".png")), p, width = width, height = height, dpi = 180, bg = "white")
  ggsave(file.path(run_dir, "figures/pdf", paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(run_dir, "figures/svg", paste0(stem, ".svg")), p, width = width, height = height, bg = "white")
}
theme_set(theme_minimal(base_size = 11))
a <- read.csv(file.path(run_dir, "audit", "samples_by_altitude.csv")); save_three(ggplot(a, aes(factor(altitude), n_samples)) + geom_col(fill="#376795") + labs(x="Altitude (m)", y="Samples", title="Sample distribution by altitude"), "01_sample_distribution_altitude")
cvg <- read.csv(file.path(run_dir, "audit", "subject_altitude_coverage.csv")); save_three(ggplot(cvg, aes(factor(altitude), Subject, fill=n_samples)) + geom_tile(color="white") + scale_fill_gradient(low="white", high="#B54A35") + labs(x="Altitude (m)", y="Subject", title="Subject-by-altitude coverage"), "02_subject_altitude_coverage")
sp <- read.csv(file.path(run_dir, "audit", "existing_split_audit.csv")); save_three(ggplot(sp, aes(factor(altitude), n_samples, fill=split)) + geom_col(position="dodge") + labs(x="Altitude (m)", y="Samples", title="Legacy train/test composition"), "03_existing_split_composition")
cm <- read.csv(file.path(run_dir, "predictions", "loso_confusion_matrix.csv")); save_three(ggplot(cm, aes(predicted, observed, fill=Freq)) + geom_tile(color="white") + geom_text(aes(label=Freq)) + labs(title="LOSO confusion matrix"), "04_loso_confusion_matrix")
p <- read.csv(file.path(run_dir, "predictions", "loso_predictions.csv"))
levels <- unlist(cfg$data$altitude_levels); pc <- classification_metrics(p$observed, p$predicted, levels)$per_class |> pivot_longer(c(precision,recall), names_to="metric", values_to="value")
save_three(ggplot(pc, aes(class,value,color=metric,group=metric))+geom_point()+geom_line()+coord_cartesian(ylim=c(0,1))+labs(x="Altitude (m)",title="LOSO per-class precision and recall"), "05_per_class_precision_recall")
fm <- read.csv(file.path(run_dir, "predictions", "loso_fold_metrics.csv")); save_three(ggplot(fm,aes(reorder(heldout_subject,balanced_accuracy),balanced_accuracy))+geom_col(fill="#376795")+coord_flip()+coord_cartesian(ylim=c(0,1))+labs(x="Held-out subject",title="LOSO balanced accuracy by subject"), "06_fold_balanced_accuracy")
b <- read.csv(file.path(run_dir,"comparisons","baseline_model_comparison.csv"))
loso <- read.csv(file.path(run_dir,"predictions","loso_overall_metrics.csv"))
b <- bind_rows(
  b |> transmute(model, balanced_accuracy = balanced_accuracy_mean, macro_f1 = macro_f1_mean),
  tibble(model=paste0("random_forest_", cfg$analysis$mode),balanced_accuracy=loso$balanced_accuracy,macro_f1=loso$macro_f1)
) |> pivot_longer(c(balanced_accuracy,macro_f1),names_to="metric",values_to="value")
save_three(ggplot(b,aes(model,value,fill=metric))+geom_col(position="dodge")+coord_cartesian(ylim=c(0,1))+labs(x=NULL,title="Comparison with baselines"), "07_baseline_comparison")
for (item in list(c("class_balance_strategy_comparison.csv","08_class_balance_comparison","strategy"), c("feature_strategy_comparison.csv","09_feature_strategy_comparison","feature_strategy"), c("tree_count_stability.csv","10_tree_count_stability","trees"))) {
  path <- file.path(run_dir,"comparisons",item[[1]]); if (file.exists(path)) {z<-read.csv(path); save_three(ggplot(z,aes(x=.data[[item[[3]]]],y=balanced_accuracy,group=1))+geom_point()+geom_line()+labs(x=item[[3]],title=item[[2]]),item[[2]])}
}
fs <- read.csv(file.path(run_dir,"features","feature_selection_frequency.csv")) |>
  arrange(desc(eligible_folds), desc(selection_frequency_when_eligible)) |>
  slice_head(n=20)
save_three(
  ggplot(fs,aes(reorder(feature,eligible_folds),eligible_folds))+
    geom_col(fill="#376795")+coord_flip()+
    labs(x=NULL,y="Eligible LOSO folds",title="Training-fold feature eligibility"),
  "11_feature_selection_frequency",8,7
)
