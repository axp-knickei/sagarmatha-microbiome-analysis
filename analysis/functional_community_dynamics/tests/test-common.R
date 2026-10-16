source("analysis/functional_community_dynamics/R/common.R")

x <- matrix(c(1, 1, 2, 2, 2, 1), nrow = 2, byrow = TRUE)
stopifnot(all(abs(rowSums(close_rows(x)) - 1) < 1e-12))
stopifnot(all(is.finite(as.numeric(profile_distance(x, "hellinger")))))
stopifnot(all(is.finite(as.numeric(profile_distance(x, "bray_curtis")))))

zero_x <- matrix(c(1, 0, 2, 0, 3, 1), nrow = 2, byrow = TRUE)
replaced <- multiplicative_replace(zero_x)
stopifnot(all(replaced > 0), all(abs(rowSums(replaced) - 1) < 1e-10))
stopifnot(all(is.finite(as.numeric(profile_distance(zero_x, "aitchison")))))

no_change <- matrix(rep(c(0.2, 0.3, 0.5), 4), nrow = 4, byrow = TRUE)
stopifnot(max(as.numeric(profile_distance(no_change, "hellinger"))) == 0)

taxonomy_change <- rbind(c(0.8, 0.1, 0.1), c(0.1, 0.8, 0.1))
functional_stable <- rbind(c(0.5, 0.5), c(0.5, 0.5))
stopifnot(as.numeric(profile_distance(taxonomy_change, "hellinger")) > as.numeric(profile_distance(functional_stable, "hellinger")))

matched_change <- rbind(c(0.9, 0.1), c(0.1, 0.9))
stopifnot(abs(as.numeric(profile_distance(matched_change, "hellinger")) - as.numeric(profile_distance(matched_change, "hellinger"))) < 1e-12)

stopifnot(identical(safe_bh(c(0.01, 0.04)), stats::p.adjust(c(0.01, 0.04), method = "BH")))
permutation_metadata <- data.table(subject_id = rep(c("A", "B"), each = 2))
permutation_control <- restricted_permutations(permutation_metadata, 9)
stopifnot(identical(as.character(permute::getBlocks(permutation_control)), permutation_metadata$subject_id))

dispersion_only <- rbind(c(-1, 0), c(1, 0), c(-4, 0), c(4, 0))
dispersion_group <- factor(c("A", "A", "B", "B"))
stopifnot(all(abs(rowsum(dispersion_only, dispersion_group) / 2) < 1e-12))
stopifnot(mean(sqrt(rowSums(dispersion_only[dispersion_group == "B", , drop = FALSE]^2))) >
          mean(sqrt(rowSums(dispersion_only[dispersion_group == "A", , drop = FALSE]^2))))

recovery_meta <- data.table(
  sample_id = c("A_T1", "A_T5", "A_T10"), subject_id = "A",
  timepoint = factor(c("T1", "T5", "T10"), levels = c("T1", "T5", "T10"), ordered = TRUE),
  sample_date = as.IDate(c("2026-01-01", "2026-01-08", "2026-01-15"))
)
recovery_distance <- as.dist(matrix(c(0, 2, 1, 2, 0, 1, 1, 1, 0), 3, 3,
                                    dimnames = list(recovery_meta$sample_id, recovery_meta$sample_id)))
recovery_test <- recovery_records(recovery_distance, recovery_meta, "KO", "hellinger")
stopifnot(nrow(recovery_test) == 1L, abs(recovery_test$recovery_fraction_from_peak - 0.5) < 1e-12)

sensitivity_input <- rbindlist(lapply(sprintf("S%02d", 1:4), function(id) {
  tax <- c(0.2, 0.4, 0.6, 0.8)
  fun <- tax * 1.2
  data.table(
    subject_id = id,
    taxonomic_standardized_distance = tax,
    functional_standardized_distance = fun,
    log_functional_to_taxonomic_ratio = log(fun / tax)
  )
}))
sensitivity_test <- participant_transition_sensitivity(sensitivity_input, "KO", "primary")
stopifnot(
  nrow(sensitivity_test$participant) == 4L,
  nrow(sensitivity_test$summary) == 2L,
  all(sensitivity_test$summary$status == "fitted"),
  sensitivity_test$summary[
    sensitivity_family == "participant_spearman_transition_coupling",
    participant_summary
  ] == 1,
  abs(sensitivity_test$summary[
    sensitivity_family == "participant_median_relative_movement",
    ratio
  ] - 1.2) < 1e-12
)

cat("Functional community dynamics unit tests passed\n")
