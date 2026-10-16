restricted_permutations <- function(metadata, block, nperm) {
  control <- permute::how(nperm = nperm)
  permute::setBlocks(control) <- factor(metadata[[block]])
  control
}

tidy_adonis <- function(result, metric, model, scheme, permutations, n) {
  out <- as.data.frame(result)
  out$hypothesis <- rownames(out)
  rownames(out) <- NULL
  names(out) <- sub("Pr\\(>F\\)", "p_value", names(out))
  data.frame(metric = metric, preprocessing_branch = metric, model = model,
             hypothesis = out$hypothesis, permutation_scheme = scheme,
             permutations = permutations, df = out$Df, sum_squares = out$SumOfSqs,
             r_squared = out$R2, statistic = out$F, p_value = out$p_value,
             sample_count = n, stringsAsFactors = FALSE)
}
