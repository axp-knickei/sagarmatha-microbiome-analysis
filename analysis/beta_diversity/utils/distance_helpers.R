total_sum_scale <- function(x) x / rowSums(x)
hellinger_transform <- function(x) sqrt(total_sum_scale(x))

multiplicative_replace <- function(x, fraction = 0.65) {
  x <- total_sum_scale(x)
  replaced <- t(apply(x, 1, function(row) {
    zeros <- row == 0; k <- sum(zeros)
    if (!k) return(row)
    positive <- row[!zeros]
    delta <- min(fraction * min(positive), 0.5 / k)
    row[zeros] <- delta
    row[!zeros] <- positive * ((1 - k * delta) / sum(positive))
    row
  }))
  dimnames(replaced) <- dimnames(x)
  replaced
}

clr_transform <- function(x) {
  if (any(x <= 0)) stop("CLR requires strictly positive values after zero replacement.", call. = FALSE)
  logged <- log(x)
  logged - rowMeans(logged)
}

quantitative_jaccard_ruzička <- function(x) {
  # vegan's quantitative Jaccard is 2B/(1+B), algebraically equal to Ružička dissimilarity.
  vegan::vegdist(x, method = "jaccard", binary = FALSE)
}

distance_summary <- function(distance) {
  values <- as.numeric(distance)
  data.frame(minimum = min(values), median = stats::median(values), mean = mean(values),
             maximum = max(values), missing_count = sum(is.na(values)))
}

distance_long <- function(distance) {
  m <- as.matrix(distance); idx <- which(lower.tri(m), arr.ind = TRUE)
  data.frame(sample_id_1 = rownames(m)[idx[, 1]], sample_id_2 = colnames(m)[idx[, 2]], distance = m[idx])
}
