# ============================================
# recommend_glmm_family.R
# ============================================
# Applies a simple GLMM-family decision table to the MAG response variable and
# writes both a machine-readable CSV and a readable report.

parse_args <- function(args) {
  parsed <- list(
    input = "data/processed/cleaned_data.rds",
    response = "Y",
    table_output = "results/tables/glmm_family_recommendation.csv",
    report_output = "results/logs/glmm_family_recommendation.txt"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/recommend_glmm_family.R [options]\n\n",
        "Options:\n",
        "  --input PATH          Processed RDS input. Default: data/processed/cleaned_data.rds\n",
        "  --response COLUMN     Response column. Default: Y\n",
        "  --table-output PATH   CSV output. Default: results/tables/glmm_family_recommendation.csv\n",
        "  --report-output PATH  Text report output. Default: results/logs/glmm_family_recommendation.txt\n",
        "  --help                Show this help.\n",
        sep = ""
      )
      quit(status = 0)
    } else if (arg %in% c("--input", "--response", "--table-output", "--report-output")) {
      if (i == length(args)) {
        stop("Missing value for argument: ", arg, call. = FALSE)
      }
      value <- args[[i + 1]]
      if (arg == "--input") {
        parsed$input <- value
      } else if (arg == "--response") {
        parsed$response <- value
      } else if (arg == "--table-output") {
        parsed$table_output <- value
      } else if (arg == "--report-output") {
        parsed$report_output <- value
      }
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  parsed
}

fmt <- function(x, digits = 6) {
  if (is.na(x)) {
    return(NA_character_)
  }
  format(x, scientific = TRUE, digits = digits)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input, call. = FALSE)
}

df <- readRDS(args$input)
if (!args$response %in% names(df)) {
  stop("Response column not found: ", args$response, call. = FALSE)
}

y <- df[[args$response]]
if (!is.numeric(y)) {
  stop("Response column must be numeric: ", args$response, call. = FALSE)
}

y <- y[!is.na(y)]
n <- length(y)
n_zero <- sum(y == 0)
n_one <- sum(y == 1)
n_positive <- sum(y > 0)
n_strict_beta <- sum(y > 0 & y < 1)
n_binary_only <- sum(y %in% c(0, 1))
n_fractional_non_binary <- sum(!(y %in% c(0, 1)))
n_integer_like <- sum(abs(y - round(y)) < .Machine$double.eps^0.5)

y_min <- min(y, na.rm = TRUE)
y_max <- max(y, na.rm = TRUE)
y_mean <- mean(y, na.rm = TRUE)
y_variance <- stats::var(y, na.rm = TRUE)
zero_fraction <- n_zero / n
one_fraction <- n_one / n
strict_beta_fraction <- n_strict_beta / n
variance_mean_ratio <- if (y_mean > 0) y_variance / y_mean else NA_real_

is_binary <- n > 0 && n_binary_only == n
is_count <- n > 0 && all(y >= 0) && n_integer_like == n
is_bounded_0_1 <- n > 0 && y_min >= 0 && y_max <= 1
is_strictly_between_0_1 <- n > 0 && y_min > 0 && y_max < 1
is_fractional_proportion <- is_bounded_0_1 && n_fractional_non_binary > 0 && !is_count
has_many_zeros <- zero_fraction >= 0.10

family_rows <- data.frame(
  family = c("Gaussian", "Poisson", "Negative binomial", "Binomial", "Beta"),
  for_what_data = c(
    "Continuous data; basically LMM-like",
    "Counts, assumes mean approximately variance",
    "Counts with overdispersion",
    "Binary or success/failure data",
    "Proportions between 0 and 1, excluding exact 0 and 1"
  ),
  suitable_as_primary_for_full_raw_y = c(
    is_fractional_proportion && !has_many_zeros,
    is_count && !is_binary && !is.na(variance_mean_ratio) && variance_mean_ratio <= 2,
    is_count && !is_binary && !is.na(variance_mean_ratio) && variance_mean_ratio > 2,
    is_binary,
    is_strictly_between_0_1
  ),
  notes = c(
    if (is_fractional_proportion) {
      "Raw Y is continuous/fractional, but many zeros and skew make raw Gaussian unattractive; transformed LMM is safer as a screen."
    } else {
      "Raw Y is not an ordinary unbounded continuous response."
    },
    if (is_count) {
      paste0("Count-like check passed; variance/mean ratio = ", fmt(variance_mean_ratio), ".")
    } else {
      "Not suitable: Y has fractional non-integer values, not counts."
    },
    if (is_count) {
      paste0("Count-like check passed; variance/mean ratio = ", fmt(variance_mean_ratio), ".")
    } else {
      "Not suitable: Y has fractional non-integer values, not counts."
    },
    if (is_binary) {
      "Suitable only if modeling detection as Y > 0 or true 0/1 outcome."
    } else {
      "Not suitable for abundance: Y is not binary 0/1. It can be used for a separate detection model, detected = Y > 0."
    },
    if (is_strictly_between_0_1) {
      "Suitable for raw Y under the table because all values are strictly between 0 and 1."
    } else if (is_bounded_0_1 && n_zero > 0) {
      "Conceptually closest GLMM family for proportions, but not directly valid because Y contains exact zeros."
    } else {
      "Not suitable because Y is not strictly bounded inside (0, 1)."
    }
  ),
  stringsAsFactors = FALSE
)

direct_family <- family_rows$family[family_rows$suitable_as_primary_for_full_raw_y]
direct_family_text <- if (length(direct_family) > 0) paste(direct_family, collapse = "; ") else "None of the listed GLMM families directly fits the full raw Y values"

recommended_primary <- if (is_fractional_proportion && has_many_zeros) {
  "Use transformed Gaussian/LMM for the current abundance screen: log10(Y + pseudocount) ~ time + (1 | Subject)."
} else if (is_strictly_between_0_1) {
  "Use Beta GLMM for raw proportional abundance."
} else if (is_binary) {
  "Use Binomial/Logistic GLMM."
} else if (is_count && !is.na(variance_mean_ratio) && variance_mean_ratio > 2) {
  "Use Negative binomial GLMM."
} else if (is_count) {
  "Use Poisson GLMM."
} else {
  "Manual model review required."
}

recommended_followup <- if (is_fractional_proportion && has_many_zeros) {
  "For a more rigorous GLMM-style follow-up, use a two-part model: Binomial GLMM for detected = Y > 0, then Beta/Gamma/lognormal-style model on positive abundance only; or use a zero-inflated beta/hurdle microbiome model."
} else {
  "No special follow-up recommendation beyond model diagnostics."
}

summary_row <- data.frame(
  response_column = args$response,
  n = n,
  min = y_min,
  max = y_max,
  max_as_percent = y_max * 100,
  mean = y_mean,
  variance = y_variance,
  variance_mean_ratio = variance_mean_ratio,
  n_zero = n_zero,
  zero_fraction = zero_fraction,
  n_one = n_one,
  one_fraction = one_fraction,
  n_positive = n_positive,
  n_strictly_between_0_1 = n_strict_beta,
  strict_between_0_1_fraction = strict_beta_fraction,
  n_fractional_non_binary = n_fractional_non_binary,
  is_binary = is_binary,
  is_count = is_count,
  is_bounded_0_1 = is_bounded_0_1,
  is_strictly_between_0_1 = is_strictly_between_0_1,
  has_many_zeros = has_many_zeros,
  direct_family_from_table = direct_family_text,
  recommended_primary = recommended_primary,
  recommended_followup = recommended_followup,
  stringsAsFactors = FALSE
)

output <- merge(
  family_rows,
  summary_row,
  by = NULL
)

dir.create(dirname(args$table_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report_output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(output, args$table_output, row.names = FALSE)

report_connection <- file(args$report_output, open = "wt")
sink(report_connection, split = TRUE)

cat("GLMM Family Recommendation\n")
cat("==========================\n\n")
cat("Input file:", args$input, "\n")
cat("Response column:", args$response, "\n\n")

cat("Observed response summary\n")
cat("-------------------------\n")
cat("n:", n, "\n")
cat("Range:", fmt(y_min), "to", fmt(y_max), "\n")
cat("Maximum as percent:", sprintf("%.4f%%", y_max * 100), "\n")
cat("Mean:", fmt(y_mean), "\n")
cat("Variance:", fmt(y_variance), "\n")
cat("Variance / mean:", fmt(variance_mean_ratio), "\n")
cat("Zeros:", n_zero, sprintf("(%.2f%%)", 100 * zero_fraction), "\n")
cat("Ones:", n_one, sprintf("(%.2f%%)", 100 * one_fraction), "\n")
cat("Strictly inside (0, 1):", n_strict_beta, sprintf("(%.2f%%)", 100 * strict_beta_fraction), "\n")
cat("Fractional non-binary values:", n_fractional_non_binary, "\n\n")

cat("Family table decision\n")
cat("---------------------\n")
for (i in seq_len(nrow(family_rows))) {
  cat(family_rows$family[[i]], "\n", sep = "")
  cat("  Suitable as primary family for full raw Y:", family_rows$suitable_as_primary_for_full_raw_y[[i]], "\n")
  cat("  Note:", family_rows$notes[[i]], "\n\n")
}

cat("Direct family from listed table:\n")
cat(direct_family_text, "\n\n")

cat("Recommended primary analysis:\n")
cat(recommended_primary, "\n\n")

cat("Recommended follow-up:\n")
cat(recommended_followup, "\n\n")

cat("Bottom line:\n")
cat(
  "The raw MAG abundance is proportional relative abundance on a 0-1 scale, ",
  "but it contains many exact zeros. Beta is the closest listed GLMM family for ",
  "proportional data, but standard beta GLMM excludes exact 0 and 1, so it is ",
  "not directly valid for the full raw response. The current transformed ",
  "Gaussian/LMM screen is therefore the most practical primary option from this table, ",
  "with a two-part or zero-inflated beta-style model as a stronger follow-up.\n",
  sep = ""
)

sink()
close(report_connection)

cat("Saved GLMM family table to ", args$table_output, "\n", sep = "")
cat("Saved GLMM family report to ", args$report_output, "\n", sep = "")
