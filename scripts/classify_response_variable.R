# ============================================
# classify_response_variable.R
# ============================================
# Classifies the MAG response variable using the project decision tree and
# writes a compact summary table plus a readable decision report.

parse_args <- function(args) {
  parsed <- list(
    input = "data/processed/cleaned_data.rds",
    response = "Y",
    table_output = "results/tables/response_variable_classification.csv",
    report_output = "results/logs/response_variable_classification.txt"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      cat(
        "Usage: Rscript scripts/classify_response_variable.R [options]\n\n",
        "Options:\n",
        "  --input PATH          Processed RDS input. Default: data/processed/cleaned_data.rds\n",
        "  --response COLUMN     Response column. Default: Y\n",
        "  --table-output PATH   CSV summary output. Default: results/tables/response_variable_classification.csv\n",
        "  --report-output PATH  Text report output. Default: results/logs/response_variable_classification.txt\n",
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

format_number <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  format(x, scientific = TRUE, digits = 6)
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

y_non_missing <- y[!is.na(y)]
y_positive <- y_non_missing[y_non_missing > 0]

n_total <- length(y)
n_missing <- sum(is.na(y))
n_non_missing <- length(y_non_missing)
n_zero <- sum(y_non_missing == 0)
n_one <- sum(y_non_missing == 1)
n_binary_values <- sum(y_non_missing %in% c(0, 1))
n_other_than_0_1 <- sum(!(y_non_missing %in% c(0, 1)))
n_positive <- length(y_positive)
n_unique <- length(unique(y_non_missing))
n_unique_positive <- length(unique(y_positive))

y_min <- min(y_non_missing, na.rm = TRUE)
y_max <- max(y_non_missing, na.rm = TRUE)
y_mean <- mean(y_non_missing, na.rm = TRUE)
y_variance <- stats::var(y_non_missing, na.rm = TRUE)
y_median <- stats::median(y_non_missing, na.rm = TRUE)
zero_fraction <- if (n_non_missing > 0) n_zero / n_non_missing else NA_real_
positive_median <- if (n_positive > 0) stats::median(y_positive, na.rm = TRUE) else NA_real_

is_binary <- n_non_missing > 0 && all(y_non_missing %in% c(0, 1))
is_integer_count <- n_non_missing > 0 && all(y_non_missing >= 0) && all(abs(y_non_missing - round(y_non_missing)) < .Machine$double.eps^0.5)
is_bounded_0_1 <- n_non_missing > 0 && y_min >= 0 && y_max <= 1
is_bounded_0_100 <- n_non_missing > 0 && y_min >= 0 && y_max <= 100 && y_max > 1
has_many_zeros <- !is.na(zero_fraction) && zero_fraction >= 0.10
has_continuous_fractional_values <- n_other_than_0_1 > 0 && !is_integer_count

decision_path <- if (is_binary) {
  "Binary yes/no -> Logistic GLMM"
} else if (is_integer_count) {
  if (!is.na(y_variance) && y_mean > 0 && y_variance > 2 * y_mean) {
    "Count -> Variance much larger than mean -> Negative binomial GLMM"
  } else {
    "Count -> Variance similar to mean -> Poisson GLMM"
  }
} else if (is_bounded_0_1 && has_many_zeros) {
  "Proportion between 0 and 1 -> Relative abundance with many zeros -> Transformation + LMM, or specialized microbiome model"
} else if (is_bounded_0_1) {
  "Proportion between 0 and 1 -> Consider beta/binomial-family model depending on denominator information"
} else if (has_continuous_fractional_values) {
  "Continuous value -> Assess transformation normality; if very skewed/many zeros, consider GLMM or two-part model"
} else {
  "Response type requires manual review"
}

recommended_primary <- if (is_bounded_0_1 && has_many_zeros && has_continuous_fractional_values) {
  "log10(Y + pseudocount) LMM as practical primary screen; consider two-part or specialized microbiome model as follow-up"
} else if (is_binary) {
  "Logistic GLMM"
} else if (is_integer_count && !is.na(y_variance) && y_mean > 0 && y_variance > 2 * y_mean) {
  "Negative binomial GLMM"
} else if (is_integer_count) {
  "Poisson GLMM"
} else {
  "Manual review"
}

summary <- data.frame(
  response_column = args$response,
  n_total = n_total,
  n_missing = n_missing,
  n_non_missing = n_non_missing,
  n_unique = n_unique,
  min = y_min,
  max = y_max,
  max_as_percent = y_max * 100,
  mean = y_mean,
  variance = y_variance,
  median = y_median,
  n_zero = n_zero,
  zero_fraction = zero_fraction,
  n_one = n_one,
  n_positive = n_positive,
  n_unique_positive = n_unique_positive,
  positive_min = if (n_positive > 0) min(y_positive, na.rm = TRUE) else NA_real_,
  positive_median = positive_median,
  positive_max = if (n_positive > 0) max(y_positive, na.rm = TRUE) else NA_real_,
  n_binary_values = n_binary_values,
  n_other_than_0_1 = n_other_than_0_1,
  is_binary_0_1 = is_binary,
  is_integer_count = is_integer_count,
  is_bounded_0_1 = is_bounded_0_1,
  is_bounded_0_100_percent = is_bounded_0_100,
  has_many_zeros = has_many_zeros,
  has_continuous_fractional_values = has_continuous_fractional_values,
  decision_path = decision_path,
  recommended_primary_model = recommended_primary,
  stringsAsFactors = FALSE
)

dir.create(dirname(args$table_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report_output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, args$table_output, row.names = FALSE)

report_connection <- file(args$report_output, open = "wt")
sink(report_connection, split = TRUE)

cat("Response Variable Classification\n")
cat("================================\n\n")
cat("Input file:", args$input, "\n")
cat("Response column:", args$response, "\n\n")

cat("Observed scale\n")
cat("--------------\n")
cat("Total observations:", n_total, "\n")
cat("Non-missing observations:", n_non_missing, "\n")
cat("Unique values:", n_unique, "\n")
cat("Range:", format_number(y_min), "to", format_number(y_max), "\n")
cat("Maximum as percent:", sprintf("%.4f%%", y_max * 100), "\n")
cat("Zeros:", n_zero, sprintf("(%.2f%%)", 100 * zero_fraction), "\n")
cat("Ones:", n_one, "\n")
cat("Values other than 0/1:", n_other_than_0_1, "\n")
cat("Positive values:", n_positive, "\n")
cat("Unique positive values:", n_unique_positive, "\n")
cat("Positive range:", format_number(summary$positive_min), "to", format_number(summary$positive_max), "\n\n")

cat("Decision-tree classification\n")
cat("----------------------------\n")
cat("Binary 0/1:", is_binary, "\n")
cat("Integer count:", is_integer_count, "\n")
cat("Bounded 0-1:", is_bounded_0_1, "\n")
cat("Bounded 0-100 percent scale:", is_bounded_0_100, "\n")
cat("Many zeros:", has_many_zeros, "\n")
cat("Continuous fractional values:", has_continuous_fractional_values, "\n\n")

cat("Decision path:\n")
cat(decision_path, "\n\n")
cat("Recommended modeling interpretation:\n")
cat(recommended_primary, "\n\n")

cat("Plain-language conclusion:\n")
cat(
  "The response variable is continuous relative abundance stored as a 0-1 fraction, ",
  "not 0-100 percent and not binary detection. Because many observations are zero, ",
  "the current log10(Y + pseudocount) LMM is a practical screening model, while a ",
  "two-part or specialized microbiome model would be a stronger follow-up analysis.\n",
  sep = ""
)

sink()
close(report_connection)

cat("Saved response classification table to ", args$table_output, "\n", sep = "")
cat("Saved response classification report to ", args$report_output, "\n", sep = "")
