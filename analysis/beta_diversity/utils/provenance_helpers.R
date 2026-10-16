timestamp_beta <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

beta_log <- function(config, stage, status, message) {
  ensure_beta_dirs(config)
  line <- paste(timestamp_beta(), stage, status, message, sep = "\t")
  cat(line, "\n")
  cat(line, "\n", file = beta_output_dir(config, "logs", paste0(stage, ".log")), append = TRUE)
}

with_stage_log <- function(config, stage, code) {
  start <- Sys.time()
  beta_log(config, stage, "STARTED", paste("config=", config$.config_path))
  tryCatch({
    value <- force(code)
    beta_log(config, stage, "COMPLETED", sprintf("elapsed_seconds=%.3f", as.numeric(difftime(Sys.time(), start, units = "secs"))))
    value
  }, error = function(e) {
    beta_log(config, stage, "FAILED", conditionMessage(e))
    stop(e)
  })
}

sha256_file <- function(path) {
  require_beta_packages("digest")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_yaml_manifest <- function(config, path, values) {
  require_beta_packages("yaml")
  manifest <- c(list(generated_at = timestamp_beta(), config = config$.config_path,
                     git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "unknown"),
                     seed = config$analysis$random_seed), values)
  yaml::write_yaml(manifest, path)
}

record_run_provenance <- function(config, args) {
  out <- beta_output_dir(config, "provenance")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  inputs <- unlist(config$inputs[c("abundance_file", "metadata_file", "taxonomy_file", "tree_file")], use.names = TRUE)
  inputs <- inputs[!is.na(inputs) & nzchar(inputs)]
  checksums <- data.frame(role = names(inputs), path = inputs,
                          sha256 = vapply(inputs, sha256_file, character(1)), stringsAsFactors = FALSE)
  write_tsv(checksums, file.path(out, "input_checksums.tsv"))
  writeLines(c(paste("timestamp:", timestamp_beta()), paste("command:", paste(commandArgs(), collapse = " ")),
               paste("seed:", config$analysis$random_seed)), file.path(out, "invocation.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
  versions <- data.frame(package = c("R", "vegan", "ape", "permute", "ggplot2", "yaml", "digest", "lme4"),
                         version = c(as.character(getRversion()), vapply(c("vegan", "ape", "permute", "ggplot2", "yaml", "digest", "lme4"), function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, character(1))))
  write_tsv(versions, file.path(out, "software_versions.tsv"))
}
