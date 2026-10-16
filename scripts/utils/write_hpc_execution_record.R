# Reusable, fail-safe HPC execution recorder.
# Source this file from analysis scripts and call write_hpc_execution_record()
# after primary outputs and session information have been written.

safe_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) NA_character_ else value
}

safe_system <- function(command, args = character()) {
  tryCatch(
    {
      value <- suppressWarnings(system2(command, args = args, stdout = TRUE, stderr = FALSE))
      status <- attr(value, "status")
      if (!is.null(status) && !identical(status, 0L)) {
        return(NA_character_)
      }
      value <- value[nzchar(value)]
      if (length(value) == 0) NA_character_ else paste(value, collapse = " | ")
    },
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
}

safe_info <- function(name) {
  value <- tryCatch(Sys.info()[[name]], error = function(e) NA_character_)
  if (is.null(value) || is.na(value) || !nzchar(value)) NA_character_ else as.character(value)
}

safe_package_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(NA_character_)
  }
  tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
}

collapse_field <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  x <- as.character(x)
  x <- x[nzchar(x) & !is.na(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

format_time_field <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(format(x, "%Y-%m-%d %H:%M:%S %Z"))
  }
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) {
    return(NA_character_)
  }
  as.character(x)
}

detect_scheduler <- function() {
  if (!is.na(safe_env("PBS_JOBID"))) {
    return("OpenPBS")
  }
  if (!is.na(safe_env("SLURM_JOB_ID"))) {
    return("SLURM")
  }
  "local_or_unknown"
}

detect_scheduler_file <- function(kind) {
  job_id <- safe_env("PBS_JOBID")
  job_name <- safe_env("PBS_JOBNAME")
  if (is.na(job_id) || is.na(job_name)) {
    return(NA_character_)
  }

  candidates <- c(
    file.path("logs", paste0(job_name, ".", job_id, ".", kind)),
    file.path("logs", paste0(job_name, ".", job_id, ".nebula.", kind))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) existing[[1]] else NA_character_
}

markdown_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("\\|", "\\\\|", x)
}

write_record_markdown <- function(record, path) {
  lines <- c(
    "# HPC Execution Record",
    "",
    paste("Written:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "field | value",
    "--- | ---"
  )
  for (name in names(record)) {
    lines <- c(lines, paste(markdown_escape(name), markdown_escape(record[[name]]), sep = " | "))
  }
  writeLines(lines, path)
}

write_hpc_execution_record <- function(
  analysis_name,
  script_path,
  result_dir,
  input_data_file = NA,
  window_definition_file = NA,
  config_files = character(),
  number_of_tests = NA,
  fdr_family = NA,
  main_results_file = NA,
  summary_log_file = NA,
  session_info_file = NA,
  notes = NA,
  start_time = NA
) {
  tryCatch(
    {
      completion_time <- Sys.time()
      runtime_seconds <- NA_real_
      if (inherits(start_time, "POSIXt")) {
        runtime_seconds <- as.numeric(difftime(completion_time, start_time, units = "secs"))
      }

      key_packages <- c("lme4", "lmerTest", "tidyverse", "dplyr", "readr", "yaml", "broom.mixed", "ggplot2")
      package_versions <- stats::setNames(
        vapply(key_packages, safe_package_version, character(1)),
        key_packages
      )
      package_version_text <- paste(
        paste(names(package_versions), package_versions, sep = "="),
        collapse = "; "
      )

      conda_environment <- safe_env("CONDA_DEFAULT_ENV")
      if (is.na(conda_environment)) {
        conda_environment <- safe_env("CONDA_ENV_NAME")
      }

      conda_module <- safe_env("CONDA_MODULE_NAME")
      loaded_modules <- safe_env("LOADEDMODULES")
      if (is.na(conda_module)) {
        conda_module <- loaded_modules
      }

      record <- data.frame(
        project = basename(getwd()),
        analysis_name = analysis_name,
        script_path = script_path,
        script = script_path,
        scheduler = detect_scheduler(),
        job_id = safe_env("PBS_JOBID"),
        job_name = safe_env("PBS_JOBNAME"),
        pbs_o_workdir = safe_env("PBS_O_WORKDIR"),
        pbs_queue = safe_env("PBS_QUEUE"),
        pbs_nodefile = safe_env("PBS_NODEFILE"),
        hostname = safe_info("nodename"),
        hostname_command = safe_system("hostname"),
        user = safe_info("user"),
        working_directory = getwd(),
        git_branch = safe_system("git", c("rev-parse", "--abbrev-ref", "HEAD")),
        git_commit = safe_system("git", c("rev-parse", "HEAD")),
        git_status_short = safe_system("git", c("status", "--short")),
        start_time = format_time_field(start_time),
        completion_time = format_time_field(completion_time),
        runtime_seconds = runtime_seconds,
        command_args = collapse_field(commandArgs()),
        input_data_file = collapse_field(input_data_file),
        window_definition_file = collapse_field(window_definition_file),
        config_files = collapse_field(config_files),
        conda_module = conda_module,
        conda_environment = conda_environment,
        conda_prefix = safe_env("CONDA_PREFIX"),
        loaded_modules = loaded_modules,
        R_version = R.version.string,
        key_R_packages = package_version_text,
        number_of_tests = number_of_tests,
        FDR_family = fdr_family,
        fdr_family = fdr_family,
        result_directory = result_dir,
        main_results_file = collapse_field(main_results_file),
        summary_log_file = collapse_field(summary_log_file),
        session_info_file = collapse_field(session_info_file),
        scheduler_output_file = detect_scheduler_file("out"),
        scheduler_error_file = detect_scheduler_file("err"),
        notes = collapse_field(notes),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      reproducibility_dir <- file.path(result_dir, "reproducibility")
      dir.create(reproducibility_dir, recursive = TRUE, showWarnings = FALSE)
      csv_path <- file.path(reproducibility_dir, "hpc_execution_record.csv")
      md_path <- file.path(reproducibility_dir, "hpc_execution_record.md")

      utils::write.csv(record, csv_path, row.names = FALSE, na = "")
      write_record_markdown(record[1, , drop = TRUE], md_path)

      message("HPC execution record written to: ", csv_path)
      invisible(record)
    },
    error = function(e) {
      warning("Failed to write HPC execution record; continuing analysis: ", conditionMessage(e), call. = FALSE)
      invisible(NULL)
    }
  )
}
