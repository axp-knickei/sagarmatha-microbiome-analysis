project_root <- Sys.getenv("NETMOSS_PROJECT_ROOT")

test_that("NetMoss PBS wrappers are valid OpenPBS shell scripts", {
  scripts <- file.path(project_root, c(
    "hpc/update_netmoss_env.pbs",
    "hpc/job_netmoss_readiness.pbs",
    "hpc/job_netmoss_analysis.pbs",
    "hpc/job_netmoss_altitude_analysis.pbs",
    "hpc/job_netmoss_altitude_stability.pbs",
    "hpc/job_netmoss_altitude_stability_smoke_test.pbs",
    "hpc/job_netmoss_publication.pbs",
    "hpc/job_netmoss_smoke_test.pbs"
  ))
  expect_true(all(file.exists(scripts)))
  for (script in scripts) {
    lines <- readLines(script, warn = FALSE)
    expect_true(any(grepl("^#PBS -N ", lines)))
    expect_true(any(grepl("^#PBS -l ", lines)))
    expect_true(any(grepl("PBS_O_WORKDIR", lines, fixed = TRUE)))
    expect_false(any(grepl("\\b(sbatch|srun|squeue)\\b", lines)))
    status <- system2("bash", c("-n", script), stdout = FALSE, stderr = FALSE)
    expect_equal(status, 0L)
  }
})

test_that("full analysis job requires explicit exploratory acknowledgement", {
  lines <- readLines(file.path(project_root, "hpc/job_netmoss_analysis.pbs"), warn = FALSE)
  expect_true(any(grepl("NETMOSS_ACKNOWLEDGE_EXPLORATORY", lines, fixed = TRUE)))
  expect_true(any(grepl("--acknowledge-exploratory", lines, fixed = TRUE)))
})

test_that("altitude analysis job uses its separate acknowledged scenario", {
  lines <- readLines(file.path(project_root, "hpc/job_netmoss_altitude_analysis.pbs"), warn = FALSE)
  expect_true(any(grepl("NETMOSS_ALTITUDE_ACKNOWLEDGE_EXPLORATORY", lines, fixed = TRUE)))
  expect_true(any(grepl("config/netmoss_altitude_analysis.yml", lines, fixed = TRUE)))
  expect_true(any(grepl("--acknowledge-exploratory", lines, fixed = TRUE)))
})

test_that("altitude stability job is isolated and explicitly acknowledged", {
  lines <- readLines(file.path(project_root, "hpc/job_netmoss_altitude_stability.pbs"), warn = FALSE)
  expect_true(any(grepl("NETMOSS_STABILITY_ACKNOWLEDGE_EXPLORATORY", lines, fixed = TRUE)))
  expect_true(any(grepl("config/netmoss_altitude_stability.yml", lines, fixed = TRUE)))
  expect_true(any(grepl("--acknowledge-exploratory", lines, fixed = TRUE)))
})
