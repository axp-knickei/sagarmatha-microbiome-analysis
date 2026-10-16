#!/bin/bash

# Static and executable smoke tests for the matched-timepoint companion jobs.

set -euo pipefail

output_root="${1:-/tmp/blackice_netmoss_matched_timepoint_smoke}"
mkdir -p "${output_root}"

if [ -n "${NETMOSS_RSCRIPT:-}" ]; then
  read -r -a rscript_cmd <<<"${NETMOSS_RSCRIPT}"
else
  rscript_cmd=(Rscript)
fi

pbs_scripts=(
  hpc/job_netmoss_matched_timepoint_companion.pbs
  hpc/job_netmoss_matched_timepoint_smoke_test.pbs
)

for script in "${pbs_scripts[@]}"; do
  test -s "${script}"
  bash -n "${script}"
  grep -q '^#PBS -N ' "${script}"
  grep -q '^#PBS -l ' "${script}"
  grep -q 'PBS_O_WORKDIR' "${script}"
  if grep -Eq '\b(sbatch|srun|squeue)\b' "${script}"; then
    echo "ERROR: Slurm command found in OpenPBS script: ${script}" >&2
    exit 1
  fi
done

"${rscript_cmd[@]}" -e '
  invisible(parse(file="scripts/netmoss/matched_timepoint_helpers.R"))
  invisible(parse(file="scripts/netmoss/run_matched_timepoint_companion.R"))
  invisible(parse(file="scripts/netmoss/summarize_matched_timepoint_companion.R"))
  cat("Matched-timepoint R parse checks passed\n")
'

"${rscript_cmd[@]}" scripts/netmoss/tests/testthat-matched-timepoint.R

"${rscript_cmd[@]}" scripts/netmoss/run_matched_timepoint_companion.R \
  --config config/netmoss_matched_timepoint_companion.yml \
  --output-dir "${output_root}" \
  --mode standard \
  --threads 2 \
  --test-mode \
  --force

"${rscript_cmd[@]}" scripts/netmoss/summarize_matched_timepoint_companion.R \
  --config config/netmoss_matched_timepoint_companion.yml \
  --output-dir "${output_root}" \
  --force

required_outputs=(
  "${output_root}/tables/timepoint_pair_support.csv"
  "${output_root}/tables/all_genus_pair_results.csv"
  "${output_root}/tables/all_participant_genus_changes.csv"
  "${output_root}/tables/predeclared_candidates.csv"
  "${output_root}/reproducibility/analysis_manifest.yml"
  "${output_root}/figures/matched_participant_support.pdf"
  "${output_root}/figures/top_ranked_evidence_map.svg"
  "${output_root}/report/matched_timepoint_companion_report.md"
)

for output in "${required_outputs[@]}"; do
  test -s "${output}"
done

"${rscript_cmd[@]}" -e '
  support <- read.csv(commandArgs(TRUE)[1])
  stopifnot(nrow(support) == 6L)
  stopifnot(all(support$matched_participants == 4L))
  stopifnot(!any(support$netmoss_network_feasible))
' "${output_root}/tables/timepoint_pair_support.csv"

echo "All matched-timepoint companion OpenPBS smoke tests passed."
