#!/bin/bash

set -euo pipefail
output_root="${1:-/tmp/blackice_netmoss_altitude_stability_smoke}"
mkdir -p "${output_root}"

if [ -n "${NETMOSS_RSCRIPT:-}" ]; then
  read -r -a rscript_cmd <<<"${NETMOSS_RSCRIPT}"
else
  rscript_cmd=(Rscript)
fi

pbs_scripts=(
  hpc/job_netmoss_altitude_stability.pbs
  hpc/job_netmoss_altitude_stability_smoke_test.pbs
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
  invisible(parse(file="scripts/netmoss/altitude_stability_helpers.R"))
  invisible(parse(file="scripts/netmoss/run_altitude_stability.R"))
  invisible(parse(file="scripts/netmoss/summarize_altitude_stability.R"))
  cat("Altitude-stability R parse checks passed\n")
'
"${rscript_cmd[@]}" scripts/netmoss/tests/testthat-altitude-stability.R
"${rscript_cmd[@]}" scripts/netmoss/run_altitude_stability.R \
  --config config/netmoss_altitude_stability.yml \
  --output-dir "${output_root}" \
  --test-mode \
  --bootstrap 3 \
  --threads 2 \
  --acknowledge-exploratory \
  --force
"${rscript_cmd[@]}" scripts/netmoss/summarize_altitude_stability.R \
  --config config/netmoss_altitude_stability.yml \
  --output-dir "${output_root}" \
  --force

repeat_root="${output_root}_repeat"
"${rscript_cmd[@]}" scripts/netmoss/run_altitude_stability.R \
  --config config/netmoss_altitude_stability.yml \
  --output-dir "${repeat_root}" \
  --test-mode \
  --bootstrap 3 \
  --threads 2 \
  --acknowledge-exploratory \
  --force
for table in scenario_summary.csv replicate_summary.csv edge_stability.csv taxon_rank_stability.csv cross_rank_concordance.csv; do
  cmp "${output_root}/tables/${table}" "${repeat_root}/tables/${table}"
done

required_outputs=(
  "${output_root}/tables/scenario_summary.csv"
  "${output_root}/tables/cohort_membership.csv"
  "${output_root}/tables/edge_stability.csv"
  "${output_root}/tables/module_stability.csv"
  "${output_root}/tables/taxon_rank_stability.csv"
  "${output_root}/tables/participant_influence.csv"
  "${output_root}/tables/cross_rank_concordance.csv"
  "${output_root}/reproducibility/analysis_manifest.yml"
  "${output_root}/figures/scenario_network_stability.pdf"
  "${output_root}/figures/primary_taxon_rank_stability.svg"
  "${output_root}/report/netmoss_altitude_stability_report.md"
)
for output in "${required_outputs[@]}"; do
  test -s "${output}"
done

echo "All NetMoss altitude stability smoke tests passed."
