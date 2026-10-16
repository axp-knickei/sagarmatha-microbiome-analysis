#!/bin/bash

# Static and executable smoke tests for every NetMoss OpenPBS wrapper.

set -euo pipefail

output_root="${1:-/tmp/blackice_netmoss_hpc_smoke}"
mkdir -p "$output_root"

if [ -n "${NETMOSS_RSCRIPT:-}" ]; then
  read -r -a rscript_cmd <<<"${NETMOSS_RSCRIPT}"
else
  rscript_cmd=(Rscript)
fi

pbs_scripts=(
  hpc/update_netmoss_env.pbs
  hpc/job_netmoss_readiness.pbs
  hpc/job_netmoss_analysis.pbs
  hpc/job_netmoss_altitude_analysis.pbs
  hpc/job_netmoss_altitude_stability.pbs
  hpc/job_netmoss_altitude_stability_smoke_test.pbs
  hpc/job_netmoss_publication.pbs
  hpc/job_netmoss_smoke_test.pbs
)

for script in "${pbs_scripts[@]}"; do
  test -s "$script"
  bash -n "$script"
  grep -q '^#PBS -N ' "$script"
  grep -q '^#PBS -l ' "$script"
  grep -q 'PBS_O_WORKDIR' "$script"
  if grep -Eq '\b(sbatch|srun|squeue)\b' "$script"; then
    echo "ERROR: Slurm command found in OpenPBS script: $script" >&2
    exit 1
  fi
done

bash -n hpc/netmoss_common.sh
"${rscript_cmd[@]}" -e 'invisible(parse(file="scripts/check_netmoss_readiness.R")); invisible(parse(file="scripts/netmoss/run_netmoss_pipeline.R")); invisible(parse(file="scripts/netmoss/create_publication_outputs.R")); invisible(parse(file="scripts/netmoss/validate_environment.R")); cat("R parse checks passed\n")'
"${rscript_cmd[@]}" scripts/netmoss/validate_environment.R

"${rscript_cmd[@]}" scripts/check_netmoss_readiness.R \
  --abundance data/test_beta_diversity/gut_test.csv \
  --metadata data/test_beta_diversity/sample_info_test.csv \
  --taxonomy data/test_beta_diversity/GDTB_test.csv \
  --output-dir "$output_root/readiness" \
  --force

"${rscript_cmd[@]}" scripts/netmoss/run_netmoss_pipeline.R \
  --config config/netmoss_analysis.yml \
  --output-dir "$output_root/prepared" \
  --test-mode \
  --prepare-only \
  --force

"${rscript_cmd[@]}" scripts/netmoss/create_publication_outputs.R \
  --config config/netmoss_analysis.yml \
  --output-dir "$output_root/publication" \
  --test-mode \
  --force

required_outputs=(
  "$output_root/readiness/netmoss_readiness_summary.txt"
  "$output_root/prepared/reproducibility/netmoss_manifest.yml"
  "$output_root/prepared/networks/group_a_network.csv"
  "$output_root/prepared/networks/group_b_network.csv"
  "$output_root/publication/figures/Figure_1_NetMoss_network_rewiring.pdf"
  "$output_root/publication/figures/Figure_1_NetMoss_network_rewiring.svg"
  "$output_root/publication/figures/Figure_1_NetMoss_network_rewiring.tiff"
  "$output_root/publication/tables/Supplementary_Data_1_NetMoss_taxon_results.csv"
  "$output_root/publication/tables/Supplementary_Data_2_NetMoss_network_edges.csv"
  "$output_root/publication/text/NetMoss_methods_draft.txt"
)

for output in "${required_outputs[@]}"; do
  test -s "$output"
done

echo "All NetMoss OpenPBS smoke tests passed."
