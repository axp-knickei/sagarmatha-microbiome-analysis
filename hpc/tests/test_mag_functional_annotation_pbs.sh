#!/usr/bin/env bash
# Static checks for the MAG functional OpenPBS wrapper and launcher.

set -euo pipefail
bash -n hpc/job_mag_functional_annotation.pbs
bash -n scripts/run_mag_functional_annotation.sh
grep -q '^#PBS -N ' hpc/job_mag_functional_annotation.pbs
grep -q '^#PBS -l ' hpc/job_mag_functional_annotation.pbs
grep -q 'PBS_O_WORKDIR' hpc/job_mag_functional_annotation.pbs
grep -q 'MAGFUNC_FROM_STAGE' hpc/job_mag_functional_annotation.pbs
grep -q 'MAGFUNC_TO_STAGE' hpc/job_mag_functional_annotation.pbs
grep -q 'MAGFUNC_TO_STAGE:-10' hpc/job_mag_functional_annotation.pbs
grep -q '08_decompose_shared_unique_carriers.R' scripts/run_mag_functional_annotation.sh
grep -q '09_analyze_shared_component_drivers.R' scripts/run_mag_functional_annotation.sh
grep -q '10_audit_candidate_annotations.R' scripts/run_mag_functional_annotation.sh
MAGFUNC_RSCRIPT=true bash scripts/run_mag_functional_annotation.sh \
  --config config/mag_functional_annotation.yml --from-stage 10 --to-stage 10
if grep -Eq '\b(sbatch|srun|squeue)\b' hpc/job_mag_functional_annotation.pbs; then
  echo "ERROR: Slurm command found in OpenPBS wrapper" >&2
  exit 1
fi
Rscript -e 'files <- list.files("analysis/mag_functional_annotation", pattern="[.]R$", recursive=TRUE, full.names=TRUE); for (f in files) parse(file=f)'
echo "MAG functional annotation OpenPBS static checks passed."
