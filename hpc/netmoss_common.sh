#!/bin/bash

# Shared OpenPBS/Conda initialization for NetMoss jobs. Source; do not execute.

setup_netmoss_conda() {
  local conda_module="${CONDA_MODULE_NAME:-intelpython/2022.2.0}"
  local conda_env="${CONDA_ENV_NAME:-glmm-hpc}"

  if command -v module >/dev/null 2>&1; then
    module purge || true
    module load "${conda_module}"
  fi

  if ! command -v conda >/dev/null 2>&1; then
    local profile
    for profile in \
      "$HOME/miniconda3/etc/profile.d/conda.sh" \
      "$HOME/anaconda3/etc/profile.d/conda.sh" \
      "$HOME/mambaforge/etc/profile.d/conda.sh" \
      "$HOME/miniforge3/etc/profile.d/conda.sh"
    do
      if [ -f "$profile" ]; then
        # shellcheck disable=SC1090
        source "$profile"
        break
      fi
    done
  fi

  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found; set CONDA_MODULE_NAME or initialize Conda." >&2
    return 1
  fi

  export NETMOSS_CONDA_ENV="$conda_env"
  export NETMOSS_RSCRIPT="conda run --no-capture-output -n ${conda_env} Rscript"
  echo "Conda environment: ${conda_env}"
  conda run -n "${conda_env}" Rscript --version
}

print_netmoss_job_context() {
  echo "Host: $(hostname)"
  echo "Working directory: $(pwd)"
  echo "PBS job ID: ${PBS_JOBID:-local}"
  echo "Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
}
