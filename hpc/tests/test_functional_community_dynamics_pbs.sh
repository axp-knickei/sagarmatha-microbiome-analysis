#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
bash -n hpc/job_functional_community_dynamics.pbs
bash -n scripts/run_functional_community_dynamics.sh
grep -q 'FUNCDYN_FROM_STAGE' hpc/job_functional_community_dynamics.pbs
grep -q 'FUNCDYN_TO_STAGE' hpc/job_functional_community_dynamics.pbs
grep -q 'FROM_STAGE=0' scripts/run_functional_community_dynamics.sh
grep -q 'TO_STAGE=4' scripts/run_functional_community_dynamics.sh
grep -q 'must be 0-4' scripts/run_functional_community_dynamics.sh
echo "Functional community dynamics PBS wrapper tests passed"
