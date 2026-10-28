#!/usr/bin/env bash
# Run the three figure benchmarks at the same time, as independent R
# processes with their own logs. Each script uses a third of half the cores
# (CORES at its top), so together they stay within half the machine.
#
#   nohup bash scripts/run_figure_benchmarks.sh \
#     > "logs/figures_$(date +%F).log" 2>&1 &
#
# OpenBLAS 0.3.20 (Ubuntu 22.04) selects a faulty kernel on the Intel Xeon
# w9-3575X; the variables below must be set before R starts.
set -u
cd "$(dirname "$0")/.."
export OPENBLAS_CORETYPE=SkylakeX
export OPENBLAS_NUM_THREADS=1
mkdir -p logs
day=$(date +%F)

run() {
  local script=$1 log=$2
  local status=0
  echo "$(date '+%F %T') start ${script}"
  Rscript --no-save --no-restore "scripts/${script}" > "logs/${log}" 2>&1 ||
    status=$?
  echo "$(date '+%F %T') end ${script} (exit ${status})"
}

run fig05_single_cell_hierarchical.R "fig05_${day}_single_cell.log" &
run fig03_variance_covariance_driven.R "fig03_${day}_covariance_driven.log" &
run fig02_bivariate_toy.R "fig02_${day}_bivariate_toy.log" &
wait
