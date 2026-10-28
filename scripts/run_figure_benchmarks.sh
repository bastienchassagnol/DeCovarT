#!/usr/bin/env bash
# Run figure benchmarks at the same time, as independent R processes with
# their own logs. Each script uses a third of half the cores (CORES at its
# top), so three together stay within half the machine.
#
#   nohup bash scripts/run_figure_benchmarks.sh \
#     > "logs/figures_$(date +%F).log" 2>&1 &
#   # a subset: bash scripts/run_figure_benchmarks.sh fig02 fig05
#
# OpenBLAS 0.3.20 (Ubuntu 22.04) selects a faulty kernel on the Intel Xeon
# w9-3575X; the variables below must be set before R starts.
set -u
cd "$(dirname "$0")/.."
export OPENBLAS_CORETYPE=SkylakeX
export OPENBLAS_NUM_THREADS=1
mkdir -p logs
day=$(date +%F)

declare -A scripts=(
  [fig02]=fig02_bivariate_toy
  [fig03]=fig03_variance_covariance_driven
  [fig05]=fig05_single_cell_hierarchical
)
declare -A labels=(
  [fig02]=bivariate_toy
  [fig03]=covariance_driven
  [fig05]=single_cell
)

run() {
  local figure=$1
  local script="${scripts[$figure]}.R"
  local status=0
  echo "$(date '+%F %T') start ${script}"
  Rscript --no-save --no-restore "scripts/${script}" \
    > "logs/${figure}_${day}_${labels[$figure]}.log" 2>&1 || status=$?
  echo "$(date '+%F %T') end ${script} (exit ${status})"
}

figures=("$@")
if [ ${#figures[@]} -eq 0 ]; then
  figures=(fig05 fig03 fig02)
fi
for figure in "${figures[@]}"; do
  if [ -z "${scripts[$figure]+x}" ]; then
    echo "Unknown figure ${figure}; expected fig02, fig03 or fig05." >&2
    exit 2
  fi
done
for figure in "${figures[@]}"; do
  run "$figure" &
done
wait
