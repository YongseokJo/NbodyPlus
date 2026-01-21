#!/usr/bin/env bash
set -euo pipefail

# Common helpers for workflow scripts.

workflow_repo_root() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  cd "$script_dir/../.." && pwd
}

workflow_timestamp() {
  date +%Y%m%d_%H%M%S
}

workflow_die() {
  echo "ERROR: $*" >&2
  exit 1
}

workflow_load_config() {
  local repo_root
  repo_root="$(workflow_repo_root)"

  # shellcheck disable=SC1091
  source "$repo_root/workflow/config.sh"
  if [[ -f "$repo_root/workflow/config.local.sh" ]]; then
    # shellcheck disable=SC1091
    source "$repo_root/workflow/config.local.sh"
  fi

  # Explicit overrides from submit.sh (or schedulers) for subprocesses.
  # These are intentionally namespaced to avoid accidental overrides from the
  # ambient shell environment.
  if [[ -n "${WF_SCHEDULER_OVERRIDE:-}" ]]; then
    SCHEDULER="$WF_SCHEDULER_OVERRIDE"
  fi
  if [[ -n "${WF_TEST_DIR_OVERRIDE:-}" ]]; then
    TEST_DIR="$WF_TEST_DIR_OVERRIDE"
  fi
  if [[ -n "${WF_RUN_CONFIG_OVERRIDE:-}" ]]; then
    RUN_CONFIG="$WF_RUN_CONFIG_OVERRIDE"
  fi
  if [[ -n "${WF_USE_CUDA_OVERRIDE:-}" ]]; then
    USE_CUDA="$WF_USE_CUDA_OVERRIDE"
  fi
  if [[ -n "${WF_USE_SEVN_OVERRIDE:-}" ]]; then
    USE_SEVN="$WF_USE_SEVN_OVERRIDE"
  fi
  if [[ -n "${WF_ENABLE_PROFILING_OVERRIDE:-}" ]]; then
    ENABLE_PROFILING="$WF_ENABLE_PROFILING_OVERRIDE"
  fi
  if [[ -n "${WF_USE_MCLUSTER_OVERRIDE:-}" ]]; then
    USE_MCLUSTER="$WF_USE_MCLUSTER_OVERRIDE"
  fi
  if [[ -n "${WF_NTASKS_OVERRIDE:-}" ]]; then
    NTASKS="$WF_NTASKS_OVERRIDE"
  fi
  if [[ -n "${WF_GPUS_OVERRIDE:-}" ]]; then
    GPUS="$WF_GPUS_OVERRIDE"
  fi
  if [[ -n "${WF_PYTHON_OVERRIDE:-}" ]]; then
    PYTHON="$WF_PYTHON_OVERRIDE"
  fi
  if [[ -n "${WF_SUMMARY_FILE_OVERRIDE:-}" ]]; then
    SUMMARY_FILE="$WF_SUMMARY_FILE_OVERRIDE"
  fi
  if [[ -n "${WF_SUMMARY_STACK_FILE_OVERRIDE:-}" ]]; then
    SUMMARY_STACK_FILE="$WF_SUMMARY_STACK_FILE_OVERRIDE"
  fi
}

workflow_setup_env() {
  # Uses MPI_CANDIDATES/CUDA_CANDIDATES/HDF5_CANDIDATES from config.

  # MPI
  for mpi_path in "${MPI_CANDIDATES[@]:-}"; do
    if [[ -x "$mpi_path/bin/mpicxx" ]]; then
      export PATH="$mpi_path/bin:$PATH"
      export LD_LIBRARY_PATH="$mpi_path/lib:${LD_LIBRARY_PATH:-}"
      export CXX="$mpi_path/bin/mpicxx"
      break
    fi
  done
  # Fallback to PATH if no candidate matched
  if [[ -z "${CXX:-}" ]] && command -v mpicxx &>/dev/null; then
    export CXX="$(command -v mpicxx)"
  fi

  # HDF5
  if [[ -z "${HDF5_DIR:-}" ]]; then
    for hdf5_path in "${HDF5_CANDIDATES[@]:-}"; do
      if [[ -f "$hdf5_path/include/H5Cpp.h" ]]; then
        export HDF5_DIR="$hdf5_path"
        break
      fi
    done
  fi
  if [[ -n "${HDF5_DIR:-}" ]]; then
    export LD_LIBRARY_PATH="$HDF5_DIR/lib:${LD_LIBRARY_PATH:-}"
  fi

  # CUDA
  if [[ "${USE_CUDA:-0}" == "1" ]]; then
    if ! command -v nvcc &>/dev/null; then
      for cuda_path in "${CUDA_CANDIDATES[@]:-}"; do
        if [[ -x "$cuda_path/bin/nvcc" ]]; then
          export CUDA_HOME="$cuda_path"
          export PATH="$CUDA_HOME/bin:$PATH"
          export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
          break
        fi
      done
    fi

    if [[ -z "${CUDAHOSTCXX:-}" ]] && command -v g++ &>/dev/null; then
      export CUDAHOSTCXX="$(command -v g++)"
    fi
  fi

  # SEVN
  if [[ "${USE_SEVN:-0}" == "1" ]]; then
    if [[ -z "${SEVN_DIR:-}" ]]; then
      for sevn_path in "${SEVN_CANDIDATES[@]:-}"; do
        if [[ -d "$sevn_path/include/sevn" ]]; then
          export SEVN_DIR="$sevn_path"
          break
        fi
      done
    fi
    if [[ -n "${SEVN_DIR:-}" ]]; then
      export LD_LIBRARY_PATH="$SEVN_DIR/lib64/sevn:${LD_LIBRARY_PATH:-}"
    fi
  fi

  # Gfortran (for McLuster)
  if [[ "${USE_MCLUSTER:-0}" == "1" ]]; then
    if ! command -v gfortran &>/dev/null; then
      for gfortran_path in "${GFORTRAN_CANDIDATES[@]:-}"; do
        # Handle glob patterns in candidates
        for expanded_path in $gfortran_path; do
          if [[ -x "$expanded_path" ]]; then
            export PATH="$(dirname "$expanded_path"):$PATH"
            break 2
          fi
        done
      done
    fi
  fi
}

workflow_sanity() {
  [[ -n "${CXX:-}" ]] || workflow_die "CXX not set (check MPI_CANDIDATES in config.sh)"
  [[ -x "${CXX}" ]] || workflow_die "CXX=$CXX is not executable"
  [[ -n "${HDF5_DIR:-}" ]] || workflow_die "HDF5_DIR not set (set HDF5_DIR or HDF5_CANDIDATES)"

  if [[ "${USE_CUDA:-0}" == "1" ]]; then
    command -v nvcc &>/dev/null || workflow_die "nvcc not found (load CUDA or set CUDA_CANDIDATES)"
    [[ -n "${CUDA_HOME:-}" ]] || workflow_die "CUDA_HOME not set (set CUDA_HOME or CUDA_CANDIDATES)"

    if [[ -n "${CUDAHOSTCXX:-}" ]]; then
      if ! echo 'int main(){return 0;}' | "$CUDAHOSTCXX" -std=c++11 -x c++ -c -o /tmp/abyss_cuda_hostcxx_test.o - 2>/dev/null; then
        workflow_die "CUDAHOSTCXX ($CUDAHOSTCXX) does not accept -std=c++11"
      fi
      rm -f /tmp/abyss_cuda_hostcxx_test.o
    fi
  fi

  if [[ "${USE_SEVN:-0}" == "1" ]]; then
    [[ -n "${SEVN_DIR:-}" ]] || workflow_die "SEVN_DIR not set (set SEVN_DIR or SEVN_CANDIDATES)"
  fi

  if [[ "${USE_MCLUSTER:-0}" == "1" ]]; then
    if ! command -v gfortran &>/dev/null; then
      echo "Warning: gfortran not found, McLuster will not be built" >&2
      echo "  (Set USE_MCLUSTER=0 to suppress this warning)" >&2
    fi
  fi
}

workflow_run_dir() {
  local repo_root tag ts
  repo_root="$(workflow_repo_root)"
  ts="$(workflow_timestamp)"
  tag="${1:-run}"
  echo "$repo_root/workflow/runs/${tag}_${ts}"
}
