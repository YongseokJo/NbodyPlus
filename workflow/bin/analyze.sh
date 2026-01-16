#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=${1:?run_dir}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

workflow_load_config
workflow_setup_env

RUN_LOG="$RUN_DIR/run.log"
BUILD_LOG="$RUN_DIR/build.log"

# Prefer run-specific settings recorded in meta.txt when present.
if [[ -f "$RUN_DIR/meta.txt" ]]; then
  meta_python="$(grep -E '^python=' "$RUN_DIR/meta.txt" | tail -n 1 | sed 's/^python=//' || true)"
  if [[ -n "${meta_python:-}" ]]; then
    PYTHON="$meta_python"
  fi
  meta_summary="$(grep -E '^summary_file=' "$RUN_DIR/meta.txt" | tail -n 1 | sed 's/^summary_file=//' || true)"
  if [[ -n "${meta_summary:-}" ]]; then
    SUMMARY_FILE="$meta_summary"
  fi
  meta_stack="$(grep -E '^summary_stack_file=' "$RUN_DIR/meta.txt" | tail -n 1 | sed 's/^summary_stack_file=//' || true)"
  if [[ -n "${meta_stack:-}" ]]; then
    SUMMARY_STACK_FILE="$meta_stack"
  fi
fi

SUMMARY_NAME="${SUMMARY_FILE:-summary.txt}"
if [[ "$SUMMARY_NAME" = /* ]]; then
  SUMMARY="$SUMMARY_NAME"
else
  SUMMARY="$RUN_DIR/$SUMMARY_NAME"
fi

{
  echo "== ABYSS workflow summary =="
  echo "date=$(date)"
  echo "run_dir=$RUN_DIR"
  if [[ -f "$RUN_DIR/meta.txt" ]]; then
    echo "-- run metadata --"
    cat "$RUN_DIR/meta.txt"
    echo ""
  fi
  if [[ -n "${USE_CUDA:-}" ]]; then
    arch_mode="CPU"
    if [[ "${USE_CUDA:-0}" == "1" ]]; then
      arch_mode="GPU"
    fi
    echo "-- architecture --"
    echo "mode=$arch_mode"
    echo "nodes=${NODES:-}"
    echo "ntasks=${NTASKS:-}"
    echo "gpus=${GPUS:-}"
    echo ""
  fi
  echo ""

  if [[ -f "$BUILD_LOG" ]]; then
    echo "-- build --"
    tail -n 20 "$BUILD_LOG" || true
    echo ""
  fi

  if [[ -f "$RUN_LOG" ]]; then
    echo "-- run --"
    if grep -q "Simulation Done" "$RUN_LOG"; then
      echo "status=ok (Simulation Done found)"
    else
      echo "status=unknown (Simulation Done not found)"
    fi

    echo ""
    echo "-- errors/warnings (tail) --"
    grep -nE "ERROR|Error|FATAL|Fatal|assert|Assertion|Segmentation|MPI_ABORT|cudaError" "$RUN_LOG" | tail -n 50 || true
    echo ""
  else
    echo "No run log found at $RUN_LOG"
  fi

  # Show produced outputs (if any)
  if [[ -d "$RUN_DIR/work" ]]; then
    echo "-- work dir artifacts --"
    (cd "$RUN_DIR/work" && find . -maxdepth 2 -type f \( -name "*.txt" -o -name "*.csv" -o -name "*.log" \) | head -n 200) || true
  fi
} > "$SUMMARY"

echo "Wrote: $SUMMARY"

# Optional: run analysis tools under tools/ with the configured Python.
# Writes per-tool logs next to the run summary.
REPO_ROOT="$(workflow_repo_root)"

PYTHON_BIN="${PYTHON:-python3}"

TOOLS_LOG="$RUN_DIR/tools_analysis.log"
{
  echo "== ABYSS tools analysis =="
  echo "date=$(date)"
  echo "python=$PYTHON_BIN"
  echo ""
} > "$TOOLS_LOG"

if declare -p ANALYZE_TOOLS &>/dev/null && (( ${#ANALYZE_TOOLS[@]} > 0 )); then
  if ! command -v "$PYTHON_BIN" &>/dev/null; then
    echo "ERROR: python not found: $PYTHON_BIN" >> "$TOOLS_LOG"
    exit 0
  fi

  for tool in "${ANALYZE_TOOLS[@]}"; do
    tool_path="$REPO_ROOT/$tool"
    if [[ ! -f "$tool_path" ]]; then
      echo "SKIP: missing tool: $tool" >> "$TOOLS_LOG"
      continue
    fi

    out_base="$RUN_DIR/$(basename "$tool" .py)"
    echo "-- running: $tool --" >> "$TOOLS_LOG"

    # Choose tool arguments based on common conventions.
    work_dir="$RUN_DIR/work"
    out_dir="$work_dir/output"
    tool_args=()
    case "$(basename "$tool")" in
      analyze_profiling.py)
        if [[ -d "$out_dir" ]]; then
          tool_args+=("$out_dir")
        else
          tool_args+=("$work_dir")
        fi
        ;;
      analyze_energy.py)
        h5_file=""
        if [[ -d "$out_dir" ]]; then
          h5_file="$(find "$out_dir" -maxdepth 2 -type f -name '*.h5' | head -n 1)"
        fi
        if [[ -z "$h5_file" ]]; then
          h5_file="$(find "$work_dir" -maxdepth 3 -type f -name '*.h5' | head -n 1)"
        fi
        if [[ -z "$h5_file" ]]; then
          echo "SKIP: no .h5 found under $work_dir (needed by analyze_energy.py)" >> "$TOOLS_LOG"
          continue
        fi
        tool_args+=("$h5_file")
        ;;
      *)
        tool_args+=("$work_dir")
        ;;
    esac

    # Conventions:
    # - pass the run work directory if the tool accepts it
    # - always capture stdout/stderr to a per-tool log
    set +e
    "$PYTHON_BIN" "$tool_path" "${tool_args[@]}" >"${out_base}.out" 2>"${out_base}.err"
    rc=$?
    set -e

    echo "rc=$rc" >> "$TOOLS_LOG"
    if [[ $rc -ne 0 ]]; then
      echo "  stderr: ${out_base}.err" >> "$TOOLS_LOG"
    else
      echo "  stdout: ${out_base}.out" >> "$TOOLS_LOG"
    fi
  done
else
  echo "No ANALYZE_TOOLS configured; summary only." >> "$TOOLS_LOG"
fi

echo "Wrote: $TOOLS_LOG"

# Optional: append concise summary from SUMMARY_TOOL
if [[ -n "${SUMMARY_TOOL:-}" ]]; then
  summary_tool_path="$REPO_ROOT/$SUMMARY_TOOL"
  if [[ -f "$summary_tool_path" ]]; then
    if command -v "$PYTHON_BIN" &>/dev/null; then
      echo "" >> "$SUMMARY"
      echo "-- concise performance/energy summary --" >> "$SUMMARY"
      set +e
      "$PYTHON_BIN" "$summary_tool_path" "$RUN_DIR" >> "$SUMMARY" 2>> "$TOOLS_LOG"
      rc=$?
      set -e
      if [[ $rc -ne 0 ]]; then
        echo "WARN: SUMMARY_TOOL failed (rc=$rc): $SUMMARY_TOOL" >> "$TOOLS_LOG"
      fi
    else
      echo "WARN: python not found for SUMMARY_TOOL: $PYTHON_BIN" >> "$TOOLS_LOG"
    fi
  else
    echo "WARN: missing SUMMARY_TOOL: $SUMMARY_TOOL" >> "$TOOLS_LOG"
  fi
fi

# Append a TSV row for stacked summaries (project-level by default).
if [[ -n "${SUMMARY_STACK_FILE:-}" && -n "${SUMMARY_TOOL:-}" ]]; then
  stack_path="$SUMMARY_STACK_FILE"
  if [[ "$stack_path" != /* ]]; then
    stack_path="$REPO_ROOT/$stack_path"
  fi

  if [[ -f "$summary_tool_path" && -x "$(command -v "$PYTHON_BIN")" ]]; then
    header="$($PYTHON_BIN "$summary_tool_path" --tsv-header-pretty 2>>"$TOOLS_LOG")"
    if [[ -z "$header" ]]; then
      echo "WARN: empty TSV header from SUMMARY_TOOL" >> "$TOOLS_LOG"
    else
      if [[ -f "$stack_path" ]]; then
        first_line="$(head -n 1 "$stack_path" | tr -d '\r')"
        if [[ "$first_line" != "$header" ]]; then
          ts="$(date +%Y%m%d_%H%M%S)"
          bak="${stack_path}.bak_${ts}"
          mv "$stack_path" "$bak"
          echo "Rotated stack file to: $bak" >> "$TOOLS_LOG"
        fi
      fi

      if [[ ! -f "$stack_path" ]]; then
        echo "$header" > "$stack_path"
      fi

      "$PYTHON_BIN" "$summary_tool_path" --tsv-row-pretty "$RUN_DIR" >> "$stack_path" 2>> "$TOOLS_LOG" || true
      echo "Wrote: $stack_path" >> "$TOOLS_LOG"
    fi
  fi
fi
