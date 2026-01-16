#!/usr/bin/env python3
"""
Summarize ABYSS performance and energy conservation for a run directory.

Usage:
    python summarize_run.py <run_dir>
    python summarize_run.py --tsv-header
    python summarize_run.py --tsv-row <run_dir>

The script looks for:
- profiling.csv under <run_dir>/work/output
- output.h5 under <run_dir>/work/output

Outputs a concise text summary to stdout, or a TSV header/row for stacking.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
from pathlib import Path

np = None
h5py = None

# Units from src/def.h (aligned with tools/analyze_energy.py)
POSITION_UNIT_PC = 4.0
TIME_UNIT_YR = 1e10
VELOCITY_UNIT_PC_YR = 4e-10
MASS_UNIT_MSUN = 0.0001424198

PC_IN_KM = 3.08567758149137e13
SEC_PER_YR = 3.1536e7
KM_S_TO_PC_YR = SEC_PER_YR / PC_IN_KM


def find_run_artifacts(run_dir: Path):
    work_dir = run_dir / "work"
    output_dir = work_dir / "output"
    profiling_csv = None
    output_h5 = None

    if output_dir.is_dir():
        prof = list(output_dir.glob("profiling.csv"))
        if prof:
            profiling_csv = prof[0]
        h5s = list(output_dir.glob("*.h5"))
        if h5s:
            output_h5 = h5s[0]

    if output_h5 is None and work_dir.is_dir():
        h5s = list(work_dir.glob("**/*.h5"))
        if h5s:
            output_h5 = h5s[0]

    return profiling_csv, output_h5


def read_meta(run_dir: Path):
    meta = {}
    meta_path = run_dir / "meta.txt"
    if not meta_path.exists():
        return meta
    for line in meta_path.read_text().splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        meta[k.strip()] = v.strip()
    return meta


def extract_tag(run_dir: Path):
    name = run_dir.name
    m = re.search(r"_(\d{8}_\d{6})$", name)
    if m:
        return name[: m.start()]
    return name


def cpu_model():
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8") as f:
            for line in f:
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return ""


def gpu_model():
    try:
        out = subprocess.check_output(["nvidia-smi", "-L"], text=True, stderr=subprocess.DEVNULL)
        if out:
            # Example: GPU 0: NVIDIA A100-SXM4-40GB (UUID: ...)
            line = out.strip().splitlines()[0]
            if ":" in line:
                return line.split(":", 1)[1].split("(", 1)[0].strip()
            return line.strip()
    except Exception:
        pass
    return ""


def gpu_arch_from_meta(meta: dict):
    """Return a trustworthy gpu_arch from meta.txt, or empty if it looks bogus."""
    gpu = (meta.get("gpu_arch", "") or "").strip()
    if not gpu:
        return ""
    bad_substrings = [
        "NVIDIA-SMI has failed",
        "VGA compatible controller:",
        "ASPEED",
    ]
    for bad in bad_substrings:
        if bad in gpu:
            return ""
    return gpu


def summarize_profiling_csv(csv_path: Path):
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        return None

    # Timer columns end with _ns
    timer_cols = [c for c in rows[0].keys() if c.endswith("_ns")]
    sums = {c: 0.0 for c in timer_cols}
    sim_time = []
    for r in rows:
        for c in timer_cols:
            val = r.get(c, "")
            if val != "":
                sums[c] += float(val)
        if "sim_time_myr" in r and r["sim_time_myr"] != "":
            sim_time.append(float(r["sim_time_myr"]))

    total_ns = sums.get("WholeRoutine_ns", sum(sums.values()))
    total_s = total_ns * 1e-9 if total_ns else 0.0

    # Build sorted list
    items = []
    for c in timer_cols:
        s = sums[c] * 1e-9
        pct = 100.0 * s / total_s if total_s > 0 else 0.0
        items.append((c.replace("_ns", ""), s, pct))
    items.sort(key=lambda x: x[1], reverse=True)

    sim_time_min = min(sim_time) if sim_time else None
    sim_time_max = max(sim_time) if sim_time else None

    return {
        "total_s": total_s,
        "items": items,
        "sim_time_min": sim_time_min,
        "sim_time_max": sim_time_max,
        "num_rows": len(rows),
    }


def compute_energy_code_units(m_msun, pos_pc, vel_kms):
    global np
    if np is None:
        try:
            import numpy as _np
            np = _np
        except Exception:
            raise RuntimeError("numpy is required for energy computations")

    m_code = m_msun / MASS_UNIT_MSUN
    pos_code = pos_pc / POSITION_UNIT_PC
    vel_pcyr = vel_kms * KM_S_TO_PC_YR
    vel_code = vel_pcyr / VELOCITY_UNIT_PC_YR

    kinetic = 0.5 * np.sum(m_code * np.sum(vel_code ** 2, axis=1))

    n = pos_code.shape[0]
    if n > 5000:
        potential = 0.0
        chunk_size = 1000
        for i in range(0, n, chunk_size):
            i_end = min(i + chunk_size, n)
            for j in range(i, n, chunk_size):
                j_end = min(j + chunk_size, n)
                diff = pos_code[i:i_end, None, :] - pos_code[None, j:j_end, :]
                r = np.linalg.norm(diff, axis=2)
                if i == j:
                    iu = np.triu_indices(i_end - i, k=1)
                    rij = r[iu]
                    inv_r = np.where(rij > 0.0, 1.0 / rij, 0.0)
                    potential -= np.sum(m_code[i:i_end][iu[0]] * m_code[i:i_end][iu[1]] * inv_r)
                else:
                    inv_r = np.where(r > 0.0, 1.0 / r, 0.0)
                    potential -= np.sum(
                        m_code[i:i_end][:, None] * m_code[j:j_end][None, :] * inv_r
                    )
    else:
        diff = pos_code[:, None, :] - pos_code[None, :, :]
        r = np.linalg.norm(diff, axis=2)
        with np.errstate(divide='ignore', invalid='ignore'):
            inv_r = np.where(r > 0.0, 1.0 / r, 0.0)
        potential = -0.5 * np.sum(m_code[:, None] * m_code[None, :] * inv_r)

    total = kinetic + potential
    return kinetic, potential, total


def summarize_energy(h5_path: Path):
    global h5py, np
    if h5py is None:
        try:
            import h5py as _h5py
            h5py = _h5py
        except Exception:
            return {"error": "h5py not available"}

    if np is None:
        try:
            import numpy as _np
            np = _np
        except Exception:
            return {"error": "numpy not available"}

    with h5py.File(h5_path, "r") as f:
        steps = [k for k in f.keys() if k.startswith("Step_")]
        if not steps:
            return {"error": "no Step_* groups found"}
        steps.sort(key=lambda x: int(x.split("_")[1]))

        times = []
        totals = []
        for k in steps:
            step = f[k]
            time_myr = float(step.attrs.get("Time_Myr", 0.0))
            masses = step["Mass_Msun"][:]
            pos = np.column_stack([step["X_pc"][:], step["Y_pc"][:], step["Z_pc"][:]])
            vel = np.column_stack([step["Vx_km_s"][:], step["Vy_km_s"][:], step["Vz_km_s"][:]])
            _, _, total = compute_energy_code_units(masses, pos, vel)
            times.append(time_myr)
            totals.append(total)

    totals = np.asarray(totals, dtype=float)
    times = np.asarray(times, dtype=float)

    e0 = totals[0] if totals.size else 0.0
    if e0 != 0.0:
        dE = (totals - e0) / abs(e0)
    else:
        dE = np.zeros_like(totals)

    return {
        "n_steps": len(totals),
        "time_min": float(times.min()) if times.size else None,
        "time_max": float(times.max()) if times.size else None,
        "energy_mean": float(totals.mean()) if totals.size else None,
        "energy_std": float(totals.std()) if totals.size else None,
        "energy_min": float(totals.min()) if totals.size else None,
        "energy_max": float(totals.max()) if totals.size else None,
        "dE_mean": float(dE.mean()) if dE.size else None,
        "dE_std": float(dE.std()) if dE.size else None,
        "dE_min": float(dE.min()) if dE.size else None,
        "dE_max": float(dE.max()) if dE.size else None,
    }


def tsv_header():
    return "\t".join(tsv_fields())


def tsv_fields():
    return [
        "tag",
        "test_name",
        "simulation_duration_myr",
        "total_wall_s",
        "dE_over_E0_mean",
        "dE_over_E0_std",
        "scheduler",
        "cpu_arch",
        "gpu_arch",
        "nodes",
        "ntasks",
        "gpus",
        "timestamp",
        "run_dir",
        "git_commit",
        "git_commit_long",
        "git_branch",
        "git_tag",
    ]


def find_repo_root(start: Path) -> Path | None:
    path = start.resolve()
    if path.is_file():
        path = path.parent
    for candidate in [path, *path.parents]:
        if (candidate / ".git").exists():
            return candidate
    return None


def git_info(run_dir: Path, meta: dict | None = None):
    meta = meta or {}
    repo_root = find_repo_root(run_dir)
    if repo_root is None:
        return (
            meta.get("git_commit", ""),
            meta.get("git_commit_long", ""),
            meta.get("git_branch", ""),
            meta.get("git_tag", ""),
        )

    base_cmd = ["git", "-C", str(repo_root)]
    try:
        commit_short = subprocess.check_output(
            base_cmd + ["rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        commit_short = ""

    try:
        commit_long = subprocess.check_output(
            base_cmd + ["rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        commit_long = ""

    try:
        branch = subprocess.check_output(
            base_cmd + ["rev-parse", "--abbrev-ref", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if branch == "HEAD":
            branch = ""
    except Exception:
        branch = ""

    try:
        tag = subprocess.check_output(
            base_cmd + ["describe", "--tags", "--exact-match"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        tag = ""

    if not tag and commit_short:
        tag = commit_short

    if not commit_short:
        commit_short = meta.get("git_commit", "")
    if not commit_long:
        commit_long = meta.get("git_commit_long", "")
    if not branch:
        branch = meta.get("git_branch", "")
    if not tag:
        tag = meta.get("git_tag", "") or commit_short

    return commit_short, commit_long, branch, tag


def tsv_header_pretty(min_width: int = 12, sep: str = "  "):
    fields = tsv_fields()
    widths = [max(len(f), min_width) for f in fields]
    cols = [f.ljust(w) for f, w in zip(fields, widths)]
    return sep.join(cols)


def build_row(run_dir: Path):
    meta = read_meta(run_dir)
    profiling_csv, output_h5 = find_run_artifacts(run_dir)
    perf = summarize_profiling_csv(profiling_csv) if profiling_csv else None
    energy = summarize_energy(output_h5) if output_h5 else None

    tag = extract_tag(run_dir)
    sim_time = None
    if energy and energy.get("time_max") is not None:
        sim_time = energy["time_max"]
    elif perf and perf.get("sim_time_max") is not None:
        sim_time = perf["sim_time_max"]

    cpu = meta.get("cpu_arch", "") or cpu_model()
    gpu = gpu_arch_from_meta(meta) or gpu_model()
    test_name = meta.get("test_dir", "")
    commit_short, commit_long, branch, tag_git = git_info(run_dir, meta)

    def sanitize(v):
        if v is None:
            return ""
        s = str(v)
        # collapse newlines and repeated whitespace to single spaces
        s = s.replace("\n", " ").replace("\r", " ")
        s = " ".join(s.split())
        return s

    row = [
        tag,
        test_name,
        format_num(sim_time),
        format_num(perf["total_s"] if perf else None),
        format_num(energy["dE_mean"] if energy and "error" not in energy else None),
        format_num(energy["dE_std"] if energy and "error" not in energy else None),
        meta.get("scheduler", ""),
        cpu,
        gpu,
        meta.get("nodes", ""),
        meta.get("ntasks", ""),
        meta.get("gpus", ""),
        meta.get("timestamp", ""),
        str(run_dir),
        commit_short,
        commit_long,
        branch,
        tag_git,
    ]
    return [sanitize(x) for x in row]


def pretty_table(rows: list[list[str]], min_width: int = 0, sep: str = "  "):
    fields = tsv_fields()
    widths = [max(len(f), min_width) for f in fields]
    for row in rows:
        for i, value in enumerate(row):
            widths[i] = max(widths[i], len(str(value)))

    numeric_fields = {
        "simulation_duration_myr",
        "total_wall_s",
        "dE_over_E0_mean",
        "dE_over_E0_std",
        "nodes",
        "ntasks",
        "gpus",
    }

    header = sep.join([f.ljust(w) for f, w in zip(fields, widths)])
    lines = [header]
    for row in rows:
        cols = []
        for field, value, width in zip(fields, row, widths):
            text = str(value)
            if field in numeric_fields and text:
                cols.append(text.rjust(width))
            else:
                cols.append(text.ljust(width))
        lines.append(sep.join(cols))
    return "\n".join(lines)


def tsv_row_pretty(run_dir: Path, min_width: int = 12, sep: str = "  "):
    # Build the same row as --tsv-row but then align using header-based widths.
    tag = extract_tag(run_dir)
    profiling_csv, output_h5 = find_run_artifacts(run_dir)
    perf = summarize_profiling_csv(profiling_csv) if profiling_csv else None
    energy = summarize_energy(output_h5) if output_h5 else None

    sim_duration = None
    if energy and energy.get("time_max") is not None:
        sim_duration = energy["time_max"]
    elif perf and perf.get("sim_time_max") is not None:
        sim_duration = perf["sim_time_max"]
    meta = read_meta(run_dir)
    cpu = meta.get("cpu_arch", "") or cpu_model()
    gpu = gpu_arch_from_meta(meta) or gpu_model()
    test_name = meta.get("test_dir", "")
    row = [
        tag,
        test_name,
        format_num(sim_duration),
        format_num(perf["total_s"] if perf else None),
        format_num(energy["dE_mean"] if energy and "error" not in energy else None),
        format_num(energy["dE_std"] if energy and "error" not in energy else None),
        meta.get("scheduler", ""),
        cpu,
        gpu,
        meta.get("nodes", ""),
        meta.get("ntasks", ""),
        meta.get("gpus", ""),
        meta.get("timestamp", ""),
        str(run_dir),
        "",
        "",
        "",
        "",
    ]

    fields = tsv_fields()
    widths = [max(len(f), min_width) for f in fields]

    commit_short, commit_long, branch, tag = git_info(run_dir, meta)
    row[-4] = commit_short
    row[-3] = commit_long
    row[-2] = branch
    row[-1] = tag

    numeric_fields = {
        "simulation_duration_myr",
        "total_wall_s",
        "dE_over_E0_mean",
        "dE_over_E0_std",
        "nodes",
        "ntasks",
        "gpus",
    }
    cols = []
    for field, value, width in zip(fields, row, widths):
        text = str(value)
        if field in numeric_fields and text:
            cols.append(text.rjust(width))
        else:
            cols.append(text.ljust(width))
    return sep.join(cols)


def format_num(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", nargs="?")
    parser.add_argument("--tsv-header", action="store_true")
    parser.add_argument("--tsv-row", action="store_true")
    parser.add_argument("--tsv-header-pretty", action="store_true")
    parser.add_argument("--tsv-row-pretty", action="store_true")
    parser.add_argument("--stack-pretty", metavar="RUNS_ROOT")
    parser.add_argument("--stack-output", metavar="OUTPUT")
    args = parser.parse_args()

    if args.tsv_header:
        print(tsv_header())
        return

    if args.tsv_header_pretty:
        print(tsv_header_pretty())
        return

    if args.stack_pretty:
        runs_root = Path(args.stack_pretty).resolve()
        if not runs_root.exists():
            print(f"ERROR: runs root not found: {runs_root}")
            sys.exit(2)

        run_dirs = sorted([p for p in runs_root.iterdir() if p.is_dir()])
        rows = [build_row(p) for p in run_dirs]
        output = pretty_table(rows)
        if args.stack_output:
            Path(args.stack_output).write_text(output + "\n")
        else:
            print(output)
        return

    if not args.run_dir:
        print("Usage: summarize_run.py <run_dir> | --tsv-header | --tsv-row <run_dir>")
        sys.exit(2)

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.exists():
        print(f"ERROR: run_dir not found: {run_dir}")
        sys.exit(2)

    meta = read_meta(run_dir)
    profiling_csv, output_h5 = find_run_artifacts(run_dir)
    perf = summarize_profiling_csv(profiling_csv) if profiling_csv else None
    energy = summarize_energy(output_h5) if output_h5 else None

    if args.tsv_row:
        tag = extract_tag(run_dir)
        sim_duration = None
        if energy and energy.get("time_max") is not None:
            sim_duration = energy["time_max"]
        elif perf and perf.get("sim_time_max") is not None:
            sim_duration = perf["sim_time_max"]
        cpu = meta.get("cpu_arch", "") or cpu_model()
        gpu = gpu_arch_from_meta(meta) or gpu_model()
        test_name = meta.get("test_dir", "")
        row = [
            tag,
            test_name,
            format_num(sim_duration),
            format_num(perf["total_s"] if perf else None),
            format_num(energy["dE_mean"] if energy and "error" not in energy else None),
            format_num(energy["dE_std"] if energy and "error" not in energy else None),
            meta.get("scheduler", ""),
            cpu,
            gpu,
            meta.get("nodes", ""),
            meta.get("ntasks", ""),
            meta.get("gpus", ""),
            meta.get("timestamp", ""),
            str(run_dir),
        ]
        commit_short, commit_long, branch, tag = git_info(run_dir, meta)
        # sanitize and append git info
        def _sanitize(v):
            if v is None:
                return ""
            s = str(v)
            s = s.replace("\n", " ").replace("\r", " ")
            s = " ".join(s.split())
            return s

        row.extend([_sanitize(commit_short), _sanitize(commit_long), _sanitize(branch), _sanitize(tag)])
        print("\t".join(row))
        return

    if args.tsv_row_pretty:
        print(tsv_row_pretty(run_dir))
        return

    print("[performance]")
    if perf:
        print(f"profiling_csv={profiling_csv}")
        print(f"rows={perf['num_rows']}")
        print(f"total_wall_s={perf['total_s']:.6f}")
        if perf["sim_time_min"] is not None and perf["sim_time_max"] is not None:
            print(f"sim_time_myr={perf['sim_time_min']:.6f}..{perf['sim_time_max']:.6f}")
        print("timer_sums_s:")
        for name, secs, pct in perf["items"]:
            print(f"  - {name}: {secs:.6f} ({pct:.2f}%)")
    else:
        print("profiling_csv=missing")

    print("")
    print("[energy]")
    if energy:
        print(f"output_h5={output_h5}")
        if "error" in energy:
            print(f"error={energy['error']}")
        else:
            print(f"steps={energy['n_steps']}")
            if energy["time_min"] is not None and energy["time_max"] is not None:
                print(f"time_myr={energy['time_min']:.6f}..{energy['time_max']:.6f}")
            print(f"energy_mean={energy['energy_mean']:.6e}")
            print(f"energy_std={energy['energy_std']:.6e}")
            print(f"energy_min={energy['energy_min']:.6e}")
            print(f"energy_max={energy['energy_max']:.6e}")
            print(f"dE_over_E0_mean={energy['dE_mean']:.6e}")
            print(f"dE_over_E0_std={energy['dE_std']:.6e}")
            print(f"dE_over_E0_min={energy['dE_min']:.6e}")
            print(f"dE_over_E0_max={energy['dE_max']:.6e}")
    else:
        print("output_h5=missing")


if __name__ == "__main__":
    main()
