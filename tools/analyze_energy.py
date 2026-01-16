#!/usr/bin/env python3
import argparse
import glob
import os
import re
import sys

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    sys.stderr.write("ERROR: matplotlib is required to plot results.\n")
    raise

# Units from src/def.h
POSITION_UNIT_PC = 4.0
TIME_UNIT_YR = 1e10
VELOCITY_UNIT_PC_YR = 4e-10
MASS_UNIT_MSUN = 0.0001424198

PC_IN_KM = 3.08567758149137e13
SEC_PER_YR = 3.1536e7
KM_S_TO_PC_YR = SEC_PER_YR / PC_IN_KM


def parse_snapshot(path):
    time_myr = None
    masses = []
    pos = []
    vel = []
    header_found = False

    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Time"):
                match = re.search(r"Time\s*=\s*([0-9eE+\-.]+)", line)
                if match:
                    time_myr = float(match.group(1))
                continue
            if line.startswith("PID"):
                header_found = True
                continue
            if not header_found:
                continue

            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                m_msun = float(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
                vx = float(parts[5])
                vy = float(parts[6])
                vz = float(parts[7])
            except ValueError:
                continue
            masses.append(m_msun)
            pos.append((x, y, z))
            vel.append((vx, vy, vz))

    if time_myr is None:
        raise ValueError(f"Missing time header in {path}")
    if not masses:
        raise ValueError(f"No particle data parsed in {path}")

    return (
        time_myr,
        np.array(masses, dtype=np.float64),
        np.array(pos, dtype=np.float64),
        np.array(vel, dtype=np.float64),
    )


def compute_energy_code_units(m_msun, pos_pc, vel_kms):
    m_code = m_msun / MASS_UNIT_MSUN
    pos_code = pos_pc / POSITION_UNIT_PC
    vel_pcyr = vel_kms * KM_S_TO_PC_YR
    vel_code = vel_pcyr / VELOCITY_UNIT_PC_YR

    kinetic = 0.5 * np.sum(m_code * np.sum(vel_code ** 2, axis=1))

    n = pos_code.shape[0]
    diff = pos_code[:, None, :] - pos_code[None, :, :]
    r = np.linalg.norm(diff, axis=2)
    iu = np.triu_indices(n, k=1)
    rij = r[iu]
    inv_r = np.where(rij > 0.0, 1.0 / rij, 0.0)
    potential = -np.sum(m_code[iu[0]] * m_code[iu[1]] * inv_r)

    return kinetic, potential, kinetic + potential


def main():
    parser = argparse.ArgumentParser(
        description="Analyze energy conservation for ABYSS output files."
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        default="test/test1/output",
        help="Directory containing output_*.txt files.",
    )
    parser.add_argument(
        "--plot",
        default=None,
        help="Path to save plot (default: <output_dir>/energy_analysis.png).",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Optional path to save CSV table.",
    )
    args = parser.parse_args()

    pattern = os.path.join(args.output_dir, "output_*.txt")
    files = glob.glob(pattern)
    if not files:
        raise SystemExit(f"No output files found at {pattern}")

    def file_key(path):
        match = re.search(r"output_(\d+)\.txt", os.path.basename(path))
        return int(match.group(1)) if match else -1

    files = sorted(files, key=file_key)

    times = []
    kinetic_list = []
    potential_list = []
    total_list = []

    for path in files:
        time_myr, masses, pos, vel = parse_snapshot(path)
        kinetic, potential, total = compute_energy_code_units(masses, pos, vel)
        times.append(time_myr)
        kinetic_list.append(kinetic)
        potential_list.append(potential)
        total_list.append(total)

    times = np.array(times)
    kinetic_list = np.array(kinetic_list)
    potential_list = np.array(potential_list)
    total_list = np.array(total_list)

    e0 = total_list[0]
    residual = (total_list - e0) / abs(e0)
    residual_abs = np.abs(residual)
    residual_plot = np.where(residual_abs > 0.0, residual_abs, np.nan)

    header = (
        f"{'idx':>4} {'time_myr':>12} {'kinetic':>18} {'potential':>18} "
        f"{'total':>18} {'dE/E0':>14}"
    )
    print(header)
    print("-" * len(header))
    for idx, (t, k, u, e, de) in enumerate(
        zip(times, kinetic_list, potential_list, total_list, residual)
    ):
        print(f"{idx:4d} {t:12.6f} {k:18.8e} {u:18.8e} {e:18.8e} {de:14.6e}")

    if args.csv:
        with open(args.csv, "w") as handle:
            handle.write("idx,time_myr,kinetic,potential,total,dE_over_E0\n")
            for idx, (t, k, u, e, de) in enumerate(
                zip(times, kinetic_list, potential_list, total_list, residual)
            ):
                handle.write(f"{idx},{t},{k},{u},{e},{de}\n")

    plot_path = args.plot or os.path.join(args.output_dir, "energy_analysis.png")
    fig, (ax_energy, ax_residual) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    ax_energy.plot(times, kinetic_list, label="Kinetic")
    ax_energy.plot(times, potential_list, label="Potential")
    ax_energy.plot(times, total_list, label="Total")
    ax_energy.set_ylabel("Energy (code units)")
    ax_energy.grid(True, alpha=0.3)
    ax_energy.legend()

    ax_residual.plot(times, residual_plot, label="|dE/E0|")
    ax_residual.set_yscale("log")
    ax_residual.set_xlabel("Time (Myr)")
    ax_residual.set_ylabel("|dE/E0|")
    ax_residual.grid(True, which="both", alpha=0.3)
    ax_residual.legend()

    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    print(f"Plot saved to: {plot_path}")


if __name__ == "__main__":
    main()
