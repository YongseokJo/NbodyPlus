"""
Analyze energy conservation for ABYSS HDF5 output files.

This script reads HDF5 output from ABYSS simulations and computes:
- Kinetic, potential, and total energy at each timestep
- Energy conservation (dE/E0) over time
- Creates plots of energy evolution

Usage:
    python analyze_energy.py <output.h5> [--plot <path>] [--csv <path>]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import h5py

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    sys.stderr.write("ERROR: matplotlib is required to plot results.\n")
    raise

# Optional high-performance backends (priority: C++ > numba > scipy > numpy)
HAS_CPP = False
HAS_NUMBA = False
HAS_SCIPY = False

try:
    from potential_energy import compute_potential_energy as _potential_cpp
    HAS_CPP = True
except ImportError:
    pass

if not HAS_CPP:
    try:
        from numba import njit, prange
        HAS_NUMBA = True
    except ImportError:
        pass

if not HAS_CPP and not HAS_NUMBA:
    try:
        from scipy.spatial.distance import pdist
        HAS_SCIPY = True
    except ImportError:
        pass

# Units from src/def.h
POSITION_UNIT_PC = 4.0
TIME_UNIT_YR = 1e10
VELOCITY_UNIT_PC_YR = 4e-10
MASS_UNIT_MSUN = 0.0001424198

PC_IN_KM = 3.08567758149137e13
SEC_PER_YR = 3.1536e7
KM_S_TO_PC_YR = SEC_PER_YR / PC_IN_KM


# -----------------------------------------------------------------------------
# High-performance potential energy backends
# -----------------------------------------------------------------------------

if HAS_NUMBA:
    @njit(parallel=True, fastmath=True)
    def _potential_numba(pos, m):
        """Numba-accelerated potential energy computation with parallel loops."""
        n = pos.shape[0]
        potential = 0.0
        # Use parallel reduction over chunks
        for i in prange(n - 1):
            local_sum = 0.0
            xi, yi, zi = pos[i, 0], pos[i, 1], pos[i, 2]
            mi = m[i]
            for j in range(i + 1, n):
                dx = xi - pos[j, 0]
                dy = yi - pos[j, 1]
                dz = zi - pos[j, 2]
                r = np.sqrt(dx * dx + dy * dy + dz * dz)
                if r > 0.0:
                    local_sum += mi * m[j] / r
            potential += local_sum
        return -potential


def _potential_scipy(pos, m):
    """Scipy-accelerated potential energy using pdist."""
    n = pos.shape[0]
    # pdist returns condensed distance matrix (upper triangle, row-major)
    distances = pdist(pos, metric='euclidean')
    # Build mass products for upper triangle pairs
    # For pdist, the index mapping is: pair (i,j) where i<j maps to
    # index = n*i - i*(i+1)/2 + j - i - 1
    # But we can use broadcasting: for all pairs (i,j) with i<j
    iu = np.triu_indices(n, k=1)
    mass_products = m[iu[0]] * m[iu[1]]
    inv_r = np.where(distances > 0.0, 1.0 / distances, 0.0)
    return -np.sum(mass_products * inv_r)


def _potential_numpy(pos, m):
    """Standard numpy computation (fallback)."""
    n = pos.shape[0]
    if n > 5000:
        # Chunked computation for memory efficiency
        potential = 0.0
        chunk_size = 1000
        for i in range(0, n, chunk_size):
            i_end = min(i + chunk_size, n)
            for j in range(i, n, chunk_size):
                j_end = min(j + chunk_size, n)
                diff = pos[i:i_end, None, :] - pos[None, j:j_end, :]
                r = np.linalg.norm(diff, axis=2)
                if i == j:
                    iu = np.triu_indices(i_end - i, k=1)
                    rij = r[iu]
                    inv_r = np.where(rij > 0.0, 1.0 / rij, 0.0)
                    potential -= np.sum(m[i:i_end][iu[0]] * m[i:i_end][iu[1]] * inv_r)
                else:
                    inv_r = np.where(r > 0.0, 1.0 / r, 0.0)
                    m_i = m[i:i_end][:, None]
                    m_j = m[j:j_end][None, :]
                    potential -= np.sum(m_i * m_j * inv_r)
    else:
        diff = pos[:, None, :] - pos[None, :, :]
        r = np.linalg.norm(diff, axis=2)
        iu = np.triu_indices(n, k=1)
        rij = r[iu]
        inv_r = np.where(rij > 0.0, 1.0 / rij, 0.0)
        potential = -np.sum(m[iu[0]] * m[iu[1]] * inv_r)
    return potential


def compute_potential(pos, m):
    """Compute potential energy using the fastest available backend."""
    if HAS_CPP:
        return _potential_cpp(pos, m)
    elif HAS_NUMBA:
        return _potential_numba(pos, m)
    elif HAS_SCIPY:
        return _potential_scipy(pos, m)
    else:
        return _potential_numpy(pos, m)


def read_hdf5_step(filename, step_name):
    """Read particle data from a specific timestep in HDF5 file."""
    with h5py.File(filename, 'r') as f:
        step = f[step_name]

        # Read metadata (convert to Python float to avoid numpy scalar issues)
        time_myr = float(step.attrs['Time_Myr'])

        # Read particle data
        masses = step['Mass_Msun'][:]
        pos = np.column_stack([
            step['X_pc'][:],
            step['Y_pc'][:],
            step['Z_pc'][:]
        ])
        vel = np.column_stack([
            step['Vx_km_s'][:],
            step['Vy_km_s'][:],
            step['Vz_km_s'][:]
        ])

        # Read energy tracking attributes if available
        energy_attrs = {}
        for attr in ['E_binary', 'E_binary_SD', 'E_merger', 'E_PN']:
            if attr in step.attrs:
                energy_attrs[attr] = step.attrs[attr]

    return time_myr, masses, pos, vel, energy_attrs


def get_timesteps(filename):
    """Get sorted list of timestep names from HDF5 file."""
    with h5py.File(filename, 'r') as f:
        steps = [key for key in f.keys() if key.startswith('Step_')]
        # Sort by step number
        steps.sort(key=lambda x: int(x.split('_')[1]))
    return steps


def compute_energy_code_units(m_msun, pos_pc, vel_kms):
    """
    Compute kinetic, potential, and total energy in code units.

    Parameters:
        m_msun: masses in solar masses
        pos_pc: positions in parsecs (N x 3)
        vel_kms: velocities in km/s (N x 3)

    Returns:
        kinetic, potential, total energy in code units
    """
    # Convert to code units
    m_code = m_msun / MASS_UNIT_MSUN
    pos_code = pos_pc / POSITION_UNIT_PC
    vel_pcyr = vel_kms * KM_S_TO_PC_YR
    vel_code = vel_pcyr / VELOCITY_UNIT_PC_YR

    # Kinetic energy: sum of 0.5 * m * v^2
    kinetic = 0.5 * np.sum(m_code * np.sum(vel_code ** 2, axis=1))

    # Potential energy: -G * sum(mi * mj / rij) for all pairs
    # In code units, G = 1
    # Use optimized backend (Numba > scipy > numpy fallback)
    potential = compute_potential(pos_code, m_code)

    return kinetic, potential, kinetic + potential


def main():
    parser = argparse.ArgumentParser(
        description="Analyze energy conservation for ABYSS HDF5 output files."
    )
    parser.add_argument(
        "hdf5_file",
        help="HDF5 output file from ABYSS simulation.",
    )
    parser.add_argument(
        "--plot",
        default=None,
        help="Path to save plot (default: <input_dir>/energy_analysis.png).",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Path to save CSV table (default: <input_dir>/energy_analysis.csv).",
    )
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Skip generating CSV output.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip generating plot (useful for headless environments).",
    )
    args = parser.parse_args()

    hdf5_path = Path(args.hdf5_file)
    if not hdf5_path.exists():
        raise SystemExit(f"Error: File '{args.hdf5_file}' not found!")

    print(f"\n{'='*60}")
    print(f"Analyzing energy conservation: {hdf5_path.name}")
    print(f"{'='*60}")
    if HAS_CPP:
        print(f"Backend: C++ OpenMP (potential_energy module)")
    elif HAS_NUMBA:
        print(f"Backend: Numba (parallel JIT)")
    elif HAS_SCIPY:
        print(f"Backend: scipy.spatial.distance.pdist")
    else:
        print(f"Backend: numpy (install numba/scipy or build C++ module for better performance)")
    print()

    # Get all timesteps
    steps = get_timesteps(args.hdf5_file)
    if not steps:
        raise SystemExit(f"No timesteps found in {args.hdf5_file}")

    print(f"Found {len(steps)} timesteps")

    times = []
    kinetic_list = []
    potential_list = []
    total_list = []
    e_binary_list = []
    e_merger_list = []

    for i, step in enumerate(steps):
        time_myr, masses, pos, vel, energy_attrs = read_hdf5_step(args.hdf5_file, step)
        kinetic, potential, total = compute_energy_code_units(masses, pos, vel)

        times.append(time_myr)
        kinetic_list.append(kinetic)
        potential_list.append(potential)
        total_list.append(total)
        e_binary_list.append(energy_attrs.get('E_binary', 0.0))
        e_merger_list.append(energy_attrs.get('E_merger', 0.0))

        if (i + 1) % 10 == 0 or i == len(steps) - 1:
            print(f"  Processed {i + 1}/{len(steps)} timesteps...")

    times = np.array(times)
    kinetic_list = np.array(kinetic_list)
    potential_list = np.array(potential_list)
    total_list = np.array(total_list)
    e_binary_list = np.array(e_binary_list)
    e_merger_list = np.array(e_merger_list)

    # Compute energy residual
    e0 = total_list[0]
    residual = (total_list - e0) / abs(e0)
    residual_abs = np.abs(residual)
    residual_plot = np.where(residual_abs > 0.0, residual_abs, np.nan)

    # Print results table
    print(f"\n{'='*60}")
    print("Energy Analysis Results")
    print(f"{'='*60}\n")

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

    # Summary statistics
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    print(f"Initial energy (E0): {e0:.8e}")
    print(f"Final energy:        {total_list[-1]:.8e}")
    print(f"Max |dE/E0|:         {np.nanmax(residual_abs):.6e}")
    print(f"Mean |dE/E0|:        {np.nanmean(residual_abs):.6e}")
    print(f"Simulation time:     {times[0]:.4f} - {times[-1]:.4f} Myr")

    # Save CSV by default (unless --no-csv is specified)
    if not args.no_csv:
        csv_path = args.csv or str(hdf5_path.parent / "energy_analysis.csv")
        with open(csv_path, "w") as handle:
            handle.write("idx,time_myr,kinetic,potential,total,dE_over_E0,E_binary,E_merger\n")
            for idx, (t, k, u, e, de, eb, em) in enumerate(
                zip(times, kinetic_list, potential_list, total_list, residual,
                    e_binary_list, e_merger_list)
            ):
                handle.write(f"{idx},{t},{k},{u},{e},{de},{eb},{em}\n")
        print(f"\nCSV saved to: {csv_path}")

    # Create plot
    if not args.no_plot:
        plot_path = args.plot or str(hdf5_path.parent / "energy_analysis.png")

        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Energy plot
        ax_energy = axes[0]
        ax_energy.plot(times, kinetic_list, label="Kinetic", linewidth=1.5)
        ax_energy.plot(times, potential_list, label="Potential", linewidth=1.5)
        ax_energy.plot(times, total_list, label="Total", linewidth=2, color='black')
        ax_energy.set_ylabel("Energy (code units)")
        ax_energy.set_title("Energy Evolution")
        ax_energy.grid(True, alpha=0.3)
        ax_energy.legend(loc='best')

        # Residual plot
        ax_residual = axes[1]
        ax_residual.plot(times, residual_plot, label="|dE/E0|", color='red', linewidth=1.5)
        ax_residual.set_yscale("log")
        ax_residual.set_xlabel("Time (Myr)")
        ax_residual.set_ylabel("|dE/E0|")
        ax_residual.set_title("Energy Conservation")
        ax_residual.grid(True, which="both", alpha=0.3)
        ax_residual.legend(loc='best')

        fig.suptitle(f"ABYSS Energy Analysis: {hdf5_path.name}", fontsize=12, fontweight='bold')
        fig.tight_layout()
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {plot_path}")
        plt.close(fig)

    print("\nDone!")


if __name__ == "__main__":
    main()
