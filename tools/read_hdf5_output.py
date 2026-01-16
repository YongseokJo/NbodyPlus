#!/usr/bin/env python3
"""
Example script for reading ABYSS HDF5 output files.

This script demonstrates how to:
1. Open and explore HDF5 output files
2. Extract particle data at different timesteps
3. Compute basic statistics
4. Create simple visualizations

Usage:
    python read_hdf5_output.py <output.h5>
"""

import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def list_timesteps(filename):
    """List all available timesteps in the HDF5 file."""
    with h5py.File(filename, 'r') as f:
        steps = sorted([key for key in f.keys() if key.startswith('Step_')])
        print(f"\nFound {len(steps)} timesteps:")
        for step in steps:
            time = f[step].attrs['Time_Myr']
            npart = f[step].attrs['NumberOfParticle']
            print(f"  {step}: t={time:.4f} Myr, N={npart} particles")
        return steps


def read_step(filename, step_name):
    """Read all data from a specific timestep."""
    with h5py.File(filename, 'r') as f:
        step = f[step_name]

        # Read metadata
        metadata = {
            'time': step.attrs['Time_Myr'],
            'n_particles': step.attrs['NumberOfParticle'],
            'E_binary': step.attrs['E_binary'],
            'E_binary_SD': step.attrs['E_binary_SD'],
            'E_merger': step.attrs['E_merger'],
            'E_PN': step.attrs['E_PN']
        }

        # Read particle data
        data = {
            'PID': step['PID'][:],
            'mass': step['Mass_Msun'][:],
            'x': step['X_pc'][:],
            'y': step['Y_pc'][:],
            'z': step['Z_pc'][:],
            'vx': step['Vx_km_s'][:],
            'vy': step['Vy_km_s'][:],
            'vz': step['Vz_km_s'][:]
        }

        # Read stellar type if available (SEVN)
        if 'Type' in step:
            data['type'] = step['Type'][:]

        return metadata, data


def compute_statistics(data):
    """Compute basic statistics of particle data."""
    r = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
    v = np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)

    stats = {
        'total_mass': np.sum(data['mass']),
        'mean_mass': np.mean(data['mass']),
        'position_rms': np.sqrt(np.mean(r**2)),
        'velocity_rms': np.sqrt(np.mean(v**2)),
        'position_range': (r.min(), r.max()),
        'velocity_range': (v.min(), v.max())
    }

    return stats


def plot_projection(data, metadata, save=False, output_dir='.'):
    """Create 2D projections of particle positions."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Plot settings
    size = data['mass'] / 10  # Scale marker size by mass
    alpha = 0.5

    # XY projection
    axes[0].scatter(data['x'], data['y'], s=size, alpha=alpha, c='blue')
    axes[0].set_xlabel('X (pc)')
    axes[0].set_ylabel('Y (pc)')
    axes[0].set_title('XY Projection')
    axes[0].axis('equal')
    axes[0].grid(True, alpha=0.3)

    # XZ projection
    axes[1].scatter(data['x'], data['z'], s=size, alpha=alpha, c='green')
    axes[1].set_xlabel('X (pc)')
    axes[1].set_ylabel('Z (pc)')
    axes[1].set_title('XZ Projection')
    axes[1].axis('equal')
    axes[1].grid(True, alpha=0.3)

    # YZ projection
    axes[2].scatter(data['y'], data['z'], s=size, alpha=alpha, c='red')
    axes[2].set_xlabel('Y (pc)')
    axes[2].set_ylabel('Z (pc)')
    axes[2].set_title('YZ Projection')
    axes[2].axis('equal')
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(f"Particle Distribution at t={metadata['time']:.4f} Myr",
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save:
        output_file = Path(output_dir) / f"projection_t{metadata['time']:.4f}.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
    else:
        plt.show()


def plot_velocity_distribution(data, metadata, save=False, output_dir='.'):
    """Plot velocity distribution."""
    v = np.sqrt(data['vx']**2 + data['vy']**2 + data['vz']**2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Histogram
    axes[0].hist(v, bins=50, alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('Velocity (km/s)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Velocity Distribution')
    axes[0].axvline(np.mean(v), color='r', linestyle='--',
                    label=f'Mean: {np.mean(v):.2f} km/s')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Velocity components
    axes[1].hist(data['vx'], bins=30, alpha=0.5, label='Vx', edgecolor='black')
    axes[1].hist(data['vy'], bins=30, alpha=0.5, label='Vy', edgecolor='black')
    axes[1].hist(data['vz'], bins=30, alpha=0.5, label='Vz', edgecolor='black')
    axes[1].set_xlabel('Velocity (km/s)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Velocity Components')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.suptitle(f"Velocity Analysis at t={metadata['time']:.4f} Myr",
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save:
        output_file = Path(output_dir) / f"velocity_t{metadata['time']:.4f}.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
    else:
        plt.show()


def main():
    if len(sys.argv) < 2:
        print("Usage: python read_hdf5_output.py <output.h5>")
        sys.exit(1)

    filename = sys.argv[1]

    if not Path(filename).exists():
        print(f"Error: File '{filename}' not found!")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Reading ABYSS HDF5 output: {filename}")
    print(f"{'='*60}")

    # List all timesteps
    steps = list_timesteps(filename)

    if not steps:
        print("No timesteps found in file!")
        sys.exit(1)

    # Read first and last timestep
    print(f"\n{'='*60}")
    print("First timestep:")
    print(f"{'='*60}")
    metadata_0, data_0 = read_step(filename, steps[0])
    print(f"Time: {metadata_0['time']:.4f} Myr")
    print(f"Number of particles: {metadata_0['n_particles']}")
    print(f"E_binary: {metadata_0['E_binary']:.6e}")
    print(f"E_merger: {metadata_0['E_merger']:.6e}")

    stats_0 = compute_statistics(data_0)
    print(f"\nStatistics:")
    print(f"  Total mass: {stats_0['total_mass']:.2f} Msun")
    print(f"  Mean mass: {stats_0['mean_mass']:.2f} Msun")
    print(f"  Position RMS: {stats_0['position_rms']:.4f} pc")
    print(f"  Velocity RMS: {stats_0['velocity_rms']:.2f} km/s")

    if len(steps) > 1:
        print(f"\n{'='*60}")
        print("Last timestep:")
        print(f"{'='*60}")
        metadata_n, data_n = read_step(filename, steps[-1])
        print(f"Time: {metadata_n['time']:.4f} Myr")
        print(f"Number of particles: {metadata_n['n_particles']}")
        print(f"E_binary: {metadata_n['E_binary']:.6e}")
        print(f"E_merger: {metadata_n['E_merger']:.6e}")

        stats_n = compute_statistics(data_n)
        print(f"\nStatistics:")
        print(f"  Total mass: {stats_n['total_mass']:.2f} Msun")
        print(f"  Mean mass: {stats_n['mean_mass']:.2f} Msun")
        print(f"  Position RMS: {stats_n['position_rms']:.4f} pc")
        print(f"  Velocity RMS: {stats_n['velocity_rms']:.2f} km/s")

    # Create plots
    print(f"\n{'='*60}")
    print("Creating visualizations...")
    print(f"{'='*60}")

    plot_projection(data_0, metadata_0, save=False)
    plot_velocity_distribution(data_0, metadata_0, save=False)

    print("\nDone!")


if __name__ == "__main__":
    main()
