#!/usr/bin/env python3
"""
Verify energy of McLuster-generated initial conditions.

Computes virial ratio Q = 2K/|U| where:
- K = kinetic energy = 0.5 * sum(m * v^2)
- U = potential energy = -G * sum(mi*mj/rij)

For virial equilibrium: Q = 1.0 (McLuster default Q=0.5 means 2K/|U|=1)

Usage:
    python verify_energy.py <ic_file> [tolerance]

Exit codes:
    0 = PASS (|Q - 1.0| < tolerance)
    1 = FAIL (|Q - 1.0| >= tolerance)
"""

import sys
import numpy as np

# Unit definitions from src/def.h
POSITION_UNIT_PC = 4.0       # position_unit = 4 pc
MASS_UNIT_MSUN = 0.0001424198  # mass_unit in Msun where G=1

# Conversion factors
PC_IN_KM = 3.08567758149137e13
SEC_PER_YR = 3.1536e7
KM_S_TO_PC_YR = SEC_PER_YR / PC_IN_KM  # km/s -> pc/yr


def compute_virial_ratio(filename):
    """
    Compute 2K/|U| for particles in ABYSS IC file.

    ABYSS IC format: x y z vx vy vz mass
    Units in file: kpc, km/s, 1e-9 Msun

    Returns: (kinetic, potential, virial_ratio)
    """
    data = np.loadtxt(filename)
    n = len(data)

    # ABYSS format: x y z vx vy vz mass (columns 0-6)
    pos_kpc = data[:, 0:3]      # kpc
    vel_kms = data[:, 3:6]      # km/s
    mass_1e9 = data[:, 6]       # 1e-9 Msun units

    # Convert to code units where G=1
    # mass: 1e-9 Msun -> Msun -> code units
    mass = mass_1e9 * 1e9 / MASS_UNIT_MSUN

    # position: kpc -> pc -> code units
    pos = pos_kpc * 1000.0 / POSITION_UNIT_PC

    # velocity: km/s -> pc/yr -> code units
    # velocity_unit = 4e-10 pc/yr, so code = vel_pc_yr / 4e-10
    vel_pc_yr = vel_kms * KM_S_TO_PC_YR
    VELOCITY_UNIT = 4e-10  # from def.h
    vel = vel_pc_yr / VELOCITY_UNIT

    # Kinetic energy: 0.5 * sum(m * v^2)
    v_squared = np.sum(vel**2, axis=1)
    kinetic = 0.5 * np.sum(mass * v_squared)

    # Potential energy: -sum(mi*mj/rij) for i<j
    # Use vectorized upper triangle calculation
    potential = 0.0
    for i in range(n - 1):
        dr = pos[i] - pos[i+1:]  # shape: (n-i-1, 3)
        r = np.sqrt(np.sum(dr**2, axis=1))
        # Avoid division by zero
        r = np.where(r > 0, r, 1e-10)
        potential -= np.sum(mass[i] * mass[i+1:] / r)

    # Virial ratio: 2K/|U|
    if abs(potential) < 1e-20:
        return kinetic, potential, float('inf')

    virial_ratio = 2.0 * kinetic / abs(potential)
    return kinetic, potential, virial_ratio


def main():
    if len(sys.argv) < 2:
        print("Usage: python verify_energy.py <ic_file> [tolerance]")
        sys.exit(1)

    ic_file = sys.argv[1]
    tolerance = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4

    K, U, Q = compute_virial_ratio(ic_file)

    print(f"Kinetic energy:   {K:.10e}")
    print(f"Potential energy: {U:.10e}")
    print(f"Virial ratio (2K/|U|): {Q:.10f}")
    print(f"Expected: 1.0 (tolerance: {tolerance})")

    deviation = abs(Q - 1.0)
    if deviation > tolerance:
        print(f"FAIL: |Q - 1.0| = {deviation:.6e} > {tolerance}")
        sys.exit(1)

    print("PASS: Energy check within tolerance")
    sys.exit(0)


if __name__ == "__main__":
    main()
