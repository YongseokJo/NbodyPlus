# Plan 04-01 Summary: Add SoA Helper Functions to ParticleData

## Status: Complete

## Commits

| Hash | Description |
|------|-------------|
| 938a7bc | feat(04-01): add SoA helper functions to ParticleData |

## What Was Built

Added helper functions to ParticleData class for bridging between AoS Particle struct and SoA layout:

### Task 1: sync_from_particle()
- Copies all particle fields from Particle struct to ParticleData at index i
- Handles position, velocity, mass, all acceleration arrays, timestep fields, neighbor info, and flags

### Task 2: sync_to_particle()
- Copies all fields from ParticleData at index i back to Particle struct
- Maintains exact numerical precision (no conversions)

### Task 3: Bulk 3-vector accessors
- `get_position_vec()`, `get_velocity_vec()`, `get_new_position_vec()`, `get_new_velocity_vec()`
- `set_position_vec()`, `set_velocity_vec()`, `set_new_position_vec()`, `set_new_velocity_vec()`
- `get_acc_total_array()`, `set_acc_total_array()`, and similar for reg/irr

### Task 4: Acceleration accumulation helpers
- `zero_all_acc()`, `zero_acc_total()`, `zero_acc_reg()`, `zero_acc_irr()`
- `add_to_acc_total()`, `add_to_acc_reg()`, `add_to_acc_irr()`
- `add_to_acc_total_vec()`, `add_to_acc_reg_vec()`, `add_to_acc_irr_vec()`

## Files Modified

- `src/particle_data.h` — Added method declarations and inline implementations
- `src/particle_data.cpp` — Added sync_from_particle and sync_to_particle implementations

## Deviations

None.
