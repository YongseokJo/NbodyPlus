# Plan 04-03 Summary: Update Prediction and Correction Routines

## Status: Complete

## Commits

| Hash | Description |
|------|-------------|
| 5d8bcbe | feat(04-03): add SoA-compatible prediction/correction routines |

## What Was Built

Added SoA-compatible free functions for particle prediction and correction that operate on ParticleData directly.

### Task 1: predict_second_order()
```cpp
void predict_second_order(const ParticleData& data, size_t i, double dt,
                          double pos_out[3], double vel_out[3])
```
- Second-order Taylor expansion prediction using acceleration and jerk
- Extracts position, velocity, acc_total from ParticleData
- Outputs predicted position/velocity to arrays

### Task 2: correct_fourth_order()
```cpp
void correct_fourth_order(ParticleData& data, size_t i, double dt,
                          const double pos[3], const double vel[3],
                          const double a[3][4])
```
- Applies 4th order corrections using snap (a[dim][2]) and crackle (a[dim][3])
- Writes corrected values to new_position/new_velocity in ParticleData

### Task 3: update_particle_state()
```cpp
void update_particle_state(ParticleData& data, size_t i)
```
- Copies new_position/new_velocity to position/velocity
- Finalizes the integration step

### Task 4: update_neighbor_radius()
```cpp
void update_neighbor_radius(ParticleData& data, size_t i, int target_neighbors)
```
- Adjusts neighbor_radius_sq based on neighbor count vs target
- Same polynomial adjustment logic as Particle::update_radius()

### Task 5: Backward compatibility
- Existing Particle:: methods unchanged
- New free functions can be called from either Particle context or pure SoA context

## Files Modified

- `src/Particle/update_particle.cpp` — Added free function implementations
- `src/global.h` — Added function declarations

## Deviations

None. Numerical logic matches exactly the existing Particle:: methods.
