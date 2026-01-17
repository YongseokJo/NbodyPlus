# Plan 04-02 Summary: Update Timestep Routines for SoA Compatibility

## Status: Complete

## Commits

| Hash | Description |
|------|-------------|
| 91066de | feat(04-02): add SoA-compatible timestep routines |

## What Was Built

Added SoA-compatible overloads for timestep calculation routines that accept ParticleData& and index instead of requiring Particle struct access.

### Task 1: Overloaded versions with ParticleData
- `getNewTimeStepReg(const ParticleData& data, size_t i)`
- `getNewTimeStepIrr(const ParticleData& data, size_t i)`
- `getNewTimeStep(const ParticleData& data, size_t i)`

### Task 2: Free function helpers
- `extract_velocity()` - Extract velocity array from ParticleData
- `extract_acc_total()` - Extract total acceleration array
- `extract_acc_reg()` - Extract regular acceleration array
- `extract_acc_irr()` - Extract irregular acceleration array

### Task 3: Backward compatibility
- Existing Particle:: methods (calculate_time_step_irr, calculate_time_step_reg) unchanged
- They continue to call the raw array versions directly
- New SoA versions delegate to existing raw array functions

### Task 4: Header consolidation
- Moved timestep function declarations to global.h
- Removed local declarations from update_particle.cpp

## Files Modified

- `src/timestep_routines.cpp` — Added SoA overloads and extract helpers
- `src/global.h` — Added timestep function declarations
- `src/Particle/update_particle.cpp` — Removed local declarations

## Deviations

None. The existing raw array interface was already suitable for SoA compatibility; we just added overloads that extract arrays from ParticleData.
