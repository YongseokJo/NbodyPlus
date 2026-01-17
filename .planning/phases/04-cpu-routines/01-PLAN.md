# Plan 04-01: Add SoA Helper Functions to ParticleData

## Frontmatter

```yaml
wave: 1
depends_on: []
files_modified:
  - src/particle_data.h
  - src/particle_data.cpp
autonomous: true
```

## Objective

Add helper functions to ParticleData that facilitate bulk operations and SoA-compatible force calculations. These will be used by subsequent plans to bridge between the existing AoS Particle struct and SoA operations.

## Context

The current ParticleData class has individual field accessors but lacks methods for:
1. Copying data from/to Particle struct
2. Bulk operations on acceleration arrays
3. Convenience methods for force calculation patterns

## Tasks

<task id="1">
Add sync_from_particle() method to copy data from a Particle struct to ParticleData at index i
</task>

<task id="2">
Add sync_to_particle() method to copy data from ParticleData at index i back to Particle struct
</task>

<task id="3">
Add bulk accessor methods for getting position/velocity as 3-vectors
</task>

<task id="4">
Add acceleration accumulation helpers (add_to_acc_*, zero_acc_*)
</task>

## Verification

- [ ] New methods compile without errors
- [ ] sync_from_particle followed by sync_to_particle preserves all field values
- [ ] Acceleration helpers correctly accumulate values

## must_haves

- sync_from_particle() copies all fields needed for force calculation
- sync_to_particle() copies results back after force calculation
- Methods maintain exact numerical values (no precision loss)
