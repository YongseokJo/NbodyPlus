# Plan 01: Add SoA Sync After Initialization

## Frontmatter

```yaml
phase: 6
plan: 01
wave: 1
depends_on: []
files_modified:
  - src/root_routines.cpp
autonomous: true
```

## Objective

After `InitializationRoutines()` completes, sync all particle data from the AoS `particles[]` array to the SoA `particle_data` container. This ensures the SoA is populated with initial particle state before the simulation loop begins.

## Context

- `InitializationRoutines()` is called in `RootRoutines()` at line 78 of `src/root_routines.cpp`
- After initialization, particles have their positions, velocities, masses, accelerations, timesteps, and neighbor info set
- The SoA container needs this data for worker processes that use MPI shared memory
- The `sync_all_from_particles()` method exists in `ParticleData` class (added in Phase 5)

## Tasks

<task id="1">
Read `src/root_routines.cpp` to confirm the location of `InitializationRoutines()` call and understand the surrounding code.
</task>

<task id="2">
Add include for `particle_data.h` if not already present at top of file.
</task>

<task id="3">
Add SoA sync call immediately after `InitializationRoutines()` returns:

```cpp
/* Initialization */
InitializationRoutines(queue_scheduler, workers);

// Sync AoS to SoA after initialization completes
particle_data.sync_all_from_particles(particles, last_particle_index + 1);
```
</task>

<task id="4">
Verify the code compiles without errors.
</task>

## Verification

- [ ] `particle_data.sync_all_from_particles()` is called after `InitializationRoutines()`
- [ ] Code compiles successfully with `build.sh`
- [ ] No runtime errors during initialization

## must_haves

- [ ] SoA is populated from AoS after particle initialization
- [ ] Satisfies IO-02: "Reading from input files populates SoA correctly"

---
*Generated: 2026-01-17*
