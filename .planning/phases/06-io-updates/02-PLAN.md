# Plan 02: Add SoA Sync After Checkpoint Restore

## Frontmatter

```yaml
phase: 6
plan: 02
wave: 1
depends_on: []
files_modified:
  - src/read_write.cpp
autonomous: true
```

## Objective

Add SoA sync at the end of `readCheckpoint()` function so that when a simulation is restored from checkpoint, the SoA container is populated from the restored AoS data.

## Context

- `readCheckpoint()` exists in `src/read_write.cpp` but is not yet called anywhere
- The function reads checkpoint data into the AoS `particles[]` array
- After restore, SoA needs to be synced to match the restored state
- This is infrastructure for when checkpoint restore is enabled

## Tasks

<task id="1">
Read `src/read_write.cpp` around line 796 to understand the `readCheckpoint()` function structure.
</task>

<task id="2">
Add include for `particle_data.h` at top of file if not already present.
</task>

<task id="3">
Add SoA sync call at the end of `readCheckpoint()`, just before the success return:

```cpp
// Sync restored AoS data to SoA
particle_data.sync_all_from_particles(particles, last_particle_index + 1);

std::cout << "Checkpoint loaded successfully: " << saveCount << " particles restored." << std::endl;
```

The sync should happen after all particles are restored but before printing the success message.
</task>

<task id="4">
Verify the code compiles without errors.
</task>

## Verification

- [ ] `readCheckpoint()` calls `particle_data.sync_all_from_particles()` after restoring particles
- [ ] Code compiles successfully
- [ ] Sync happens before the function returns

## must_haves

- [ ] Checkpoint restore populates SoA from restored AoS data
- [ ] Prepares for future checkpoint/restart functionality

---
*Generated: 2026-01-17*
