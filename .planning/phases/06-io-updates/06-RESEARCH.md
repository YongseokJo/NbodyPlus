# Phase 6: I/O Updates — Research

## Summary

Phase 6 requires updating file I/O to work with the SoA layout. The codebase has two main I/O locations:
1. `src/read_write.cpp` — HDF5 output and input reading
2. `src/initialization_routines.cpp` — Particle initialization after reading

## Current I/O Implementation Analysis

### Input: `readData()` in `read_write.cpp`

**Current flow:**
1. Opens input file (text format, 7 columns: mass, x, y, z, vx, vy, vz)
2. Reads into temporary `data[][]` array
3. Calls `particles[row].initialize(data[row], row)` — AoS direct write
4. Normalizes particles via `particles[i].normalize_particle()`

**SoA change needed:**
- After `particles[].initialize()` populates AoS, sync to SoA via `particle_data.sync_from_particle()`
- Or populate SoA directly and sync to AoS (but AoS still needed for SDAR compatibility)

**Recommended approach:** Keep AoS initialization, add SoA sync at end.

### Output: `writeParticle()` in `read_write.cpp`

**Current flow:**
1. Iterates over `particles[]` array (AoS)
2. Uses `ptcl->predict_particle_second_order()` to get predicted pos/vel
3. Accesses `ptcl->pid`, `ptcl->mass`, etc. directly
4. Collects group members from `particles[ptcl->members[j]]`
5. Writes to HDF5 datasets: PID, Mass_Msun, X_pc, Y_pc, Z_pc, Vx_km_s, Vy_km_s, Vz_km_s

**SoA consideration:**
- Output uses predicted values, not raw stored values
- The prediction method `predict_particle_second_order()` is on the AoS Particle struct
- Current code already collects into std::vector before HDF5 write

**Decision:** Output can stay as-is since it reads from AoS (which is synced from SoA via exit-only pattern). No changes strictly required for IO-01.

### Checkpoint: `writeCheckpoint()` and `readCheckpoint()`

**Current flow:**
- Writes all particle fields to HDF5 for restart
- Reads back and restores particle state

**SoA consideration:**
- Checkpoint writes from AoS `particles[]`
- Checkpoint reads into AoS `particles[]`
- After read, need to sync to SoA

**Recommended:** Add SoA sync after `readCheckpoint()`.

### Initialization: `InitializationRoutines()` in `initialization_routines.cpp`

**Current flow:**
- Called after `readData()`
- Uses `particles[]` AoS throughout
- Runs various TASK_INIT_* to compute accelerations
- Time step initialization uses `particles[i]` directly

**SoA consideration:**
- Initialization tasks use worker processes that access MPI shared memory
- MPI shared memory is now SoA (`ParticleDataMPI`)
- The `particles[]` AoS is still used by root for orchestration

**Decision:** After initialization completes, ensure SoA is populated. The worker tasks already use SoA via the bridge pattern from Phase 4.

## Key Findings

### What needs to change for IO-01 (HDF5 output from SoA):

1. **Option A (minimal):** Keep current output code unchanged. AoS is authoritative for FewBody operations and gets synced. Output reads from AoS which is always valid.

2. **Option B (SoA-native):** Modify `writeParticle()` to read directly from SoA:
   - Replace `ptcl->mass` with `particle_data.get_mass(i)`
   - Replace `ptcl->pid` with `particle_data.get_pid(i)`
   - Handle prediction using SoA fields

**Recommendation:** Option A is safer. Output already works correctly since AoS is maintained.

### What needs to change for IO-02 (Input file reading populates SoA):

1. After `readData()` completes, call `particle_data.sync_all_from_particles(particles, num_particles)` to populate SoA from the initialized AoS data.

2. After `readCheckpoint()` completes, call the same sync function.

3. This should happen in `main.cpp` or wherever these functions are called from.

## File Locations

| File | Current State | Changes Needed |
|------|---------------|----------------|
| `src/read_write.cpp` | Uses AoS for output | Minimal (add SoA sync after readData/readCheckpoint) |
| `src/initialization_routines.cpp` | Uses AoS | None (workers already use SoA via bridge) |
| `src/main.cpp` | Calls readData() | Add SoA sync call after initialization |

## Dependencies

- `ParticleData::sync_all_from_particles()` — added in Phase 5
- Global `particle_data` instance — from Phase 2

## Implementation Approach

### Plan 1: Add SoA sync after input reading

**Files:** `src/main.cpp` or wherever readData() is called
**Changes:**
1. After `readData()` returns successfully
2. After `InitializationRoutines()` completes
3. Call `particle_data.sync_all_from_particles(particles, last_particle_index + 1)`

### Plan 2: Add SoA sync after checkpoint reading

**Files:** `src/main.cpp` or wherever checkpoint is restored
**Changes:**
1. After `readCheckpoint()` returns true
2. Call `particle_data.sync_all_from_particles(particles, last_particle_index + 1)`

## Verification

- IO-01: Run simulation, check HDF5 output contains correct data
- IO-02: Run simulation from fresh start, verify SoA is populated correctly
- Checkpoint: Save and restore checkpoint, verify simulation continues correctly

## Risks

- **Low risk:** Changes are additive sync calls, not modifications to existing I/O logic
- **Compatibility:** AoS remains authoritative for SDAR operations, so no data loss

## Research Complete

The I/O phase is straightforward:
1. Keep existing AoS-based I/O code
2. Add SoA sync points after data is read into AoS
3. Output continues to work because AoS is maintained

---
*Generated: 2026-01-17*
