# Phase 18: Particle Type Breakdown - Research

## Overview

Phase 18 aims to separate timing for different particle types (CM vs regular) and few-body operations to understand which cause more computational work.

## Existing Infrastructure

### Relevant TimerIDs in profiler.h
The profiler already has several relevant timer categories:
- `IrregularCMLoop` (line 33) - Already instrumented in compute_acceleration.cpp:160-208
- `FewBodyTermination` (line 47) - Already instrumented in irregular_routines.cpp:324-544
- `FewBodySearch` (line 48) - Uses old performance system, not PROFILE_START/STOP
- `FewBodyInitialization` (line 49) - Already instrumented in irregular_routines.cpp:608-751
- `FewBodyIntegration` (line 50) - NOT instrumented yet

### Existing Instrumentation Gaps
1. **FewBodySearch** - Uses legacy `performance.FewBodySearch` in irregular_routines.cpp
2. **FewBodyIntegration** - No PROFILE_START/STOP calls, needs instrumentation
3. **CM vs Regular particle separation** - No counter for particle types

## CM Particle Processing

### Detection Points
From `compute_acceleration.cpp`:
- Line 357: `if (this->is_cm_particle)` - identifies CM particles
- Line 418: `if (!ptcl->is_cm_particle)` - checks neighbor particles
- Line 639: `if (ptcl->is_cm_particle)` - CM particle neighbor handling

### IrregularCMLoop Timer
Currently times lines 160-208 in compute_acceleration.cpp which handles:
- CM particle neighbor iteration
- Computing forces from CM particle members

### Regular vs CM Distinction
The `compute_acceleration_irr()` function processes both regular and CM particles:
- Regular particles: standard neighbor loop
- CM particles: additional loop to iterate over member particles

## Few-Body Operations

### FewBodySearch (irregular_routines.cpp:551-604)
Binary search algorithm for new group formation:
- Uses `formBinaries()` function
- Currently timed with legacy performance struct

### FewBodyInitialization (irregular_routines.cpp:608-751)
Group initialization:
- Creates SDAR integrators
- Sets up group data structures
- Already instrumented with PROFILE_START/STOP

### FewBodyTermination (irregular_routines.cpp:324-544)
Group termination handling:
- Processes group dissolution
- Updates particle states
- Already instrumented with PROFILE_START/STOP

### FewBodyIntegration (NOT INSTRUMENTED)
SDAR integration in `fb_integration.cpp`:
- `Group::ARIntegration()` function
- Called via TASK_AR_INTEGRATION in worker loop
- Needs PROFILE_START/STOP instrumentation

## Requirements Analysis

### PTYPE-01: Separate timing for CM particles vs regular particles
**Current state:** IrregularCMLoop timer exists but doesn't separate counts
**Gap:** Need to track:
- Count of CM particles processed per interval
- Count of regular particles processed per interval
- Separate timing for CM vs regular force computation

### PTYPE-02: Track CM particle frequency and compute time ratio
**Gap:** Need new counters:
- `interval_cm_particle_count_`
- `interval_regular_particle_count_`
- Ratio calculation: `cm_time / regular_time`

### PTYPE-03: Few-body search time
**Current state:** Uses legacy `performance.FewBodySearch`
**Gap:** Need to migrate to PROFILE_START/STOP(TimerID::FewBodySearch)

### PTYPE-04: Few-body initialization time
**Current state:** Already instrumented with PROFILE_START/STOP
**Gap:** Need to ensure stats are output in CSV/JSON

### PTYPE-05: Few-body integration time
**Current state:** TimerID exists but NOT instrumented
**Gap:** Need PROFILE_START/STOP in fb_integration.cpp

## Implementation Strategy

### Plan 18-01: Particle Type Counter Infrastructure
- Add `interval_cm_particle_count_` and `interval_regular_particle_count_` to Profiler
- Add `PROFILE_PARTICLE_TYPE(is_cm)` macro
- Add `ParticleTypeBreakdown` struct for aggregated stats

### Plan 18-02: Instrument Particle Types
- Call PROFILE_PARTICLE_TYPE in compute_acceleration_irr()
- Ensure timing is captured per particle type

### Plan 18-03: Migrate FewBodySearch to Profiler
- Replace `performance.FewBodySearch +=` with PROFILE_START/STOP
- Update irregular_routines.cpp

### Plan 18-04: Instrument FewBodyIntegration
- Add PROFILE_START/STOP in ARIntegration()
- Location: fb_integration.cpp

### Plan 18-05: Particle Type Statistics Output
- Add CM/regular breakdown to console output
- Add CSV columns for particle type stats
- Add JSON section for particle type breakdown
- Compute and output ratios

## Key Files to Modify

1. `src/profiler.h`
   - Add particle type counters
   - Add PROFILE_PARTICLE_TYPE macro
   - Add output methods for new stats

2. `src/Particle/compute_acceleration.cpp`
   - Add PROFILE_PARTICLE_TYPE calls

3. `src/irregular_routines.cpp`
   - Migrate FewBodySearch to PROFILE_START/STOP

4. `src/FewBody/fb_integration.cpp`
   - Add PROFILE_START/STOP for FewBodyIntegration

## Dependencies

- Phase 15 complete (OnlineStats infrastructure)
- Phase 17 complete (per-particle timing)
- Existing few-body timers partially instrumented
