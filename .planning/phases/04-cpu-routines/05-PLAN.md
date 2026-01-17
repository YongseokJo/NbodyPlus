# Plan 04-05: Update Routine Orchestration and Integration

## Frontmatter

```yaml
wave: 4
depends_on: [02-PLAN.md, 04-PLAN.md]
files_modified:
  - src/regular_routines.cpp
  - src/irregular_routines.cpp
autonomous: true
```

## Objective

Update the main routine orchestration files to work with SoA-aware subroutines:
- `regular_routines.cpp` — Regular force dispatch
- `irregular_routines.cpp` — Irregular force loop with skip list

## Context

These files orchestrate the force calculation workflow:
1. Build particle lists (RegularList, IrregularList)
2. Dispatch to queue scheduler or direct calculation
3. Handle FewBody group formation/termination
4. Update particle states

Current pattern uses `particles[index]` to access Particle structs. With SoA subroutines from Plan 04-04, we need to ensure data flows correctly.

## Tasks

<task id="1">
Update regular_routines.cpp:
- Ensure RegularRoutines() passes particle data correctly to subroutines
- Update any direct particle field access to use ParticleData when beneficial
</task>

<task id="2">
Update irregular_routines.cpp:
- IrregularRoutines() main loop
- Skip list creation/update still uses particles[] but force calculation uses SoA
- Update particle state access patterns where hot
</task>

<task id="3">
Add ParticleData initialization in main simulation loop (if not already present):
- Ensure sync_from_particle() is called before force calculation
- Ensure sync_to_particle() is called after results computed
</task>

<task id="4">
Update GPU integration path (if CPU-GPU hybrid):
- Ensure data flows from ParticleData to GPU (using Phase 3 infrastructure)
</task>

## Verification

- [ ] Full simulation completes without errors
- [ ] Energy conservation within baseline tolerance (≤3.2e-5)
- [ ] Performance comparable or improved

## must_haves

- Simulation produces correct output
- Energy conservation dE/E0 ≤ 3.2e-5 (baseline: 3.20665e-05)
- No regression in wall time
- FewBody/SDAR integration still functional
