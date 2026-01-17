# Phase 4 Plan Verification

## Phase Goal
Update all CPU force calculation and integration routines to use SoA access patterns.

## Requirements Coverage

| Requirement | Plan(s) | Covered |
|-------------|---------|---------|
| CPU-01: Regular force SoA | 04-PLAN | ✓ |
| CPU-02: Irregular force SoA | 04-PLAN | ✓ |
| CPU-03: Prediction routines | 03-PLAN | ✓ |
| CPU-04: Correction routines | 03-PLAN | ✓ |
| CPU-05: Timestep routines | 02-PLAN | ✓ |

## Wave Structure

| Wave | Plans | Purpose |
|------|-------|---------|
| 1 | 01, 02 | Foundation (helpers, timestep) |
| 2 | 03 | Prediction/correction |
| 3 | 04 | Force calculations (core) |
| 4 | 05 | Integration/orchestration |

## Dependency Analysis

```
01-PLAN (SoA helpers)
   ↓
03-PLAN (prediction/correction) ←── depends on helpers
   ↓
04-PLAN (force calculations) ←── depends on helpers + prediction
   ↓
05-PLAN (orchestration) ←── depends on force calculations

02-PLAN (timestep) ─── independent, can run parallel with 01
   ↓
05-PLAN (orchestration)
```

## Risk Assessment

| Plan | Risk | Mitigation |
|------|------|------------|
| 01 | Low | Pure additions, no existing code changed |
| 02 | Low | Wrappers preserve backward compatibility |
| 03 | Medium | Free functions must match method behavior exactly |
| 04 | High | Core physics — verify bit-identical results |
| 05 | Medium | Integration testing required |

## Energy Conservation Verification

Critical metric: dE/E0 ≤ 3.2e-5 (baseline: 3.20665e-05)

Verification points:
1. After Plan 03: Run prediction/correction unit test
2. After Plan 04: Compare acceleration values to baseline
3. After Plan 05: Full simulation with energy check

## Gaps Identified

1. **FewBody/SDAR fields** — Not addressed in Phase 4 (deferred to Phase 5)
2. **MPI sync** — Not addressed (deferred to Phase 2)
3. **Full ParticleData ownership** — Code still primarily uses Particle struct

## Conclusion

Plans cover all requirements with appropriate wave ordering. The bridge approach (keeping Particle methods, adding SoA free functions) minimizes risk while enabling SoA benefits in hot paths.

---
*Verified: 2026-01-17*
