# Phase 21: Instrumentation Fixes — Verification Report

**Phase:** 21
**Name:** Instrumentation Fixes
**Date:** 2026-01-19
**Status:** partial_success (architectural limitation)

## Goal Verification

**Phase Goal:** Fix missing profiling data identified in v2.2 analysis

### Requirements Status

| Requirement | Description | Status | Notes |
|-------------|-------------|--------|-------|
| FIX-01 | Fix neighbor count instrumentation | ⚠ Blocked | MPI aggregation design incompatible |
| FIX-02 | Fix worker compute time tracking | ⚠ Blocked | MPI aggregation design incompatible |
| FIX-03 | Fix CM particle type breakdown | ⚠ Blocked | MPI aggregation design incompatible |
| FIX-04 | Cache statistics fallback | ✓ Complete | Status messages working |

## Architectural Issue Discovered

The MPI aggregation approach using `MPI_Reduce` is **fundamentally incompatible** with the root/worker architecture:

### The Problem

```
ROOT (rank 0)                    WORKERS (ranks 1-N)
     |                                |
root_routines.cpp               worker_routines.cpp
     |                                |
     | aggregateAcrossRanks()         | (different code path)
     |                                |
     v                                |
  MPI_Reduce() ─────────────────> NEVER CALLED
     |                                |
  DEADLOCK!                           |
```

- `MPI_Reduce` is a **collective operation** requiring ALL ranks to participate
- Root and workers execute different code paths with no shared synchronization point
- Enabling `-DUSE_MPI` causes deadlock on first output

### What Works (Root-Side Data)

These metrics are collected on root and work correctly:
- Queue dispatch statistics (depth, starvation events)
- Worker assignment (particles per worker)
- Timer statistics (IrregularForce, RegularGPU, etc.)

### What Doesn't Work (Worker-Side Data)

These require aggregation from workers to root:
- Neighbor count profiling
- Worker compute times
- Particle type breakdown (CM vs regular)

## Recommended Solution

**Phase 21.5: Redesign Aggregation** should use point-to-point MPI:

```cpp
// Instead of MPI_Reduce (collective):
// Have workers send data to root at a sync point

// Option 1: Piggyback on existing worker→root messages
// Workers already send force results; add profiling data

// Option 2: Dedicated profiling message at output time
// Root broadcasts "send profiling data" signal
// Workers respond with MPI_Send, root collects with MPI_Recv
```

## Current Code State

- `USE_MPI` is **NOT defined** (intentional - prevents deadlock)
- Aggregation code exists but is inactive
- Root-side profiling works correctly
- Simulation runs without issues

## Verification Summary

- **Code Implementation:** ✓ Complete (aggregation code exists)
- **Runtime Behavior:** ⚠ Partial (root data only, worker data not aggregated)
- **Simulation Stability:** ✓ Runs correctly

## Next Steps

1. **Accept partial solution:** Use root-side profiling data (queue stats, worker assignment)
2. **Future work (Phase 21.5):** Redesign aggregation using point-to-point MPI

---
*Verification completed: 2026-01-19*
*Architectural limitation identified - requires redesign for full functionality*
