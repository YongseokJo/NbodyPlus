# Summary: Plan 11-01 — Add Async Profiler Timers

**Status:** Complete
**Commit:** 6c7489d

## What Was Built

Added four new timer IDs to the profiler for measuring async MPI operations:
- `MPIIsend` — Time in MPI_Isend calls
- `MPIIrecv` — Time in MPI_Irecv calls
- `MPIWaitany` — Time in MPI_Waitany calls
- `MPITestany` — Time in MPI_Testany calls

## Requirements Satisfied

- [x] PROF-01: Timer for MPI_Isend operations
- [x] PROF-02: Timer for MPI_Irecv operations
- [x] PROF-03: Timer for MPI_Waitany operations
- [x] PROF-04: Timer for MPI_Testany operations

## Files Modified

| File | Change |
|------|--------|
| src/profiler.h | Added timer IDs to enum (lines 67-71) and name mappings (lines 390-394) |

## Verification

Grep confirms all four timers present in enum and name array.
