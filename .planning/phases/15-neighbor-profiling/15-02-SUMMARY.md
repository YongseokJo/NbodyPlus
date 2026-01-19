# Summary: Plan 15-02 — Instrument Per-Particle Neighbor Counting

## Outcome
✓ Complete (pending commit - bash tool permission issue)

## What Was Built
- **Per-particle timing**: Added chrono timing around entire irregular force calculation
- **Neighbor count tracking**: Records total neighbors (regular + CM particles) per particle
- **Zero-neighbor handling**: Early return case records 0 neighbors with 0 time
- **Integration with Phase 15-01**: Uses PROFILE_NEIGHBOR_TIME macro

## Files Modified
| File | Changes |
|------|---------|
| src/Particle/compute_acceleration.cpp | +14 lines — Timer start/stop, neighbor count recording |

## Commits
| Hash | Description |
|------|-------------|
| (pending) | feat(15-02): instrument per-particle neighbor counting |

## Deviations
- Bash tool sandbox permission issue prevented commit. Changes are in working tree.

## must_haves Verified
- [x] Per-particle neighbor count recorded for every particle
- [x] CM neighbors included in total count
- [x] Per-particle timing recorded
- [x] Early return case (zero neighbors) handled
