# Summary: 12-05 Overlap Measurement and Profiling

## Commit
ce86a38

## Deliverables
- TimerID::AsyncWindow - measures total async window duration
- TimerID::OverlapWork - measures local work (CM iteration) during window
- Profiler instrumentation around async window
- Overlap effectiveness percentage in performance summary
- Overlap timers in aggregated summary output

## Files Modified
- `src/profiler.h` (lines 73-75, 399-401, 474-488, 703-707)
- `src/irregular_routines.cpp` (lines 163-164, 173, 205, 218-219)

## Requirements Satisfied
- OVLP-02: MPI_Testany used for opportunistic completion checks
- OVLP-03: Profiler measures overlap effectiveness

## Profiler Metrics Added
1. **AsyncWindow**: Total time in async window (sends to completion)
2. **OverlapWork**: Time spent on CM particle iteration during window
3. **MPITestany**: Time in non-blocking completion tests

## Output Example
```
AsyncWindow                    0.42        5.2%       1234       340.51       892.15
OverlapWork                    0.08        1.0%       567        141.24       356.78
MPITestany                     0.15        1.9%      5678         26.42        89.34
Async Overlap Effectiveness: 19.0% of async window spent on local work
```

## Technical Notes
- Overlap effectiveness = (OverlapWork time) / (AsyncWindow time) * 100%
- Higher percentage means more useful work during communication latency
- Low percentage expected in simple test cases (few CM particles)
- Real benefit visible in production runs with many binaries
