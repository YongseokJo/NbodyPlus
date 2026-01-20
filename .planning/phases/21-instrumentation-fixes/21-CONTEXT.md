# Phase 21: Instrumentation Fixes - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix the four missing profiling metrics identified in v2.2 analysis:
1. Neighbor count (showing 0 despite 140K IrregularForce calls)
2. Worker compute time (showing 0.000s despite ~120s total)
3. CM particle tracking (showing 0 — may be correct if no binaries)
4. Cache statistics fallback (perf_event unavailable on HPC)

This phase fixes broken instrumentation. Profiler reorganization is deferred to Phase 21.5.

</domain>

<decisions>
## Implementation Decisions

### Diagnosis Approach
- Claude decides most efficient approach (read code paths, add debug logging, or trace execution)
- Claude decides fix order based on logical dependencies

### Fix Approach
- Fix properly, not minimal patches — quality matters even if Phase 21.5 rewrites
- Rewrite is acceptable if cleaner than patching
- Each fix should be clean and maintainable

### Cache Fallback Strategy
- Claude decides most practical fallback when perf_event unavailable
- Options include: skip silently, estimate from timing, or sample-based estimation

### Verification Method
- Re-run full profiling after fixes
- Use: `workflow/bin/submit.sh --test-dir tests/test4 --tag <descriptive_tag>`
- Confirm non-zero values for all four metrics

### Claude's Discretion
- Diagnosis approach and order
- Cache fallback implementation
- Specific fix techniques per issue

</decisions>

<specifics>
## Specific Ideas

From profiling analysis `profiling_20260119_013613`:
```
neighbor_profiling.count = 0  (expected: ~281,463)
worker_compute_time = 0.000s  (expected: ~120s distributed across 15 workers)
cm_count = 0                  (verify if test has binaries)
counters_available = false    (need fallback)
```

Key files to investigate:
- `src/Particle/compute_acceleration.cpp` — Neighbor count macro placement
- `src/queue_scheduler.h` — Worker compute time recording
- `src/irregular_routines.cpp` — CM particle detection
- `src/profiler.h` — Cache stats fallback

</specifics>

<deferred>
## Deferred Ideas

- **Profiler reorganization** — Phase 21.5: Combine and streamline redundant profiling code
- **New profiling metrics** — Future phases as needed

</deferred>

---

*Phase: 21-instrumentation-fixes*
*Context gathered: 2026-01-19*
