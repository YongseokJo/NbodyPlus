# Phase 22: Profiler Cleanup - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Remove unused profiling timers, consolidate redundant code, and streamline JSON output. This is infrastructure cleanup — no new profiling capabilities.

</domain>

<decisions>
## Implementation Decisions

### Timer Removal Criteria
- Timers that are always 0: **Keep but disable** (comment out or #ifdef for future use)
- Timers with instrumentation but never called (IrregularNeighborLoop, IrregularCMLoop): **Remove completely**
- Duplicate timers like RegularForce/RegularGPU: **Keep both, rename RegularForce → RegularCPU** for clarity
- Worker-side timers not aggregated (WorkerRecvWait, WorkerTaskDispatch, WorkerSendComplete): **Remove completely**

### JSON Output Structure
- Top-level summary section: **Claude's discretion** (determine what summary fields are most useful)
- Timer output: **Non-zero plus explicitly enabled** (output non-zero timers, plus any marked as "always output")
- Timer verbosity: **Two modes** — verbose flag for full output (all 7 fields), default to simplified
- Specialized sections (neighbor_profiling, worker_distribution, etc.): **Claude's discretion** to reorganize
- Add **schema_version field** (e.g., "schema_version": "2.3") at top of JSON

### Code Organization
- Split profiler.h: **Minimal split** — profiler.h (main) + profiler_types.h (structs/enums only)
- Singleton pattern: **Claude's discretion** (evaluate trade-offs for this codebase)
- Helper classes (OnlineStats, Histogram, WorkerDistributionTracker): **Claude's discretion** based on dependencies
- PROFILE_* macros: **Consolidate to essentials** — PROFILE_START/STOP, PROFILE_SCOPE, PROFILE_NEIGHBOR_TIME, remove the rest

### Backward Compatibility
- TimerID removal: **Remove and update immediately** — clean break, remove enum entries and all call sites
- Analysis tools (analyze_profiling.py): **Claude's discretion** whether immediate update needed
- Old profiling output files: **No backward compat needed** — old runs can be re-run if needed

### Claude's Discretion
- Whether to add summary section and what fields to include
- Specialized section reorganization
- Singleton vs instance pattern
- Helper class placement (profiler_types.h vs inline)
- Analysis tool update timing

</decisions>

<specifics>
## Specific Ideas

- Rename RegularForce → RegularCPU to clarify it's for potential CPU-only path
- schema_version should be "2.3" matching milestone version
- Disabled timers should use #ifdef blocks rather than deletion, for easy re-enablement

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 22-profiler-cleanup*
*Context gathered: 2026-01-19*
