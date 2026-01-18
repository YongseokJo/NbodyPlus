# Phase 10: Targeted Optimization - Context

**Gathered:** 2026-01-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Implement SIMD vectorization for the irregular force bottleneck (53.5% of wall time identified in Phase 9). Validate with before/after profiling using test1 case. Verify energy conservation within baseline tolerance.

**In scope:** AVX-512 vectorization, data layout restructuring, neighbor list optimization
**Deferred to v2.1:** MPI batching, GPU irregular forces

</domain>

<decisions>
## Implementation Decisions

### Optimization Strategy
- Target irregular force loop with explicit AVX-512 intrinsics
- Use compile-time preprocessor switch (`#ifdef __AVX512F__`) for AVX level selection
- Focus on AVX-512 as baseline (cluster has support)
- MPI batching deferred to v2.1 — v2.0 is compute optimization only

### SIMD Approach
- Explicit intrinsics, not compiler auto-vectorization
- AVX-512 baseline with preprocessor-controlled paths
- Scalar fallback for non-AVX systems via `#else` branch

### Data Layout
- Restructure particle data layout if needed for better SIMD utilization
- Optimize neighbor list structure for cache-friendly access patterns
- Aggressive restructuring permitted — major changes acceptable

### Validation Approach
- Use same test1 case as Phase 9 baseline for comparison
- Target: at least 20% improvement in IrregularForce time
- Energy conservation tolerance: dE/E0 < 1e-4 (same as baseline)
- Accept small floating-point differences from SIMD operation reordering

### Development Testing
- Unit tests for vectorized functions in isolation
- Energy conservation check after each significant code change
- Run full test1 simulation to validate correctness during development

### Risk Management
- Aggressive restructuring is acceptable for maximum SIMD gains
- Fallback plan: git revert and try simpler approach if issues arise
- No compile-time flag to toggle optimization — commit to the restructuring

### Claude's Discretion
- Specific AVX-512 intrinsic choices
- Loop unrolling factors
- Memory alignment strategy
- Exact neighbor list data structure design

</decisions>

<specifics>
## Specific Ideas

- Phase 9 identified IrregularForce at 53.5% of wall time as primary bottleneck
- ~105 million MPI messages per interval (QueueRun/MPISend overhead) — acknowledged but deferred
- Irregular forces evaluated ~122x more frequently than regular forces
- Compile with `-mavx512f` flag for AVX-512 support

</specifics>

<deferred>
## Deferred Ideas

- **MPI batching** — Reduce 105M messages/interval by grouping particles. Significant impact but separate effort. → v2.1
- **GPU irregular forces** — Port irregular force to GPU kernel. Major undertaking. → v2.1
- **Runtime AVX detection** — Compile-time preprocessor is sufficient for HPC clusters.

</deferred>

---

*Phase: 10-targeted-optimization*
*Context gathered: 2026-01-17*
