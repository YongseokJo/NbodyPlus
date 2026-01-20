# Phase 23: Deep Analysis - Context

**Gathered:** 2026-01-19
**Status:** Ready for planning

<domain>
## Phase Boundary

Comprehensive profiling analysis to understand bottleneck breakdown and prioritize optimizations. Run extended profiling, analyze where time goes, and produce prioritized recommendations for Phases 24-26. No code optimization in this phase — analysis and reporting only.

</domain>

<decisions>
## Implementation Decisions

### Analysis Depth
- Run 100 steps for the main profiling run
- Use mixture approach: 1 long run (100 steps) + 2 short runs (20 steps each) for variance estimation
- Profiling interval: every 10 steps (current default)

### Output Format
- Single comprehensive Markdown file (ANALYSIS.md)
- Location: `.planning/` top-level for visibility
- Include both high-level summary (3-5 buckets: Compute, Communication, Idle) AND detailed appendix with full timer breakdown
- Charts optional — text tables primary, generate PNG charts if easy

### Priority Matrix Criteria
- Rank optimizations by **impact only** (expected speedup percentage)
- Minimum threshold: **15%+ expected speedup** to pursue an optimization
- Present data for go/no-go decisions — user decides whether to proceed with Phases 24-25
- Calculate impact using Amdahl's Law (theoretical) + empirical comparison + confidence range

### Scaling Investigation
- Test two configurations: current test case + larger configuration
- User will provide two different 100K particle ICs for larger tests
- Primary focus: **wall time scaling** (how runtime grows with particle count)
- Run on SLURM cluster (submit jobs, realistic HPC environment)

### Claude's Discretion
- Exact structure of analysis report sections
- Which charts to generate (if any)
- How to present confidence ranges
- Appendix organization

</decisions>

<specifics>
## Specific Ideas

- User provides the 100K particle ICs — don't generate synthetic test cases
- Focus on actionable numbers: "X% of time in Y, reducing by Z gives W% speedup"
- Keep report readable — summary at top for quick decisions

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 23-deep-analysis*
*Context gathered: 2026-01-19*
