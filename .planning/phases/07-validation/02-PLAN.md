# Plan 07-02: Analyze Results and Compare to Baseline

## Frontmatter
```yaml
phase: 7
plan: 2
wave: 2
depends_on: [1]
files_modified: []
autonomous: true
```

## Objective
Verify energy conservation and performance metrics from the soa_final run against baseline.

## Context
**Baseline (from summary_runs.tsv):**
| Metric | Value |
|--------|-------|
| dE/E0 mean | 3.20665e-05 |
| Wall time | 28.33s |
| Git commit | 7e90c74 |

**Acceptance criteria:**
- VAL-02: dE_over_E0_mean ≤ 3.2e-5 (or within acceptable margin)
- VAL-03: status=ok (Simulation Done found)
- VAL-04: Performance improvement measurable

## Tasks

<task id="1">
<action>verify</action>
<description>Check run completed successfully</description>
<command>grep "status=" /gpfs/home/vjl4366/pkg/ABYSS/workflow/runs/soa_final_*/summary.txt</command>
<expected>status=ok (Simulation Done found)</expected>
</task>

<task id="2">
<action>verify</action>
<description>Check energy conservation</description>
<command>grep "soa_final" /gpfs/home/vjl4366/pkg/ABYSS/summary_runs.tsv</command>
<expected>dE_over_E0_mean column shows value ≤ 3.2e-5</expected>
</task>

<task id="3">
<action>compare</action>
<description>Compare soa_final vs baseline metrics</description>
<command>head -2 /gpfs/home/vjl4366/pkg/ABYSS/summary_runs.tsv && grep "soa_final" /gpfs/home/vjl4366/pkg/ABYSS/summary_runs.tsv</command>
<notes>
Columns to compare:
- dE_over_E0_mean: should be ≤ 3.20665e-05
- total_wall_s: compare against 28.33s baseline
</notes>
</task>

<task id="4">
<action>verify</action>
<description>Review energy analysis output</description>
<command>cat /gpfs/home/vjl4366/pkg/ABYSS/workflow/runs/soa_final_*/analyze_energy.out</command>
<expected>Summary shows Max and Mean |dE/E0| values</expected>
</task>

## Verification Criteria

- [ ] Simulation completed with status=ok
- [ ] dE_over_E0_mean ≤ 3.2e-5 (baseline level)
- [ ] total_wall_s captured (for performance comparison)
- [ ] No errors in run.log

## must_haves
- VAL-02 satisfied: energy conservation matches or beats baseline
- VAL-03 satisfied: simulation completes successfully
- VAL-04 data captured: wall time available for comparison
