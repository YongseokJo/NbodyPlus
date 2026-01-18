# Plan 07-01: Run Final Validation Test

## Frontmatter
```yaml
phase: 7
plan: 1
wave: 1
depends_on: []
files_modified: []
autonomous: true
```

## Objective
Run the final SoA validation test using the existing workflow infrastructure and capture metrics in summary_runs.tsv.

## Context
- Baseline metrics captured: dE/E0 = 3.20665e-05, wall time = 28.33s
- All intermediate phase runs (phase4, phase5, phase6) show acceptable energy conservation
- Workflow already handles compile → run → analyze pipeline

## Tasks

<task id="1">
<action>run</action>
<description>Submit final validation run with tag `soa_final`</description>
<command>cd /gpfs/home/vjl4366/pkg/ABYSS && workflow/bin/submit.sh --tag soa_final --scheduler slurm</command>
<notes>This will compile current code, run simulation on tests/test1, and analyze results. Run completes asynchronously via SLURM.</notes>
</task>

<task id="2">
<action>verify</action>
<description>Check job submission</description>
<command>squeue -u $USER</command>
<expected>Job appears in queue with name abyss_soa_final</expected>
</task>

<task id="3">
<action>wait</action>
<description>Wait for job completion</description>
<notes>Job typically completes in ~5-10 minutes including compile time. Check periodically with `squeue -u $USER` or wait for SLURM email notification.</notes>
</task>

## Verification Criteria

- [ ] Job submitted successfully (submit.log created)
- [ ] Run directory created under workflow/runs/soa_final_*
- [ ] Job completes without error

## must_haves
- Final validation run submitted with tag `soa_final`
- Run uses same test configuration as baseline (tests/test1, 16 ranks, 1 GPU)
