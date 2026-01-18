# Plan 03: Verify HDF5 Output Works with SoA

## Frontmatter

```yaml
phase: 6
plan: 03
wave: 2
depends_on: [01]
files_modified: []
autonomous: true
```

## Objective

Verify that HDF5 output (`writeParticle()`) continues to work correctly after the SoA integration. The output reads from AoS which is maintained via exit-only sync pattern, so no code changes should be needed - just verification.

## Context

- `writeParticle()` in `src/read_write.cpp` reads from `particles[]` AoS array
- The AoS is kept in sync with SoA via exit-only pattern (Phase 5)
- HDF5 datasets written: PID, Mass_Msun, X_pc, Y_pc, Z_pc, Vx_km_s, Vy_km_s, Vz_km_s
- Output uses predicted values via `predict_particle_second_order()`

## Tasks

<task id="1">
Review `writeParticle()` function to confirm it reads from AoS and doesn't need changes for SoA compatibility.
</task>

<task id="2">
Run a test simulation to generate HDF5 output:
```bash
./build.sh && cd tests/test1 && ./test.sh
```
</task>

<task id="3">
Verify the output HDF5 file contains valid particle data:
```bash
h5dump -H output/*.h5  # Check structure
```
</task>

<task id="4">
If the test produces output, verify the energy analysis succeeds:
```bash
python3 tools/analyze_energy.py tests/test1/output/*.h5
```
</task>

## Verification

- [ ] HDF5 output file is created
- [ ] Output contains expected datasets (PID, Mass, positions, velocities)
- [ ] Particle data in output is physically reasonable
- [ ] Energy analysis tool can read the output

## must_haves

- [ ] IO-01: "output.h5 contains correct particle data"
- [ ] No regression in output quality

---
*Generated: 2026-01-17*
