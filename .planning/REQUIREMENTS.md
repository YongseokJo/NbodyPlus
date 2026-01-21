# Requirements: ABYSS v3.0 McCluster Integration

**Defined:** 2026-01-20
**Core Value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations

## v3.0 Requirements

Requirements for McLuster IC generator integration. Each maps to roadmap phases.

### Build System

- [x] **BUILD-01**: McLuster source compiled with gfortran + gcc (SSE/BSE enabled)
- [x] **BUILD-02**: `mcluster_sse` binary produced alongside ABYSS binary
- [x] **BUILD-03**: Build system detects Fortran compiler availability
- [x] **BUILD-04**: Clean/rebuild targets work for mcluster

### Configuration

- [x] **CONFIG-01**: `[mcluster]` section parsed from TOML config
- [x] **CONFIG-02**: N (number of stars) parameter supported
- [x] **CONFIG-03**: M (total mass) parameter supported (alternative to N)
- [x] **CONFIG-04**: P (density profile) parameter: Plummer, King, etc.
- [x] **CONFIG-05**: R (half-mass radius in pc) parameter
- [x] **CONFIG-06**: f (IMF selection) parameter
- [x] **CONFIG-07**: Z (metallicity) parameter
- [x] **CONFIG-08**: b (binary fraction) parameter
- [x] **CONFIG-09**: e (stellar evolution epoch in Myr) parameter
- [x] **CONFIG-10**: `generate_only` flag for IC-only mode

### Runtime Integration

- [ ] **RUNTIME-01**: ABYSS main.cpp detects `[mcluster]` section presence
- [ ] **RUNTIME-02**: McLuster subprocess spawned with correct arguments
- [ ] **RUNTIME-03**: Wait for McLuster completion before simulation
- [ ] **RUNTIME-04**: McLuster output captured and validated
- [ ] **RUNTIME-05**: Exit after IC generation if `generate_only = true`
- [ ] **RUNTIME-06**: Proceed to simulation if `generate_only = false` (default)
- [ ] **RUNTIME-07**: Existing IC file mode preserved (no `[mcluster]` section)

### Output Compatibility

- [ ] **OUTPUT-01**: McLuster output format matches ABYSS nbody.dat expectations
- [ ] **OUTPUT-02**: Units conversion if needed (N-body vs astrophysical)
- [ ] **OUTPUT-03**: Generated IC file placed in correct location for ABYSS

### Verification

- [ ] **VERIFY-01**: End-to-end test: config -> McLuster -> ABYSS simulation
- [ ] **VERIFY-02**: Generate-only test: config -> McLuster -> IC file (no simulation)
- [ ] **VERIFY-03**: Run-only test: existing IC -> ABYSS simulation (regression)
- [ ] **VERIFY-04**: Energy conservation check with McLuster-generated ICs

## v3.1+ Requirements

Deferred to future release.

### Extended Features

- **EXT-01**: Additional McLuster parameters (fractal dimension, mass segregation)
- **EXT-02**: Multiple IC generation in batch mode
- **EXT-03**: IC visualization/validation tools
- **EXT-04**: Parameter sweep support

## Out of Scope

| Feature | Reason |
|---------|--------|
| McLuster GUI | Command-line workflow sufficient |
| McLuster source modifications | Use upstream code as-is |
| Non-SSE McLuster build | SSE/BSE stellar evolution required |
| MPI optimization (v2.4 work) | Separate milestone, deferred |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| BUILD-01 | Phase 25 | Complete |
| BUILD-02 | Phase 25 | Complete |
| BUILD-03 | Phase 25 | Complete |
| BUILD-04 | Phase 25 | Complete |
| CONFIG-01 | Phase 26 | Complete |
| CONFIG-02 | Phase 26 | Complete |
| CONFIG-03 | Phase 26 | Complete |
| CONFIG-04 | Phase 26 | Complete |
| CONFIG-05 | Phase 26 | Complete |
| CONFIG-06 | Phase 26 | Complete |
| CONFIG-07 | Phase 26 | Complete |
| CONFIG-08 | Phase 26 | Complete |
| CONFIG-09 | Phase 26 | Complete |
| CONFIG-10 | Phase 26 | Complete |
| RUNTIME-01 | Phase 27 | Pending |
| RUNTIME-02 | Phase 27 | Pending |
| RUNTIME-03 | Phase 27 | Pending |
| RUNTIME-04 | Phase 27 | Pending |
| RUNTIME-05 | Phase 27 | Pending |
| RUNTIME-06 | Phase 27 | Pending |
| RUNTIME-07 | Phase 27 | Pending |
| OUTPUT-01 | Phase 28 | Pending |
| OUTPUT-02 | Phase 28 | Pending |
| OUTPUT-03 | Phase 28 | Pending |
| VERIFY-01 | Phase 29 | Pending |
| VERIFY-02 | Phase 29 | Pending |
| VERIFY-03 | Phase 29 | Pending |
| VERIFY-04 | Phase 29 | Pending |

**Coverage:**
- v3.0 requirements: 28 total
- Mapped to phases: 28
- Unmapped: 0

---
*Requirements defined: 2026-01-20*
*Last updated: 2026-01-20 - Phase 26 complete (14/28 requirements done)*
