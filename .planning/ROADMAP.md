# Roadmap: v3.0 McCluster Integration

**Created:** 2026-01-20
**Milestone:** v3.0
**Starting Phase:** 25 (continuing from v2.3)
**Phases:** 5 (25-29)
**Requirements:** 28 v3.0 requirements

## Overview

Integrating McLuster IC generator into ABYSS enables seamless end-to-end simulation workflows. This milestone adds build system support for McLuster compilation, extends the TOML config parser with `[mcluster]` section, implements subprocess orchestration in main.cpp, handles output format compatibility, and verifies the integration with comprehensive tests.

## Phase Structure

### Phase 25: Build System Integration

**Goal:** McLuster binary (mcluster_sse) compiles as part of ABYSS build system

**Dependencies:** None (foundation phase)

**Requirements:**
- BUILD-01: McLuster source compiled with gfortran + gcc (SSE/BSE enabled)
- BUILD-02: `mcluster_sse` binary produced alongside ABYSS binary
- BUILD-03: Build system detects Fortran compiler availability
- BUILD-04: Clean/rebuild targets work for mcluster

**Plans:** 1 plan

Plans:
- [x] 25-01-PLAN.md - Create root Makefile with mcluster integration

**Success Criteria:**
1. User runs `make` and mcluster_sse binary appears in expected location
2. User without gfortran gets clear error message (not cryptic build failure)
3. User runs `make clean` and mcluster artifacts are removed
4. User can rebuild mcluster independently with `make mcluster`

---

### Phase 26: Config Parser Extension

**Goal:** ABYSS TOML parser recognizes and validates `[mcluster]` configuration section

**Dependencies:** Phase 25 (need binary to know valid parameter ranges)

**Requirements:**
- CONFIG-01: `[mcluster]` section parsed from TOML config
- CONFIG-02: N (number of stars) parameter supported
- CONFIG-03: M (total mass) parameter supported (alternative to N)
- CONFIG-04: P (density profile) parameter: Plummer, King, etc.
- CONFIG-05: R (half-mass radius in pc) parameter
- CONFIG-06: f (IMF selection) parameter
- CONFIG-07: Z (metallicity) parameter
- CONFIG-08: b (binary fraction) parameter
- CONFIG-09: e (stellar evolution epoch in Myr) parameter
- CONFIG-10: `generate_only` flag for IC-only mode

**Plans:** 2 plans

Plans:
- [x] 26-01-PLAN.md - Add keys() to toml.hpp and create MclusterConfig struct
- [x] 26-02-PLAN.md - Implement mcluster config parsing and validation

**Success Criteria:**
1. User adds `[mcluster]` section to config and ABYSS parses without error
2. User specifies N=1000 and value appears in parsed parameters
3. User specifies invalid parameter (e.g., negative N) and gets validation error
4. User can specify either N or M (mutually exclusive) for cluster size
5. User omits `[mcluster]` section and existing behavior is unchanged

---

### Phase 27: Runtime Integration

**Goal:** ABYSS main.cpp orchestrates McLuster subprocess based on config presence

**Dependencies:** Phase 26 (need parsed config to drive runtime behavior)

**Requirements:**
- RUNTIME-01: ABYSS main.cpp detects `[mcluster]` section presence
- RUNTIME-02: McLuster subprocess spawned with correct arguments
- RUNTIME-03: Wait for McLuster completion before simulation
- RUNTIME-04: McLuster output captured and validated
- RUNTIME-05: Exit after IC generation if `generate_only = true`
- RUNTIME-06: Proceed to simulation if `generate_only = false` (default)
- RUNTIME-07: Existing IC file mode preserved (no `[mcluster]` section)

**Plans:** 2 plans

Plans:
- [x] 27-01-PLAN.md - Create mcluster_runner module (subprocess, validation, transform)
- [x] 27-02-PLAN.md - Integrate mcluster runner into main.cpp

**Success Criteria:**
1. User runs ABYSS with `[mcluster]` config and sees McLuster subprocess output
2. User sees ABYSS wait for McLuster to complete before starting simulation
3. User with `generate_only = true` sees McLuster run then ABYSS exit cleanly
4. User without `[mcluster]` section sees existing IC-file behavior (no McLuster)
5. User sees clear error message if McLuster subprocess fails

---

### Phase 28: Output Format Handling

**Goal:** McLuster output is compatible with ABYSS nbody.dat format expectations

**Dependencies:** Phase 27 (need subprocess working to test output handling)

**Requirements:**
- OUTPUT-01: McLuster output format matches ABYSS nbody.dat expectations
- OUTPUT-02: Units conversion if needed (N-body vs astrophysical)
- OUTPUT-03: Generated IC file placed in correct location for ABYSS

**Plans:** 1 plan

Plans:
- [x] 28-01-PLAN.md - Fix unit conversion in transformMclusterOutput()

**Success Criteria:**
1. User runs integrated workflow and ABYSS reads McLuster-generated IC without errors
2. User sees IC file in expected location (configurable or standard path)
3. User with astrophysical units in McLuster config gets proper N-body unit conversion

---

### Phase 29: Verification and Testing

**Goal:** Integration verified with comprehensive end-to-end tests

**Dependencies:** Phase 28 (need complete integration to verify)

**Requirements:**
- VERIFY-01: End-to-end test: config -> McLuster -> ABYSS simulation
- VERIFY-02: Generate-only test: config -> McLuster -> IC file (no simulation)
- VERIFY-03: Run-only test: existing IC -> ABYSS simulation (regression)
- VERIFY-04: Energy conservation check with McLuster-generated ICs

**Success Criteria:**
1. User runs test suite and sees all three modes pass (generate+run, generate-only, run-only)
2. User sees energy conservation within tolerance for McLuster-generated cluster
3. User can verify regression test confirms existing IC file mode unchanged
4. User has example config demonstrating `[mcluster]` usage

---

## Progress

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 25 | Build System Integration | BUILD-01, BUILD-02, BUILD-03, BUILD-04 | Complete |
| 26 | Config Parser Extension | CONFIG-01 through CONFIG-10 | Complete |
| 27 | Runtime Integration | RUNTIME-01 through RUNTIME-07 | Complete |
| 28 | Output Format Handling | OUTPUT-01, OUTPUT-02, OUTPUT-03 | Complete |
| 29 | Verification and Testing | VERIFY-01, VERIFY-02, VERIFY-03, VERIFY-04 | Pending |

**Coverage:** 28/28 requirements mapped

---

## Dependency Graph

```
Phase 25 (Build)
    |
    v
Phase 26 (Config)
    |
    v
Phase 27 (Runtime)
    |
    v
Phase 28 (Output)
    |
    v
Phase 29 (Verification)  <-- CURRENT
```

Linear dependency chain - each phase unlocks the next.

---
*Roadmap created: 2026-01-20*
*Last updated: 2026-01-21 - Phase 28 complete (1 plan executed)*
