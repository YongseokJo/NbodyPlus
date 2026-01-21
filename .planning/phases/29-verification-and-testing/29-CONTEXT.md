# Phase 29: Verification and Testing - Context

**Gathered:** 2026-01-21
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify the McLuster integration works end-to-end through comprehensive tests covering all three operational modes (generate+run, generate-only, run-only) plus energy conservation validation. Tests verify IC generation pipeline, not full N-body simulations.

</domain>

<decisions>
## Implementation Decisions

### Test execution approach
- Tests live alongside existing ABYSS tests
- Shell scripts contain test logic; make targets invoke them
- Two test modes: fixture-based (quick, no gfortran needed) and live McLuster (full end-to-end)
- Tests verify IC generation only, not full simulation runs

### Energy tolerance criteria
- 1e-4 relative tolerance for energy conservation
- Check initial energy only (not drift over timesteps)
- Hard fail if energy check fails

### Test configuration
- N=1000+ particles for realistic cluster sizes
- Multiple density profiles tested (Plummer, King, etc.)
- Cover key McLuster parameters: IMF variations, metallicity options
- TOML config files checked in as fixtures in test directory

### Failure reporting
- Per-test status output (test name + PASS/FAIL as each runs)
- On failure: show error message, expected vs actual, relevant file paths
- Exit codes: 0 = all pass, 1 = any failure
- Default: run all tests, report all failures; --fail-fast flag stops on first failure

### Claude's Discretion
- How to measure/read energy from ABYSS outputs
- Exact test file organization
- Which specific IMF/metallicity combinations to test
- Fixture file naming conventions

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches for test structure and execution.

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 29-verification-and-testing*
*Context gathered: 2026-01-21*
