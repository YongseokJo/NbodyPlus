# Phase 25: Build System Integration - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

McLuster binary (mcluster_sse) compiles as part of ABYSS build system. Users run `make` and get both ABYSS and mcluster binaries. Build system handles missing Fortran compiler gracefully.

</domain>

<decisions>
## Implementation Decisions

### Binary location
- Symlink mcluster binary to same directory as ABYSS binary (build/ or bin/)
- Rename from `mcluster_sse` to `mcluster` for simpler invocation
- Original binary stays in mcluster/ directory, symlink points there

### Compiler detection
- Missing gfortran triggers warning and skips mcluster build (not hard error)
- Simple warning message: "gfortran not found, skipping mcluster" (no install hints)
- Support `DISABLE_MCLUSTER=1` flag to skip mcluster even when gfortran available

### Build targets
- Default `make` builds ABYSS + mcluster (if gfortran available)
- Explicit targets: `make mcluster`, `make mcluster-clean`, `make mcluster-rebuild`
- `make clean` removes both ABYSS and mcluster artifacts
- No separate `abyss-only` target — use `DISABLE_MCLUSTER=1` instead

### Build output
- Full compiler output by default (show all gfortran/gcc output)
- Support `QUIET=1` to suppress output for cleaner CI logs
- If mcluster build fails, stop entire build immediately
- Print clear "McLuster built successfully" message on success

### Claude's Discretion
- Exact Makefile implementation patterns
- How to detect current ABYSS binary location for symlink placement
- Specific compiler flags inherited from existing mcluster Makefile

</decisions>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 25-build-system-integration*
*Context gathered: 2026-01-20*
