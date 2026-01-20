# Phase 26: Config Parser Extension - Context

**Gathered:** 2026-01-20
**Status:** Ready for planning

<domain>
## Phase Boundary

ABYSS TOML parser recognizes and validates `[mcluster]` configuration section with parameters: N, M, P, R, f, Z, b, e, generate_only. Parser provides clear error messages and enforces parameter constraints.

</domain>

<decisions>
## Implementation Decisions

### Parameter defaults
- All optional parameters have sensible defaults (except N/M — see below)
- Default N=10000 stars (if user specifies N without value, but N or M required)
- Default density profile: Plummer (P=0)
- Remaining params (R, f, Z, b, e, generate_only): match McLuster's own defaults

### Validation behavior
- Two-phase validation: basic checks at parse time, detailed checks at runtime
- Unknown parameters in [mcluster] section: error and stop (catches typos)
- Type mismatches: error immediately (strict types, no conversion attempts)
- Range validation: Claude decides which ranges are worth checking early

### N vs M handling
- N (star count) and M (total mass) are mutually exclusive
- If both specified: M takes precedence, warn user that N was ignored
- If neither specified: error — user must provide at least one
- M units: solar masses (M_sun)

### Error messages
- Verbose with hints: show valid ranges, suggest fixes
- Include line number from TOML file for easy location
- Stop at first error (not collecting all errors)
- Suggest similar parameter names for typos ("Unknown 'nstar', did you mean 'N'?")

### Claude's Discretion
- Exact range validation checks to implement at parse time
- How to detect similar parameter names (edit distance, etc.)
- Internal data structure for storing mcluster config

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

*Phase: 26-config-parser-extension*
*Context gathered: 2026-01-20*
