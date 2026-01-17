# Phase 5: SDAR Compatibility - Context

**Gathered:** 2026-01-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Maintain SDAR library integration via proxy pattern. SDAR expects AoS particle structures — we need to bridge between our new SoA layout and what SDAR requires. This phase creates the proxy/adapter layer for FewBody operations without modifying SDAR's core expectations.

</domain>

<decisions>
## Implementation Decisions

### Proxy Scope
- Claude determines which SDAR entry points need wrapping based on code analysis
- Support both individual particle and particle group access patterns
- Claude determines minimal vs full interface based on actual SDAR usage
- Surgical changes only to FewBody code — small targeted changes okay, no major restructuring

### Sync Timing
- Claude determines optimal sync strategy (eager vs lazy)
- Sync strategy must minimize computational overhead while maintaining physical integrity
- Mixed visibility: auto-sync by default, explicit override available for control
- Claude determines dirty tracking approach based on overhead analysis

### Data Ownership
- Buffer strategy optimized for minimal overhead (likely persistent or pooled)
- Working sets are small (2-10 particles) — binary/triple interactions
- Claude determines memory location (stack/heap/thread-local) based on call patterns
- SDAR calls are sequential — no thread-safety concerns for buffers

### Error Handling
- SDAR integration errors propagate up to caller — no retry at proxy layer
- Sync failures are fatal — abort simulation on inconsistent state
- Physical validation (energy/momentum checks) behind debug flag only
- Sync logging is configurable via debug level or runtime flag

### Claude's Discretion
- Exact set of SDAR functions to wrap
- Proxy interface design (full Particle mimic vs minimal)
- Sync strategy details (eager/lazy, dirty tracking)
- Buffer allocation approach (stack/heap/pool)

</decisions>

<specifics>
## Specific Ideas

- Small working sets (2-10 particles) means compact AoS buffers
- Sequential SDAR calls simplifies buffer management — no concurrency concerns
- Goal is minimal overhead — sync costs should not dominate FewBody compute time

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-sdar-compatibility*
*Context gathered: 2026-01-17*
