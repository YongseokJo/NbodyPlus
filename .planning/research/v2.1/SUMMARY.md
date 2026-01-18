# Research Summary: MPI Communication Optimization

## Key Findings

### Stack Recommendation

**Non-blocking MPI with request management** is the recommended approach:

| Function | Purpose |
|----------|---------|
| `MPI_Isend` | Non-blocking task send |
| `MPI_Irecv` | Non-blocking result receive |
| `MPI_Waitany` | Wait for first completion |
| `MPI_Testany` | Non-blocking completion check |

Core pattern:
1. Post all receives first (one per worker)
2. Send tasks asynchronously (don't wait after each)
3. Process completions as they arrive (`MPI_Waitany`)
4. Re-post receives, send next tasks

### Features Priority

| Priority | Feature |
|----------|---------|
| P0 | Non-blocking send/recv |
| P0 | Request lifecycle management |
| P0 | Buffer safety |
| P1 | Communication-computation overlap |
| P1 | Profiling integration |
| P2 | Worker-side async (explore) |

### Architecture Changes

**Components to modify:**

1. `Worker` struct — Add request arrays, async send/receive methods
2. `QueueScheduler` — Add `waitQueueAsync()`, `testQueueAsync()`
3. `IrregularRoutines` — Restructure main loop for async

**Build order:**
- Phase 11: Async infrastructure in Worker/QueueScheduler
- Phase 12: Integration with IrregularRoutines
- Phase 13: Measurement, tuning, validation

### Critical Pitfalls

| Pitfall | Prevention |
|---------|------------|
| Buffer modification before wait | Dedicated buffers per worker |
| Request handle leaks | Track all requests, always wait/test |
| Deadlock from misordered ops | Post receives before sends |
| CM dependency breakage | Maintain special CM handling |
| No actual overlap | Profile overlap effectiveness |

### Batching Strategy (for v2.2+)

**Recommended: Dependency-Aware Batching**

1. Separate particles into "simple" and "CM-dependent"
2. Batch simple particles freely (64-128 per batch)
3. Process CM particles after neighbors complete

Expected impact: 100-200x message reduction, 10-30% speedup

## Research Files

- `STACK.md` — MPI functions, patterns, request management
- `FEATURES.md` — Table stakes, differentiators, anti-features
- `ARCHITECTURE.md` — Component changes, build order, data flow
- `PITFALLS.md` — Common mistakes, prevention, recovery
- `BATCHING.md` — Strategies for v2.2+ implementation

## Sources

- [LLNL HPC Tutorials - Non-blocking MPI](https://hpc-tutorials.llnl.gov/mpi/non_blocking/)
- [ENCCS Intermediate MPI](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/)
- [Optimizing Computation-Communication Overlap](https://www.cs.umd.edu/~bhatele/pubs/pdf/2019/ics2019.pdf)
- [MPI Tutorial](https://mpitutorial.com/)

---
*Research completed: 2026-01-18*
