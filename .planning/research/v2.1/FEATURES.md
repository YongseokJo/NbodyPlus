# Features Research: MPI Optimization Techniques

## Table Stakes (Must Have)

These are essential for any async MPI implementation:

### 1. Non-Blocking Send/Receive
- Replace `MPI_Send` with `MPI_Isend`
- Replace `MPI_Recv` with `MPI_Irecv`
- Proper request handle management

### 2. Completion Handling
- `MPI_Waitany` for processing first available result
- `MPI_Testany` for non-blocking completion checks
- Proper status handling or explicit `MPI_STATUS_IGNORE`

### 3. Buffer Safety
- No buffer modification between Isend and Wait
- Clear buffer ownership semantics
- Separate buffers per outstanding request

### 4. Request Lifecycle Management
- Initialize requests to `MPI_REQUEST_NULL`
- Track active vs completed requests
- Clean reuse of request slots

### 5. Profiling Integration
- Measure async operation timing
- Track overlap effectiveness
- Compare wait time before/after

## Differentiators (Advanced Techniques)

These provide extra benefit beyond basic async:

### 1. Communication-Computation Overlap
- Post sends, do local work, then wait
- Achieves 3-88% improvement in benchmarks
- Requires careful code restructuring

### 2. Pipelining
- Multiple outstanding sends before any wait
- Reduces idle time between task completions
- Increases in-flight message count

### 3. Adaptive Wait Strategy
- `MPI_Testany` loop during computation phases
- Fall back to `MPI_Waitany` when no local work
- Balances progress checking with work efficiency

### 4. Worker-Side Async (Optional)
- Workers use `MPI_Irecv` to pre-post receives
- Reduces receive latency
- More complex implementation

### 5. Request Batching
- Group multiple logical operations
- Reduce MPI call overhead
- Prerequisite for future message batching

## Anti-Features (Do NOT Build)

These are explicitly out of scope for v2.1:

### 1. Message Batching/Aggregation
- Grouping multiple particles per message
- Complex due to CM dependencies
- Research only, defer implementation to v2.2

### 2. Persistent Requests
- `MPI_Send_init` / `MPI_Recv_init` / `MPI_Start`
- Better for fixed patterns, ABYSS has variable patterns
- Overhead not justified

### 3. One-Sided Communication (RMA)
- `MPI_Put` / `MPI_Get` / `MPI_Accumulate`
- Major architectural change
- Different programming model

### 4. MPI-3 Neighborhood Collectives
- `MPI_Neighbor_alltoall` etc.
- Requires graph communicator setup
- Overkill for current root/worker pattern

### 5. Dynamic Process Management
- `MPI_Comm_spawn` / Sessions
- Not supported on most batch systems
- Wrong abstraction level

## Feature Priority for v2.1

| Priority | Feature | Rationale |
|----------|---------|-----------|
| P0 | Non-blocking send/recv | Foundation for all other features |
| P0 | Request management | Required for correctness |
| P0 | Completion handling | Core functionality |
| P1 | Profiling integration | Measure improvement |
| P1 | Comm-compute overlap | Main performance benefit |
| P2 | Worker-side async | Explore if beneficial |
| P2 | Adaptive wait strategy | Optimization |

## Sources

- [LLNL HPC Tutorials](https://hpc-tutorials.llnl.gov/mpi/non_blocking/)
- [ENCCS Intermediate MPI](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/)
- [Optimizing Computation-Communication Overlap](https://www.cs.umd.edu/~bhatele/pubs/pdf/2019/ics2019.pdf)

---
*Research completed: 2026-01-18*
