# Stack Research: MPI Async Communication

## Recommended Approach

**Non-blocking point-to-point with request management**

For ABYSS v2.1, use the standard MPI non-blocking primitives:

| Function | Purpose | Use Case |
|----------|---------|----------|
| `MPI_Isend` | Non-blocking send | Root sends tasks to workers |
| `MPI_Irecv` | Non-blocking receive | Root receives completion signals |
| `MPI_Waitany` | Wait for any request | Process first available result |
| `MPI_Testany` | Non-blocking test | Check for completions while doing work |
| `MPI_Waitall` | Wait for all requests | Synchronization points |

### Core Pattern

```cpp
// Root-side async pattern
MPI_Request send_requests[MAX_WORKERS];
MPI_Request recv_requests[MAX_WORKERS];

// Post all receives first (avoid deadlock)
for (int i = 0; i < num_workers; i++) {
    MPI_Irecv(&results[i], ..., worker_rank[i], ..., &recv_requests[i]);
}

// Send tasks asynchronously
for (int i = 0; i < num_tasks; i++) {
    MPI_Isend(&tasks[i], ..., target_worker, ..., &send_requests[i]);
}

// Process completions as they arrive
while (completed < total) {
    int index;
    MPI_Waitany(num_workers, recv_requests, &index, MPI_STATUS_IGNORE);
    // Process result from worker[index]
    // Post new receive for next result
    // Send next task if available
}
```

## Request Management Strategy

**Fixed-size request arrays** (recommended for ABYSS):

- Allocate `MPI_Request` arrays sized to `num_workers`
- Each worker has dedicated slots in send/recv request arrays
- Reuse slots after `MPI_Wait` completes
- Avoid dynamic allocation in hot loops

**Why not persistent requests (`MPI_Send_init`/`MPI_Recv_init`):**

- Better for repeated identical communication patterns
- ABYSS has variable-length particle lists per timestep
- Overhead of `MPI_Start` not worthwhile for variable patterns

## Overlap Strategy

To overlap communication with computation:

1. **Post receives early** — Before any sends
2. **Send multiple tasks** — Don't wait after each send
3. **Do local work** — Between sends and waits
4. **Use `MPI_Testany`** — For opportunistic completion checks

```cpp
// Overlap pattern
MPI_Isend(...);
MPI_Isend(...);
// Do local computation here
while (local_work_remaining && !all_complete) {
    do_chunk_of_local_work();
    MPI_Testany(n, requests, &index, &flag, &status);
    if (flag) process_completion(index);
}
MPI_Waitall(remaining, ...);
```

## What NOT to Do

| Anti-pattern | Why It's Bad | Instead |
|--------------|--------------|---------|
| Modify send buffer before wait | Undefined behavior, data corruption | Wait before reusing buffer |
| `MPI_Request_free` on active requests | Cannot check completion, errors become fatal | Always wait/test first |
| Blocking send then blocking recv | Potential deadlock with many messages | Use non-blocking or mixed |
| Single Isend → Wait → repeat | No overlap benefit | Batch multiple sends |
| Ignoring request status | Can't detect errors | Check status or use `MPI_STATUS_IGNORE` explicitly |

## MPI Progress Considerations

Non-blocking operations need "progress" to complete:

- **MPICH-style**: Progress happens during MPI calls (Wait, Test, etc.)
- **OpenMPI**: Similar, but may have async progress threads
- **Hardware offload**: Some InfiniBand setups allow true background progress

For ABYSS:
- Call `MPI_Testany` periodically during computation
- Don't assume background progress without testing
- Profile to verify overlap is actually happening

## Sources

- [LLNL HPC Tutorials - Non-blocking MPI](https://hpc-tutorials.llnl.gov/mpi/non_blocking/)
- [ENCCS Intermediate MPI](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/)
- [MPI Send Modes](https://www.mcs.anl.gov/research/projects/mpi/sendmode.html)
- [Cornell Virtual Workshop - Nonblocking Communication](https://cvw.cac.cornell.edu/mpip2p/nonblocking-communication/wait-test-free)

---
*Research completed: 2026-01-18*
