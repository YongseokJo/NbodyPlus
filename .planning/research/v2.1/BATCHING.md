# Batching Research: Strategies for v2.2+

This document captures research on batching strategies for future implementation.
Batching is out of scope for v2.1 but needs investigation due to complexity.

## The Batching Challenge

Current: ~105M messages for ~140K particle evaluations (750 messages/particle)

Target: Reduce to ~500K-1M messages (100-200x reduction)

### Why Batching is Complex for ABYSS

1. **CM Particle Dependencies**
   - CM particles need neighbor forces computed first
   - Can't batch CM particles with their neighbors
   - Need dependency tracking across batches

2. **Variable Work Per Particle**
   - Particle with 100 neighbors: 10x compute time vs 10 neighbors
   - Fixed batches cause load imbalance
   - Need adaptive batch sizing

3. **Callback Mechanism**
   - Current: one completion per particle
   - Batched: one completion per batch
   - Skip list updates, termination checks need refactoring

## Batching Strategies Considered

### Strategy A: Static Chunking

**Approach:** Divide particle list into fixed-size chunks (e.g., 64 particles/chunk)

**Pros:**
- Simple to implement
- Predictable memory usage
- Easy to reason about

**Cons:**
- Load imbalance (some chunks heavier than others)
- CM dependencies cross chunk boundaries
- Can't adapt to varying neighbor counts

**Verdict:** Too simple for ABYSS complexity

### Strategy B: Work-Weighted Batching

**Approach:** Batch particles so each batch has similar total work (sum of neighbor counts)

**Pros:**
- Better load balance
- Adapts to particle distribution

**Cons:**
- Requires pre-computing work estimates
- CM dependencies still complex
- Batch composition changes each timestep

**Verdict:** Better, but CM handling unclear

### Strategy C: Dependency-Aware Batching

**Approach:** 
1. Separate particles into "simple" (no CM deps) and "CM-dependent"
2. Batch simple particles freely
3. Process CM particles after their neighbors complete

**Pros:**
- Respects physics dependencies
- Simple particles get full batching benefit
- CM particles processed correctly

**Cons:**
- Two-phase processing
- CM particles may still be fine-grained
- Need to track which neighbors are done

**Verdict:** Most promising for ABYSS

### Strategy D: Worker-Pull Model

**Approach:** Workers pull batches from shared work queue rather than root pushing

**Pros:**
- Natural load balancing
- Workers never idle if work exists
- Scales with worker count

**Cons:**
- Requires shared queue (MPI-3 RMA or shared memory)
- Synchronization complexity
- Major architectural change

**Verdict:** Good long-term, too invasive for v2.2

## Recommended Approach for v2.2

**Dependency-Aware Batching (Strategy C)** with these phases:

### Phase 1: Simple Particle Batching
- Identify particles without CM dependencies
- Batch into chunks of 64-128 particles
- Send batch, receive batch completion
- ~60-80% of particles, immediate benefit

### Phase 2: CM-Aware Processing
- Track which simple particles are done
- Release CM particles when neighbors complete
- May still be fine-grained for CM particles
- Correctness preserved

### Phase 3: Adaptive Batch Sizing
- Measure work per particle (neighbor count as proxy)
- Adjust batch size to balance load
- Target: each batch ~same total work

## Implementation Considerations

### Message Format Change

**Current:**
```cpp
struct Queue {
    task_name_t task;
    int pid;         // single particle
    double next_time;
};
```

**Batched:**
```cpp
struct BatchQueue {
    task_name_t task;
    int batch_size;
    int pids[MAX_BATCH_SIZE];  // multiple particles
    double next_time;
};
```

### Callback Aggregation

**Current:** One `MPI_Recv` per particle completion

**Batched:** 
- One `MPI_Recv` per batch
- Worker sends back completion status for all particles in batch
- Root updates all particles at once

### Skip List Integration

**Current:** Update skip list per particle

**Batched:**
- Batch skip list updates
- Or: collect updates, apply in bulk after batch

## Estimated Impact

| Metric | Current | With Batching |
|--------|---------|---------------|
| Messages/interval | ~105M | ~1M |
| Message size | ~24 bytes | ~1-2 KB |
| Total bandwidth | ~2.5 GB | ~1-2 GB |
| MPI overhead | High | Low |

Expected improvement: 10-30% reduction in irregular force time

## Dependencies on v2.1

Batching builds on async infrastructure:
- Request management (v2.1) → extended for batch requests
- Profiling (v2.1) → measure batch effectiveness
- Completion handling (v2.1) → batch completion callbacks

## Open Questions for v2.2

1. Optimal batch size? (Need to measure on real workloads)
2. How to handle CM particles efficiently?
3. Should workers aggregate results or send individually?
4. Memory overhead of batch buffers?

## Sources

- [MPI Tutorial - Scatter/Gather](https://mpitutorial.com/tutorials/mpi-scatter-gather-and-allgather/)
- [Dynamic Load Balancing in MPI](https://link.springer.com/chapter/10.1007/978-3-540-77704-5_10)
- [Task Aggregation in HPC](https://ieeexplore.ieee.org/document/5470464/)

---
*Research completed: 2026-01-18*
*For implementation in v2.2+*
