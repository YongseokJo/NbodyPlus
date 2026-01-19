# Phase 17: Worker Distribution - Summary

## Completion Date
2026-01-18

## Goal
Measure how work is distributed across workers to identify imbalance.

## Requirements Addressed
- WRKR-01: Count particles processed per worker per interval ✓
- WRKR-02: Track compute time per worker per interval ✓
- WRKR-03: Identify "heavy" particles (>2σ compute time) ✓
- WRKR-04: Compute load balance ratio (max_worker_time / avg_worker_time) ✓

## Plans Executed

### Wave 1
- **17-01**: WorkerDistributionTracker infrastructure
  - Added `HeavyParticleInfo` struct
  - Added `WorkerDistributionTracker` class with particle/time tracking
  - Added Profiler methods for worker tracking
  - Added PROFILE_WORKER_* macros

### Wave 2 (Parallel)
- **17-02**: Per-worker particle count instrumentation
  - Instrumented assignQueueAuto(), assignQueueAutoRegularList(), assignWorker()
  - Added worker tracking initialization in IrregularRoutines

- **17-03**: Per-worker compute time and heavy particle tracking
  - Extended compute_acceleration_irr() to record worker compute time
  - Implemented heavy particle detection (>2σ threshold)

### Wave 3
- **17-04**: Worker distribution statistics output
  - Console: Worker distribution stats with load balance ratio
  - CSV: WorkerParticles_*, WorkerComputeTime_*, LoadBalanceRatio columns
  - JSON: worker_distribution object with nested structures

## Key Deliverables

### New Classes/Structs
- `HeavyParticleInfo`: particle_id, worker_rank, compute_time_ns, neighbor_count
- `WorkerDistributionTracker`: Per-worker particle counts and compute times

### New Macros
- `PROFILE_WORKER_ASSIGNMENT(worker_rank)` - Track particle assignments
- `PROFILE_WORKER_COMPUTE(worker_rank, compute_ns)` - Track compute time
- `PROFILE_HEAVY_PARTICLE(pid, worker, time_ns, neighbors)` - Track outliers

### New Profiler Methods
- `initializeWorkerTracking(num_workers)` - Setup worker tracking
- `recordWorkerAssignment(worker_rank)` - Count particle assignment
- `recordWorkerComputeTime(worker_rank, compute_ns)` - Accumulate compute time
- `recordHeavyParticle(...)` - Record outlier particles
- `getCurrentParticleWorkerRank()` - Get current particle's worker

### Output Enhancements
- Console: Load balance ratio with max time worker identification
- CSV: 11 new columns for worker distribution metrics
- JSON: worker_distribution object with particles_per_worker, compute_time, heavy_particles

## Files Modified
- `src/profiler.h` - WorkerDistributionTracker class, output methods
- `src/queue_scheduler.h` - PROFILE_WORKER_ASSIGNMENT calls
- `src/irregular_routines.cpp` - Worker tracking initialization
- `src/Particle/compute_acceleration.cpp` - Worker compute time recording

## Metrics Available
- Particles per worker: min, max, mean, stddev
- Compute time per worker: min, max, mean (seconds)
- Load balance ratio: max_worker_time / avg_worker_time
- Heavy particles: Count and details (PID, worker, time, neighbors)
- Imbalance warning: When ratio > 1.5
