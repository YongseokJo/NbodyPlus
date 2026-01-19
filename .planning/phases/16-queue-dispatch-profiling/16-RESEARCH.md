# Research: Phase 16 — Queue Dispatch Profiling

## Overview

Phase 16 instruments the MPI queue dispatch system to understand if root-side dispatch overhead or worker starvation contributes to load imbalance.

## Current Queue Scheduler Architecture

### Key Components

**QueueScheduler class** (`src/queue_scheduler.h`):
- `assignQueueAuto()` — Assigns particles from queue to free workers
- `runQueueAuto()` — Sends queued tasks to workers via MPI
- `waitQueue(type)` — Blocking (0) or non-blocking (1) wait for worker completion
- `callback(worker)` — Handles worker completion and reassignment

**Worker struct** (`src/worker.h`):
- `run_queue()` → `send_task()` — Sends task via MPI_Send
- `callback()` — Receives completion via MPI_Recv
- `num_queues`, `current_queue` — Queue management

### Existing Profiling

Already instrumented with:
- `TimerID::QueueAssign` — Time in assignQueueAuto()
- `TimerID::QueueRun` — Time in runQueueAuto()
- `TimerID::QueueCallback` — Time in callback()
- `TimerID::QueueWait` — Time in MPI_Probe wait
- `TimerID::MPISend` — Time in MPI_Send
- `TimerID::MPIRecv` — Time in MPI_Recv

### Data Flow

```
Root:                           Worker:
assignQueueAuto()               (waiting)
   └→ adds to WorkersToGo
runQueueAuto()
   └→ worker->runQueue()
      └→ MPI_Send(queue)  →→→   MPI_Recv(queue)
                                compute_acceleration_irr()
waitQueue(0)                    MPI_Send(TERMINATE_TAG)
   └→ MPI_Probe()         ←←←
callback()
   └→ MPI_Recv()
```

## Requirements Analysis

### QUEUE-01: Dispatch Latency (Worker Side)

**Goal:** Measure time from worker completion to receiving next task.

**Implementation:**
- Already have `WorkerRecvWait` timer in worker loop (from Phase 8)
- Need per-dispatch tracking, not just cumulative
- Add `OnlineStats` tracker for dispatch latency distribution

**Location:** Worker-side code (not in queue_scheduler.h — that's root only)
- Need to check `src/irregular_routines.cpp` for worker loop

### QUEUE-02: Root-Side Dispatch Overhead

**Goal:** Track time spent in assign vs waiting.

**Implementation:**
- Already have `QueueAssign` and `QueueWait` timers
- Need per-interval breakdown: assign_time vs wait_time
- Can use existing timer infrastructure — just need output

**New metrics:**
- `assign_time / wait_time` ratio per interval
- If assign_time >> wait_time: dispatch is bottleneck
- If wait_time >> assign_time: workers are bottleneck

### QUEUE-03: Queue Depth Over Time

**Goal:** Sample pending tasks in queue.

**Implementation:**
- `_queue_list.size()` is the pending queue depth
- Sample at regular intervals (e.g., after each batch of assignments)
- Track min/max/avg queue depth per interval

**Add:**
```cpp
class QueueDepthTracker {
    void sample(int depth);  // Called periodically
    int min_depth, max_depth;
    OnlineStats depth_stats;
};
```

### QUEUE-04: Worker Starvation Detection

**Goal:** Detect when worker waits with non-empty queue.

**Definition:** Starvation = worker idle time while `_queue_list.size() > 0`

**Implementation:**
- Track when waitQueue returns with queue still non-empty
- Count starvation events per interval
- Log event details: which worker, queue depth at time, wait duration

## Integration with Phase 15 Infrastructure

Reuse from Phase 15:
- `OnlineStats` class for statistics
- `interval_*` pattern for per-output tracking
- `resetIntervalStats()` integration
- CSV/JSON output patterns

## Key Files to Modify

1. **src/profiler.h**
   - Add QueueDepthTracker class
   - Add queue dispatch statistics to Profiler
   - Add output methods

2. **src/queue_scheduler.h**
   - Add queue depth sampling in assignQueueAuto()
   - Add starvation detection in waitQueue()
   - Track per-dispatch timing

3. **src/irregular_routines.cpp** (if worker loop exists there)
   - Worker-side dispatch latency (if not already in worker loop)

## Plan Structure

| Plan | Focus | Dependencies |
|------|-------|--------------|
| 16-01 | Queue depth tracking infrastructure | None |
| 16-02 | Root-side dispatch timing breakdown | 16-01 |
| 16-03 | Worker starvation detection | 16-01 |
| 16-04 | Queue dispatch output (CSV/JSON/console) | 16-01, 16-02, 16-03 |

## Overhead Considerations

- Queue depth sampling: O(1) per sample, minimal overhead
- Starvation detection: O(1) conditional check, minimal overhead
- All profiling conditional on PERFORMANCETRACE flag
