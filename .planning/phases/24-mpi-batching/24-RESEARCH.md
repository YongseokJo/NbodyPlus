# Phase 24 Research: MPI Message Batching

**Completed:** 2026-01-19
**Phase:** 24 (MPI Message Batching)

---

## Current Architecture Analysis

### Message Flow (Current)

```
For each particle task:
  Root → Worker: Queue struct (1 MPI_Send)
  Worker → Root: completion int (1 MPI_Send)

With 10.3M tasks/interval = 20.6M MPI messages/interval
```

### Key Data Structures

**Queue struct** (`src/queue.h`):
```cpp
struct Queue {
    task_name_t task;       // int8_t - Task type to execute
    int pid;                // Particle ID (target of task)
    double next_time;       // Next time for time-based tasks
};
```

**Worker struct** (`src/worker.h`):
- `send_task(Queue& queue)` - Sends single Queue via MPI_Send
- `callback()` - Receives completion via MPI_Recv
- Has queue buffer (`queues[MAX_QUEUE]`) but sends one at a time

### Dispatch Flow

**Root side** (`queue_scheduler.h`, `irregular_routines.cpp`):
1. `queue_scheduler.initializeIrr()` - Sets up particle list
2. `queue_scheduler.assignQueueAuto()` - Assigns one task per free worker
3. `worker->addQueue(queue)` - Adds to worker's queue
4. `worker->runQueue()` → `worker->send_task()` → `MPI_Send()`
5. `queue_scheduler.waitQueue()` - Waits for any completion
6. `worker->callback()` → `MPI_Recv()` completion

**Worker side** (`worker_routines.cpp`):
1. `MPI_Recv(&queue, ...)` - Receive single Queue
2. Switch on task type, process particle
3. `MPI_Isend(&send_value, ...)` - Signal completion

### Bottleneck Analysis

From Phase 23 analysis:
- **MPI overhead:** 28% of wall time (4.62s/interval)
- **MPISend:** 4.03s (24.1%)
- **MPIRecv:** 0.59s (3.5%)
- **QueueWait:** 0.50s (3.0%)
- **Starvation events:** 2.4M per interval

Root can't dispatch fast enough → workers starve despite work available.

---

## Implementation Approach

### Batched Queue Structure

```cpp
// New batched queue structure
struct BatchedQueue {
    task_name_t task;           // All tasks in batch share same type
    int count;                  // Number of particles (1 to MAX_BATCH_SIZE)
    int pids[MAX_BATCH_SIZE];   // Particle IDs
    double next_time;           // Shared next_time (irregular tasks)
};

// Config
constexpr int MAX_BATCH_SIZE = 100;  // Configurable upper bound
extern int batch_size;                // Runtime configurable (default: 50)
```

### MPI Type Registration

Need to create MPI_Datatype for BatchedQueue in MPI initialization.

### Root-Side Changes

**queue_scheduler.h modifications:**

1. Add batch accumulation buffer:
```cpp
std::vector<int> _batch_buffer;
task_name_t _batch_task;
double _batch_next_time;
```

2. Add batch dispatch method:
```cpp
void sendBatch(Worker* worker) {
    BatchedQueue batch;
    batch.task = _batch_task;
    batch.count = std::min((int)_batch_buffer.size(), batch_size);
    batch.next_time = _batch_next_time;
    for (int i = 0; i < batch.count; i++) {
        batch.pids[i] = _batch_buffer.back();
        _batch_buffer.pop_back();
    }
    MPI_Send(&batch, 1, batched_queue_type_mpi, worker->rank, BATCH_QUEUE_TAG, MPI_COMM_WORLD);
}
```

3. Modify `assignQueueAuto()` to batch:
   - Fill batch buffer from queue_list
   - When free worker available, send full batch
   - Handle end-of-queue partial batches

### Worker-Side Changes

**worker_routines.cpp modifications:**

1. Add batched receive path:
```cpp
case BATCH_QUEUE_TAG:
    MPI_Recv(&batch, 1, batched_queue_type_mpi, ROOT, BATCH_QUEUE_TAG, ...);
    for (int i = 0; i < batch.count; i++) {
        ptcl_id = batch.pids[i];
        // Process particle based on batch.task
    }
    // Single completion message with count
    MPI_Isend(&batch.count, 1, MPI_INT, ROOT, TERMINATE_TAG, ...);
    break;
```

2. Keep single-task path for compatibility (CM particles, etc.)

### Completion Handling

Root callback needs to track completed count vs single task:
- Batch completion: `_completed_queues += return_value` (count)
- Single completion: `_completed_queues += 1`

---

## Design Decisions

### 1. Batch Size Strategy

- **Default:** 50 particles per batch
- **Range:** 10-100 (configurable via config file)
- **Rationale:** Balance between reduced messages and latency

### 2. Task Homogeneity

- Only batch same-type tasks (TASK_IRR_FORCE with TASK_IRR_FORCE)
- CM particle tasks (TASK_AR_INTEGRATION) stay single-dispatch
- Reason: Different tasks have different processing needs

### 3. Partial Batches

- Send partial batch when:
  - Queue exhausted (< batch_size remaining)
  - All workers busy (flush to avoid starvation)
- Track partial vs full for profiling

### 4. CM Particle Handling

- Keep existing single-dispatch for TASK_AR_INTEGRATION
- CM particles go to specific workers (cm_particle_worker_map)
- Batching regular particles only

### 5. MPI Tags

- New tag: `BATCH_QUEUE_TAG = 5` for batched messages
- Keep `QUEUE_TAG = 4` for single messages (backward compat)

---

## Risk Analysis

| Risk | Mitigation |
|------|------------|
| Batch accumulation delays dispatch | Flush on timer or worker availability |
| Worker compute time variance | Profile batch processing time |
| Memory overhead | Fixed-size batch buffer (bounded) |
| Correctness (energy) | Verify energy conservation unchanged |
| Partial batch inefficiency | Track and tune threshold |

---

## Expected Outcomes

From Amdahl's Law (Phase 23):
- **Conservative (50% reduction):** 16% speedup
- **Expected (67% reduction):** 21% speedup
- **Optimistic (75% reduction):** 27% speedup

Target: Reduce 10.3M messages → 200K messages (50x reduction)

---

## File Modifications Required

| File | Changes |
|------|---------|
| `src/queue.h` | Add BatchedQueue struct |
| `src/global.h` | Add BATCH_QUEUE_TAG, batched_queue_type_mpi |
| `src/queue_scheduler.h` | Batch accumulation, sendBatch() |
| `src/worker_routines.cpp` | Batched receive/process loop |
| `src/worker.h` | Optional batch_callback() |
| `src/mpi_types.cpp` | Register BatchedQueue MPI type |
| `tests/test4/config.toml` | Add batch_size config option |

---

## Testing Strategy

1. **Unit test:** Verify BatchedQueue serialization
2. **Integration test:** Run with batch_size=1 (should behave like current)
3. **Performance test:** Compare batch_size=10, 50, 100
4. **Correctness test:** Energy conservation check

---

*Research completed: 2026-01-19*
