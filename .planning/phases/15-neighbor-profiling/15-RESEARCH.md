# Phase 15 Research: Neighbor Profiling

## Overview

This phase adds per-particle neighbor count profiling to understand load imbalance in the irregular force calculation.

## Current State Analysis

### Where Neighbor Counts Exist

1. **`Particle::num_neighbors`** — Stored per particle, represents the current neighbor count
2. **`compute_acceleration_irr()`** in `src/Particle/compute_acceleration.cpp`:
   - Line 66: Checks `this->num_neighbors == 0` for early return
   - Line 103: Tracks `neighbor_pairs` (batch.count)
   - Line 128: `neighbor_pairs = batch.count`
   - Line 152: `PROFILE_WORK(TimerID::IrregularPairsEvaluated, neighbor_pairs)`
   - Line 155: Tracks `cm_pairs` separately
   - Line 203: `PROFILE_WORK(TimerID::IrregularPairsEvaluated, cm_pairs)`

3. **CM Particles** — Some particles are center-of-mass particles that trigger additional work in the CM loop (lines 154-203)

### Existing Profiler Infrastructure

From `src/profiler.h`:

1. **TimerStats** already tracks:
   - `total_ns`, `count`, `min_ns`, `max_ns`, `last_ns`
   - `interval_total_ns`, `interval_count`, `interval_min_ns`, `interval_max_ns`
   - `work_units`, `interval_work_units`
   - Optional `Histogram` for distributions

2. **PROFILE_WORK macro** — Already used to track neighbor pairs: `PROFILE_WORK(TimerID::IrregularPairsEvaluated, neighbor_pairs)`

3. **Histogram class** — Already supports:
   - `record(duration_ns)` for timing distributions
   - `getPercentile(p)`, `getMedian()`, `toJSON()`
   - Log-scale buckets from <1us to >10s

4. **AggregatedStats** — MPI aggregation with min/max/avg across ranks

### Key Insight: What's Missing

The current profiler tracks **total neighbor pairs evaluated per interval**, but not:
- Per-particle neighbor count distribution
- Which particles are outliers
- Correlation between neighbor count and compute time
- Per-particle timing

## Implementation Strategy

### 1. New Profiler Structures Needed

```cpp
// Online statistics calculator (Welford's algorithm)
struct OnlineStats {
    long long count = 0;
    double mean = 0.0;
    double M2 = 0.0;      // For variance calculation
    long long min_val = LLONG_MAX;
    long long max_val = 0;

    void update(long long value) {
        count++;
        double delta = value - mean;
        mean += delta / count;
        M2 += delta * (value - mean);
        min_val = std::min(min_val, value);
        max_val = std::max(max_val, value);
    }

    double variance() const { return count > 1 ? M2 / (count - 1) : 0.0; }
    double stddev() const { return std::sqrt(variance()); }
};

// Per-particle profiling data
struct ParticleProfileData {
    int neighbor_count;
    long long compute_time_ns;
    bool is_outlier;
};
```

### 2. Where to Instrument

**Location 1: `compute_acceleration_irr()` entry**
```cpp
// Record neighbor count
int total_neighbors = this->num_neighbors;
PROFILE_NEIGHBOR_COUNT(total_neighbors);
```

**Location 2: Full function timing (already exists as IrregularForce)**
- Already tracked at caller level in `irregular_routines.cpp`
- Need per-particle granularity

**Location 3: After CM loop completes**
```cpp
int cm_extra_neighbors = cm_pairs;
total_neighbors += cm_extra_neighbors;
```

### 3. Outlier Detection Strategy

Using 2σ threshold with online calculation:
1. Track running mean and stddev during execution
2. After interval completes, flag particles where `count > mean + 2*stddev`
3. Store outlier info for analysis output

### 4. Correlation Calculation

For correlation between neighbor count (X) and compute time (Y):
```cpp
// Pearson correlation via online algorithm
struct CorrelationTracker {
    long long n = 0;
    double sum_x = 0, sum_y = 0;
    double sum_xy = 0;
    double sum_x2 = 0, sum_y2 = 0;

    void update(double x, double y) {
        n++;
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
        sum_y2 += y * y;
    }

    double correlation() const {
        if (n < 2) return 0.0;
        double num = n * sum_xy - sum_x * sum_y;
        double den = std::sqrt((n * sum_x2 - sum_x * sum_x) *
                               (n * sum_y2 - sum_y * sum_y));
        return den > 0 ? num / den : 0.0;
    }
};
```

### 5. Integration Points

**In `src/profiler.h`:**
- Add `OnlineStats neighbor_stats`
- Add `CorrelationTracker neighbor_time_correlation`
- Add new TimerID for per-particle timing
- Extend `writeCSV()` to include neighbor statistics

**In `src/Particle/compute_acceleration.cpp`:**
- Add timing around full `compute_acceleration_irr()`
- Record neighbor count at function entry
- Record total neighbors (including CM) at function exit

**In `src/irregular_routines.cpp`:**
- Per-particle timing already exists via worker dispatch
- Need to associate timing with neighbor count

### 6. Output Format

**CSV additions:**
```
NeighborCount_min,NeighborCount_max,NeighborCount_mean,NeighborCount_stddev,
NeighborCount_outliers,NeighborTime_correlation
```

**Histogram for neighbor counts:**
- Use existing Histogram class with integer buckets
- Buckets: 0-10, 10-50, 50-100, 100-200, 200-500, 500+

## Overhead Analysis

**Expected overhead:**
- OnlineStats::update() — 5-10 ns per call (simple arithmetic)
- CorrelationTracker::update() — 10-15 ns per call
- ~140K particles per interval → ~2-3 ms overhead
- Acceptable: <2% of 130s wall time

**Mitigation:**
- Use integer arithmetic where possible
- Batch updates if needed
- Conditional compilation via `LOAD_BALANCE_PROFILING` flag

## Files to Modify

1. **`src/profiler.h`**
   - Add OnlineStats, CorrelationTracker classes
   - Add TimerID::NeighborCount, TimerID::ParticleComputeTime
   - Add neighbor statistics output methods

2. **`src/Particle/compute_acceleration.cpp`**
   - Add per-particle timing around compute_acceleration_irr()
   - Record neighbor count at entry
   - Record total neighbors at exit

3. **`src/irregular_routines.cpp`**
   - Associate per-particle timing with neighbor counts
   - Aggregate statistics per interval

## Success Criteria

1. [ ] Neighbor count tracked for every particle
2. [ ] Per-interval statistics (min/avg/max/stddev) output to CSV
3. [ ] Outlier particles identified (>2σ neighbor count)
4. [ ] Correlation coefficient computed between neighbor count and compute time
5. [ ] Overhead < 5% of runtime
