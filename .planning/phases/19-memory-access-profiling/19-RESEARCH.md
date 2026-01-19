# Phase 19 Research: Memory Access Profiling

**Generated:** 2026-01-18
**Phase:** 19 - Memory Access Profiling
**Goal:** Measure cache behavior to determine if memory access is limiting performance

## Research Summary

Phase 19 requires integrating hardware performance counters to measure cache miss rates and estimate memory bandwidth. Two main approaches exist: PAPI (Performance API) and direct perf_event_open() syscall.

## Key Findings

### 1. perf_event_open() vs PAPI

**perf_event_open() (Recommended)**
- Native Linux syscall, no external dependencies
- Available on all modern Linux kernels (2.6.32+)
- Direct access to hardware performance counters
- No library installation required on HPC cluster
- Requires `<linux/perf_event.h>` and `<linux/hw_breakpoint.h>`

**PAPI**
- Portable abstraction layer over perf_event
- Requires separate installation/module load
- May not be available on all HPC clusters
- Adds external dependency to Makefile

**Decision:** Use perf_event_open() directly for portability. PAPI can be added as optional enhancement later.

### 2. Available Hardware Counters

Standard Intel Skylake-AVX512 counters (from Makefile: `-march=skylake-avx512`):

| Counter | Config Value | Description |
|---------|-------------|-------------|
| L1D misses | PERF_COUNT_HW_CACHE_L1D + READ + MISS | L1 data cache read misses |
| L1D accesses | PERF_COUNT_HW_CACHE_L1D + READ + ACCESS | L1 data cache read accesses |
| LL misses | PERF_COUNT_HW_CACHE_LL + READ + MISS | Last-level cache misses |
| LL accesses | PERF_COUNT_HW_CACHE_LL + READ + ACCESS | Last-level cache accesses |

Note: "LL" (Last Level) refers to L3 on Skylake. L2 requires Intel native events.

### 3. Cache Line Size and Bandwidth Estimation

**Cache Line:** 64 bytes on Skylake

**Bandwidth formula:**
```
bandwidth_GB_s = (cache_misses * cache_line_size) / time_seconds / 1e9
```

**Peak theoretical bandwidth:**
- Skylake-AVX512: ~100-150 GB/s (depends on memory channels)
- Typical achieved: 60-80% of peak

### 4. perf_event_open() C++ Integration

```cpp
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>

// Start counter
int perf_fd = syscall(__NR_perf_event_open, &pe, 0, -1, -1, 0);
ioctl(perf_fd, PERF_EVENT_IOC_RESET, 0);
ioctl(perf_fd, PERF_EVENT_IOC_ENABLE, 0);

// ... code to measure ...

// Stop and read
ioctl(perf_fd, PERF_EVENT_IOC_DISABLE, 0);
read(perf_fd, &count, sizeof(count));
close(perf_fd);
```

### 5. MPI Considerations

- Hardware counters are per-core, not per-process
- Each MPI rank measures its own core's counters
- Counters must be opened/closed within the process
- No MPI aggregation needed for raw counts (each rank reports its own)

### 6. Roofline Analysis Basics

**Operational Intensity:**
```
OI = FLOPs / Bytes_transferred
```

**Classification:**
- If `OI < machine_balance`: Memory-bound
- If `OI > machine_balance`: Compute-bound

**Machine balance** for Skylake:
- Peak FLOPs: ~2 TFLOPs/s (with AVX-512)
- Peak bandwidth: ~100 GB/s
- Balance point: ~20 FLOPs/byte

For force calculation:
- ~40-60 FLOPs per particle pair
- ~192 bytes loaded per pair (3 particles × 8 doubles × 8 bytes)
- Estimated OI: ~0.3 FLOPs/byte → **Memory-bound**

### 7. Integration Points in Existing Code

**profiler.h:**
- Add `CacheStats` struct similar to `ParticleTypeStats`
- Add `PROFILE_CACHE_START/STOP` macros
- Conditionally compiled with `#ifdef PERFORMANCETRACE`

**compute_acceleration.cpp:**
- Wrap neighbor loop with cache measurement
- Already has `particle_start_time` for per-particle timing (Phase 15)
- Add cache counters around same region

**Makefile:**
- No changes needed (perf_event is kernel-provided)

### 8. Fallback Strategy

If perf_event_open() fails (permissions, container, etc.):
1. Log warning to stderr
2. Report "-1" for cache metrics
3. Continue with timing-only profiling

This avoids breaking the build on systems without counter access.

## Recommendations

### Plan Structure (3 plans, 2 waves)

**Wave 1:**
- Plan 19-01: CacheStats infrastructure (struct, perf_event wrapper, macros)
- Plan 19-02: Instrument force loop with cache measurement

**Wave 2:**
- Plan 19-03: Memory statistics output (console/CSV/JSON, roofline classification)

### Risk Mitigation

| Risk | Mitigation |
|------|------------|
| perf_event fails | Graceful fallback with warning |
| Counter multiplexing | Measure L1D and LL in separate groups |
| Permission denied | CAP_SYS_ADMIN or paranoid=0 required |
| MPI interference | Each rank uses own counters |

## References

- [perf_event_open(2) man page](https://man7.org/linux/man-pages/man2/perf_event_open.2.html)
- [Linux perf Examples](https://www.brendangregg.com/perf.html)
- [Cache miss analysis with perf](https://www.baeldung.com/linux/analyze-cache-misses)
- [PAPI overview](https://hpc.llnl.gov/software/development-environment-software/papi-performance-application-programming-interface)
- [CPP_LPE_wrap library](https://github.com/jasonspencer/CPP_LPE_wrap)

---
*Research completed: 2026-01-18*
