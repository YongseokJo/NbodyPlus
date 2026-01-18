# Phase 10 Research: Targeted Optimization

**Researched:** 2026-01-17

## Executive Summary

The irregular force calculation in `compute_acceleration_irr()` takes 53.5% of wall time. The hot path is the neighbor loop (lines 105-168) which iterates over neighbors, predicts positions/velocities, computes distance, and accumulates acceleration and jerk.

**Key finding:** The current code structure prevents effective SIMD vectorization due to:
1. Non-contiguous neighbor access (pointer chasing through `particles[neighbors[...]]`)
2. Branch-heavy logic (active particle checks, CM particle handling)
3. Per-particle data layout in Particle struct (AoS for particle data)

**Recommended approach:** Restructure neighbor data access and vectorize the force accumulation inner loops using AVX-512 intrinsics.

---

## 1. Hot Code Analysis

### Primary Bottleneck: `compute_acceleration_irr()` in `src/Particle/compute_acceleration.cpp`

The function has three profiled sections:
- **IrregularNeighborLoop** (lines 105-168): Main neighbor loop
- **IrregularCMLoop** (lines 174-217): Center-of-mass particle loop
- **IrregularCorrection** (lines 235-283): 4th-order correction

### Neighbor Loop Structure (Critical Path)

```cpp
for (int i=0; i<this->num_neighbors; i++) {
    ptcl = &particles[neighbors[this->neighbors_offset + i]];

    if (!ptcl->is_active) { continue; }  // Branch

    r2 = 0.0; vx = 0.0;

    ptcl->predict_particle_second_order(dt, pos_neighbor, vel_neighbor);

    for (int dim=0; dim<DIM; dim++) {  // DIM=3
        x[dim] = pos_neighbor[dim] - pos[dim];
        v[dim] = vel_neighbor[dim] - vel[dim];
        r2 += x[dim]*x[dim];
        vx += v[dim]*x[dim];
    }

    m_r3 = ptcl->mass/(r2*sqrt(r2));

    for (int dim=0; dim<DIM; dim++){
        a_tmp[dim] += m_r3*x[dim];
        adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
    }
}
```

### Vectorization Blockers

| Issue | Location | Impact |
|-------|----------|--------|
| Indirect array access | `particles[neighbors[offset+i]]` | Cache misses, no prefetch |
| Branch in loop | `if (!ptcl->is_active)` | Branch misprediction |
| Method call | `predict_particle_second_order()` | Function call overhead |
| Division by r3 | `m_r3 = mass/(r2*sqrt(r2))` | Expensive FP operation |
| Scattered data | Each neighbor from different memory location | Poor cache utilization |

---

## 2. Data Structure Analysis

### Current Particle Struct (AoS)

From `src/particle.h`:
```cpp
struct Particle {
    int pid;
    int particle_index;
    double position[DIM];       // 24 bytes
    double velocity[DIM];       // 24 bytes
    double mass;                // 8 bytes
    double acc_total[DIM][4];   // 96 bytes
    // ... more fields
};
```

Total particle size: ~400+ bytes. Poor cache utilization when accessing just position/velocity/mass.

### Existing ParticleData SoA Container

From `src/particle_data.h`:
```cpp
class ParticleData {
    double* pos_x_;
    double* pos_y_;
    double* pos_z_;
    double* vel_x_;
    double* vel_y_;
    double* vel_z_;
    double* mass_;
    double* acc_total_[3][4];
    // ...
};
```

**Opportunity:** SoA container already exists! However, it's not used in the irregular force path. The code uses the AoS `particles[]` array.

### Neighbor List Access

```cpp
extern int* neighbors;  // Global neighbor index array
// Access: neighbors[particle.neighbors_offset + i]
```

Neighbors are stored contiguously per-particle, but the particle data is scattered.

---

## 3. Prediction Function Analysis

From `src/particle.h`:
```cpp
template <typename T>
void predict_particle_second_order(double dt, T pos[], T vel[]) {
    dt = dt * enzo_time_step;
    for (int dim = 0; dim < DIM; dim++) {
        pos[dim] = ((acc_total[dim][1] * dt / 3 + acc_total[dim][0]) * dt / 2 + velocity[dim]) * dt + position[dim];
        vel[dim] = (acc_total[dim][1] * dt / 2 + acc_total[dim][0]) * dt + velocity[dim];
    }
}
```

This requires: position[3], velocity[3], acc_total[3][0], acc_total[3][1].

**SIMD opportunity:** Prediction is 3 independent calculations that can be vectorized.

---

## 4. Build System Analysis

From `src/Makefile`:
```makefile
CXXFLAGS_BASE = -O2 -march=native -ftree-vectorize -std=c++11 -Wall -Wextra -Wuninitialized -g
```

- `-march=native`: Compiler can use available SIMD (cluster has AVX-512)
- `-ftree-vectorize`: Auto-vectorization enabled
- `-O2`: Standard optimization level

**Finding:** Compiler already has vectorization enabled, but the code structure prevents it.

---

## 5. AVX-512 Optimization Strategy

### AVX-512 Capabilities

- 8 doubles per register (512 bits / 64 bits)
- Key intrinsics: `_mm512_load_pd`, `_mm512_store_pd`, `_mm512_fmadd_pd`, `_mm512_div_pd`
- Masked operations: `_mm512_mask_*` for branch-free conditional updates
- Gather/scatter: `_mm512_i64gather_pd` for indirect loads

### Approach: Pre-gather neighbor data

Instead of: Loop → indirect access → compute → accumulate

Do: Gather neighbors into contiguous buffer → vectorized loop → accumulate

```cpp
// Pre-gather phase: collect neighbor data into aligned buffers
alignas(64) double nb_pos_x[BATCH_SIZE];
alignas(64) double nb_pos_y[BATCH_SIZE];
alignas(64) double nb_pos_z[BATCH_SIZE];
alignas(64) double nb_vel_x[BATCH_SIZE];
// ... etc

for (int i = 0; i < num_neighbors; i += BATCH_SIZE) {
    // Gather batch of neighbor data
    gather_neighbor_data(i, batch, nb_pos_x, ...);

    // Vectorized force calculation
    vectorized_force_accumulate(batch_count, pos, vel, nb_pos_x, ...);
}
```

### Force Calculation Vectorization

The force calculation per neighbor:
1. dx = nb_pos - pos (3 ops)
2. dv = nb_vel - vel (3 ops)
3. r2 = dx² + dy² + dz² (reduction)
4. vx = dx*dvx + dy*dvy + dz*dvz (reduction)
5. m_r3 = mass / (r2 * sqrt(r2)) (expensive)
6. a += m_r3 * dx (3 ops)
7. adot += m_r3 * (dv - 3*dx*vx/r2) (9 ops)

**Vectorization pattern:** Process 8 neighbors at once:
- Load 8 neighbor positions (gather)
- Compute 8 dx, dy, dz values
- Compute 8 r2 values
- Compute 8 m_r3 values (rsqrt approximation + refinement)
- Accumulate 8 force contributions

### rsqrt Optimization

Instead of `1/(r2*sqrt(r2))` = `r^(-3)`:
```cpp
__m512d r2_vec = ...;
__m512d r_inv = _mm512_rsqrt14_pd(r2_vec);  // 1/sqrt(r2) approximation
// Refine with Newton-Raphson if needed
__m512d r3_inv = _mm512_mul_pd(r_inv, _mm512_mul_pd(r_inv, r_inv));  // r^(-3)
__m512d m_r3 = _mm512_mul_pd(mass_vec, r3_inv);
```

---

## 6. Implementation Plan Options

### Option A: Minimal Change - Vectorize Inner Loop Only

**Scope:** Add AVX-512 intrinsics to force accumulation, keep neighbor access pattern.

**Pros:** Smaller code change, lower risk.
**Cons:** Still limited by memory access patterns.

**Expected improvement:** 10-20%

### Option B: Pre-gather + Vectorize (Recommended)

**Scope:**
1. Pre-gather active neighbor data into aligned buffers
2. Use ParticleData SoA for faster neighbor data access
3. Vectorize force calculation with AVX-512

**Pros:** Best SIMD utilization, better cache usage.
**Cons:** Larger code change, needs careful memory management.

**Expected improvement:** 30-50%

### Option C: Full SoA Path

**Scope:** Port entire irregular force path to use ParticleData SoA consistently.

**Pros:** Maximum performance potential.
**Cons:** Largest change, risk of introducing bugs.

**Expected improvement:** 40-60%

---

## 7. Recommended Implementation

Based on context decisions (aggressive restructuring permitted, explicit AVX-512 intrinsics):

**Implement Option B with elements of C:**

1. **Create neighbor data buffer** - Pre-gather position, velocity, mass for active neighbors
2. **Use SoA layout for buffer** - Aligned arrays for SIMD
3. **Implement AVX-512 force kernel** - Process 8 neighbors per iteration
4. **Handle remainder** - Masked operations for non-multiple-of-8 neighbor counts
5. **Scalar fallback** - `#ifdef __AVX512F__` guard with scalar fallback

### File Changes Required

| File | Change |
|------|--------|
| `src/Particle/compute_acceleration.cpp` | Add vectorized neighbor loop |
| `src/particle_data.h` (optional) | Add SoA prediction helper |
| `src/Makefile` | Ensure `-mavx512f` flag (already via `-march=native`) |

### Memory Alignment

AVX-512 requires 64-byte alignment for best performance:
```cpp
alignas(64) double buffer[8];  // C++11 alignas
```

---

## 8. Validation Strategy

1. **Correctness:** Compare total energy before/after optimization
2. **Tolerance:** dE/E0 < 1e-4 (same as baseline)
3. **Profiling:** Before/after timing of IrregularForce
4. **Target:** At least 20% reduction in irregular force time

---

## 9. Risk Assessment

| Risk | Mitigation |
|------|------------|
| Floating-point differences | Accept small differences due to SIMD operation reordering |
| Memory alignment issues | Use alignas(64) for all SIMD buffers |
| Neighbor count < 8 | Use masked operations for remainder |
| Active particle filtering | Pre-filter during gather phase |
| Build compatibility | Use `#ifdef __AVX512F__` for fallback |

---

## 10. Reference Code Patterns

### Pre-gather Pattern

```cpp
void gather_active_neighbors(const Particle* ptcls, const int* neighbor_indices,
                              int num_neighbors, int offset, double dt,
                              double* __restrict pos_x, double* __restrict pos_y, double* __restrict pos_z,
                              double* __restrict vel_x, double* __restrict vel_y, double* __restrict vel_z,
                              double* __restrict mass, int& active_count) {
    active_count = 0;
    for (int i = 0; i < num_neighbors; i++) {
        const Particle* nb = &ptcls[neighbor_indices[offset + i]];
        if (!nb->is_active) continue;

        // Predict and store
        double p[3], v[3];
        nb->predict_particle_second_order(dt, p, v);
        pos_x[active_count] = p[0];
        pos_y[active_count] = p[1];
        pos_z[active_count] = p[2];
        vel_x[active_count] = v[0];
        vel_y[active_count] = v[1];
        vel_z[active_count] = v[2];
        mass[active_count] = nb->mass;
        active_count++;
    }
}
```

### AVX-512 Force Kernel

```cpp
#ifdef __AVX512F__
void vectorized_force_accumulate(int count,
                                  const double pos[3], const double vel[3],
                                  const double* nb_pos_x, const double* nb_pos_y, const double* nb_pos_z,
                                  const double* nb_vel_x, const double* nb_vel_y, const double* nb_vel_z,
                                  const double* nb_mass,
                                  double a_out[3], double adot_out[3]) {
    __m512d ax = _mm512_setzero_pd();
    __m512d ay = _mm512_setzero_pd();
    __m512d az = _mm512_setzero_pd();
    __m512d adx = _mm512_setzero_pd();
    __m512d ady = _mm512_setzero_pd();
    __m512d adz = _mm512_setzero_pd();

    __m512d px = _mm512_set1_pd(pos[0]);
    __m512d py = _mm512_set1_pd(pos[1]);
    __m512d pz = _mm512_set1_pd(pos[2]);
    __m512d vx = _mm512_set1_pd(vel[0]);
    __m512d vy = _mm512_set1_pd(vel[1]);
    __m512d vz = _mm512_set1_pd(vel[2]);
    __m512d three = _mm512_set1_pd(3.0);

    int i = 0;
    for (; i + 8 <= count; i += 8) {
        // Load 8 neighbors
        __m512d nb_px = _mm512_load_pd(nb_pos_x + i);
        __m512d nb_py = _mm512_load_pd(nb_pos_y + i);
        __m512d nb_pz = _mm512_load_pd(nb_pos_z + i);
        __m512d nb_vx = _mm512_load_pd(nb_vel_x + i);
        __m512d nb_vy = _mm512_load_pd(nb_vel_y + i);
        __m512d nb_vz = _mm512_load_pd(nb_vel_z + i);
        __m512d nb_m  = _mm512_load_pd(nb_mass + i);

        // dx = nb_pos - pos
        __m512d dx = _mm512_sub_pd(nb_px, px);
        __m512d dy = _mm512_sub_pd(nb_py, py);
        __m512d dz = _mm512_sub_pd(nb_pz, pz);

        // dv = nb_vel - vel
        __m512d dvx = _mm512_sub_pd(nb_vx, vx);
        __m512d dvy = _mm512_sub_pd(nb_vy, vy);
        __m512d dvz = _mm512_sub_pd(nb_vz, vz);

        // r2 = dx*dx + dy*dy + dz*dz
        __m512d r2 = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));

        // vdotr = dx*dvx + dy*dvy + dz*dvz
        __m512d vdotr = _mm512_fmadd_pd(dx, dvx, _mm512_fmadd_pd(dy, dvy, _mm512_mul_pd(dz, dvz)));

        // r_inv = 1/sqrt(r2), r3_inv = r_inv^3
        __m512d r_inv = _mm512_rsqrt14_pd(r2);
        // Newton-Raphson refinement for double precision
        __m512d half = _mm512_set1_pd(0.5);
        __m512d three_halves = _mm512_set1_pd(1.5);
        __m512d r_inv_sq = _mm512_mul_pd(r_inv, r_inv);
        r_inv = _mm512_mul_pd(r_inv, _mm512_fnmadd_pd(half, _mm512_mul_pd(r2, r_inv_sq), three_halves));

        __m512d r3_inv = _mm512_mul_pd(r_inv, _mm512_mul_pd(r_inv, r_inv));
        __m512d m_r3 = _mm512_mul_pd(nb_m, r3_inv);

        // a += m_r3 * dx
        ax = _mm512_fmadd_pd(m_r3, dx, ax);
        ay = _mm512_fmadd_pd(m_r3, dy, ay);
        az = _mm512_fmadd_pd(m_r3, dz, az);

        // coeff = -3 * vdotr / r2
        __m512d coeff = _mm512_mul_pd(_mm512_div_pd(vdotr, r2), three);

        // adot += m_r3 * (dv - 3*dx*vdotr/r2)
        adx = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dx, dvx), adx);
        ady = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dy, dvy), ady);
        adz = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dz, dvz), adz);
    }

    // Reduce vectors to scalars
    a_out[0] = _mm512_reduce_add_pd(ax);
    a_out[1] = _mm512_reduce_add_pd(ay);
    a_out[2] = _mm512_reduce_add_pd(az);
    adot_out[0] = _mm512_reduce_add_pd(adx);
    adot_out[1] = _mm512_reduce_add_pd(ady);
    adot_out[2] = _mm512_reduce_add_pd(adz);

    // Handle remainder with scalar code
    for (; i < count; i++) {
        // ... scalar fallback
    }
}
#endif
```

---

## RESEARCH COMPLETE

Research findings documented. Ready for planning phase.

**Key decisions for planner:**
1. Use pre-gather + vectorize approach (Option B)
2. Add AVX-512 kernel with scalar fallback
3. Target the neighbor loop as primary optimization
4. Also vectorize CM loop with same pattern
5. Verify energy conservation after each change

---
*Research completed: 2026-01-17*
