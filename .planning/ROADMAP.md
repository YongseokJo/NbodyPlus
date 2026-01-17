# ABYSS SoA Conversion — Roadmap

## Milestone 1: AoS to SoA Conversion (v1.0)

### Phase 1: Core SoA Container

**Goal:** Create the `ParticleData` SoA container class with memory management.

**Requirements:** CONT-01, CONT-02, CONT-03

**Deliverables:**
- `src/particle_data.h` — SoA container declaration
- `src/particle_data.cpp` — Implementation with allocation/deallocation
- Accessor functions for all ~62 arrays (critical + important fields)
- Padding logic to avoid power-of-2 cache conflicts

**Pitfalls to address:**
- #1 Conversion overhead (design for permanent SoA storage)
- #2 Cache associativity (add padding to array sizes)

---

### Phase 2: MPI Integration

**Goal:** Convert MPI shared memory to use multiple `MPI_Win` objects.

**Requirements:** MPI-01, MPI-02, MPI-03

**Deliverables:**
- Update `src/mpi_routines.cpp` with `MPI_Win_allocate_shared` per array
- Window management (create, free) for all SoA arrays
- Verify cross-rank access works correctly

**Pitfalls to address:**
- #6 MPI window explosion (monitor creation time, test with target MPI)

---

### Phase 3: GPU Integration

**Goal:** Convert GPU data structures and kernels to SoA layout.

**Requirements:** GPU-01, GPU-02, GPU-03, GPU-04, GPU-05

**Deliverables:**
- `src/particle_data_gpu.h` — GPU SoA container declaration
- `src/particle_data_gpu.cu` — Device allocation and transfers
- Update `src/cuda/cuda_kernels.cu` — Kernel signatures and access patterns
- Update `src/cuda/cuda_acceleration.cu` — Transfer dispatch

**Pitfalls to address:**
- #4 GPU transfer fragmentation (use CUDA streams for overlap)
- #8 Mixed precision (be explicit about float vs double)

---

### Phase 4: CPU Routines

**Goal:** Update all CPU force calculation and integration routines.

**Requirements:** CPU-01, CPU-02, CPU-03, CPU-04, CPU-05

**Deliverables:**
- Update `src/Particle/compute_acceleration.cpp`
- Update `src/regular_routines.cpp`
- Update `src/irregular_routines.cpp`
- Update `src/timestep_routines.cpp`
- Update prediction/correction in `src/Particle/update_particle.cpp`

**Pitfalls to address:**
- #7 Energy conservation regression (validate after changes)

---

### Phase 5: SDAR Compatibility

**Goal:** Maintain SDAR library integration via proxy pattern.

**Requirements:** SDAR-01, SDAR-02, SDAR-03

**Deliverables:**
- `src/particle_proxy.h` — Proxy class for SDAR compatibility
- Update `src/FewBody/` files to use proxy
- Sync logic before/after SDAR operations

**Pitfalls to address:**
- #5 SDAR integration breakage (use proxy pattern, test thoroughly)

---

### Phase 6: I/O Updates

**Goal:** Update file I/O to work with SoA layout.

**Requirements:** IO-01, IO-02

**Deliverables:**
- Update `src/read_write.cpp` — HDF5 output from SoA
- Update `src/initialization_routines.cpp` — Input file reading

---

### Phase 7: Validation

**Goal:** Verify correctness and measure performance improvement.

**Requirements:** VAL-01, VAL-02, VAL-03, VAL-04

**Deliverables:**
- Run `workflow/bin/submit.sh --tag soa_final`
- Compare `summary_runs.tsv` metrics vs baseline:
  - `dE_over_E0_mean` ≤ 3.2e-5 (energy conservation)
  - `total_wall_s` improvement (performance)
- Document results

---

## Phase Summary

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 1 | Core SoA Container | CONT-01, CONT-02, CONT-03 | Complete ✓ |
| 2 | MPI Integration | MPI-01, MPI-02, MPI-03 | Complete ✓ |
| 3 | GPU Integration | GPU-01..05 | Complete ✓ |
| 4 | CPU Routines | CPU-01..05 | Pending |
| 5 | SDAR Compatibility | SDAR-01, SDAR-02, SDAR-03 | Pending |
| 6 | I/O Updates | IO-01, IO-02 | Pending |
| 7 | Validation | VAL-01..04 | Pending |

---

## Baseline Reference

From `summary_runs.tsv` (baseline_20260116_234947):
- `dE_over_E0_mean`: 3.20665e-05
- `total_wall_s`: 28.33s
- Hardware: Intel Xeon Gold 6338, NVIDIA A30

---
*Generated: 2026-01-17*
