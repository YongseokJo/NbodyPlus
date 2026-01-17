# Plan 05: Build Integration and Compilation Test

## Goal

Integrate all new GPU files into the build system and verify compilation succeeds.

## Requirement

All GPU requirements (GPU-01 through GPU-04) — build verification

## Deliverable

Updated `src/Makefile` with GPU SoA files, successful compilation

## Tasks

### Task 1: Update Makefile for particle_data_gpu.cu

Add `particle_data_gpu.cu` to the CUDA sources.

**File**: `src/Makefile`

Find the CUDA source list (likely `CUDA_SRCS` or similar) and add:
```makefile
CUDA_SRCS += particle_data_gpu.cu
```

Or if there's a pattern like:
```makefile
CUDA_OBJS = cuda/cuda_kernels.o cuda/cuda_acceleration.o cuda/cuda_routines.o
```

Add:
```makefile
CUDA_OBJS = cuda/cuda_kernels.o cuda/cuda_acceleration.o cuda/cuda_routines.o particle_data_gpu.o
```

### Task 2: Add header dependency

Ensure `particle_data_gpu.h` is in the include path and dependencies are tracked.

```makefile
particle_data_gpu.o: particle_data_gpu.cu particle_data_gpu.h cuda/cuda_defs.h def.h
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
```

### Task 3: Compile individual files

Test that each file compiles:

```bash
cd src
make cuda/cuda_kernels.o
make cuda/cuda_acceleration.o
make particle_data_gpu.o
```

### Task 4: Full build test

Run full build to verify no link errors:

```bash
cd src
make clean
make
```

### Task 5: Verify object file sizes

Check that compiled object files are reasonable sizes:

```bash
ls -la cuda/*.o particle_data_gpu.o
```

## Verification

- [ ] `particle_data_gpu.o` compiles successfully
- [ ] `cuda_kernels.o` compiles with updated kernel signature
- [ ] `cuda_acceleration.o` compiles with SoA integration
- [ ] Full build succeeds with no link errors
- [ ] No new compiler warnings (or document expected ones)

## Dependencies

- Plans 01-04 (all GPU source changes)

## Troubleshooting

### Common issues:

1. **Missing include path**: Add `-I..` or `-I.` to NVCCFLAGS
2. **Undefined BATCH_SIZE**: Ensure `cuda_defs.h` is included
3. **Link error for ParticleData**: Need to link `particle_data.o` with CUDA objects
4. **Circular dependency**: May need forward declarations

### If kernel signature change breaks callers:

Temporarily keep both old and new kernel versions:
```cpp
// Legacy wrapper (remove in Phase 4)
__global__ void compute_forces_legacy(
    const i_particle_t* d_Ip, const j_particle_t* d_Jp, ...);

// New SoA version
__global__ void compute_forces(
    const cuda_real_t* d_i_pos_x, ...);
```

---
*Generated: 2026-01-17*
