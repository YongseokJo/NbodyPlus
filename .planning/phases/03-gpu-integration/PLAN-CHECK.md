# Phase 3 Plan Verification

## Phase Goal

Convert GPU data structures and kernels to SoA layout.

## Requirements Coverage

| REQ-ID | Requirement | Plan | Status |
|--------|-------------|------|--------|
| GPU-01 | Create `ParticleDataGPU` device container | Plan 01, 02 | ✓ Covered |
| GPU-02 | Implement host-to-device transfer | Plan 02, 04 | ✓ Covered |
| GPU-03 | Implement device-to-host transfer | Plan 02, 04 | ✓ Covered |
| GPU-04 | Update `compute_forces_kernel` for SoA | Plan 03 | ✓ Covered |
| GPU-05 | Update `predict_particle_kernel` for SoA | N/A | ⚠ Not applicable |

**Note on GPU-05**: The current codebase does not have a GPU prediction kernel. Prediction is done on CPU. This requirement may be obsolete or for a future feature.

## Plan Dependencies

```
Plan 01 (Header)
    ↓
Plan 02 (Implementation)
    ↓
Plan 03 (Kernel update)
    ↓
Plan 04 (Host integration)
    ↓
Plan 05 (Build test)
```

All plans can be executed sequentially (Wave 1).

## Pitfalls Addressed

| Pitfall | Plan | Mitigation |
|---------|------|------------|
| #4 GPU Transfer Fragmentation | 02, 04 | CUDA streams for async overlap |
| #8 Mixed Precision | 02 | Explicit float/double conversion paths |

## Must-Have Checklist

- [ ] `ParticleDataGPU` struct exists with device pointers
- [ ] `allocate()` creates all device arrays
- [ ] `deallocate()` frees all device arrays
- [ ] `copy_*_to_device()` transfers from host
- [ ] `copy_results_to_host()` transfers results back
- [ ] `compute_forces_kernel` accepts SoA parameters
- [ ] Shared memory uses separate arrays (not struct)
- [ ] Kernel invocation updated in `cuda_acceleration.cu`
- [ ] Build compiles without errors

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Shared memory usage exceeds limits | Low | High | Calculate: 8 arrays × 256 threads × 8 bytes = 16KB (under 48KB limit) |
| Interface change breaks callers | Medium | Medium | Provide legacy wrapper during transition |
| Performance regression | Low | Medium | Coalesced access should improve; validate with benchmark |

## Estimated Complexity

| Plan | Complexity | Notes |
|------|------------|-------|
| 01 | Low | Header file only |
| 02 | Medium | ~150 lines of CUDA code |
| 03 | Medium | Kernel rewrite, careful attention needed |
| 04 | High | Integration with existing code, interface changes |
| 05 | Low | Build system update |

## Verdict

**Plans are APPROVED for execution.**

All GPU requirements except GPU-05 are covered. GPU-05 appears to be for a feature not present in the current codebase (GPU prediction kernel).

---
*Generated: 2026-01-17*
