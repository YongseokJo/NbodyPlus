# Plan 04-04: Update Force Calculation Routines

## Frontmatter

```yaml
wave: 3
depends_on: [01-PLAN.md, 03-PLAN.md]
files_modified:
  - src/Particle/compute_acceleration.cpp
autonomous: true
```

## Objective

Update the three main force calculation functions to use SoA access patterns in their inner loops:
- `compute_acceleration_irr()` — Irregular force calculation
- `compute_acceleration_reg()` — Regular force calculation (CPU path)
- `update_regular_particle_cuda()` — Post-GPU neighbor transition

## Context

These functions are O(N×Nneighbor) or O(N²) and represent the hot path of the simulation. Converting inner loops to SoA access provides cache efficiency benefits.

Current pattern:
```cpp
ptcl = &particles[neighbors[offset + i]];
// Access ptcl->position, ptcl->velocity, ptcl->mass
```

Target pattern:
```cpp
int j = neighbors[offset + i];
// Access particle_data.get_pos_x(j), etc. or use pointers for batch access
```

## Strategy

Use a **bridge approach** similar to Phase 3:
1. Keep function signatures (Particle:: methods)
2. At function entry, prepare SoA views/caches
3. Inner loops use SoA access
4. At function exit, write results back to Particle struct

## Tasks

<task id="1">
Create acceleration accumulator struct that uses SoA layout for temporary results:
- a_tmp[3], adot_tmp[3] as contiguous arrays for vectorization
</task>

<task id="2">
Update compute_acceleration_irr() inner neighbor loop:
- Cache neighbor data in local SoA arrays before inner loop
- Use ParticleData accessors for neighbor field access
- Maintain identical mathematical operations
</task>

<task id="3">
Update compute_acceleration_reg() all-particle loop:
- Same pattern as irregular
- Handle neighbor_radius_sq comparison using SoA accessor
</task>

<task id="4">
Update update_regular_particle_cuda() neighbor transition loops:
- Use SoA access for hashTable operations
</task>

<task id="5">
Verify energy conservation by comparing acc_total values before/after refactoring
</task>

## Verification

- [ ] Force calculations produce bit-identical acceleration values
- [ ] Neighbor lists are correctly updated
- [ ] No change to force calculation physics

## must_haves

- acc_total, acc_regular, acc_irregular values unchanged
- Neighbor detection logic unchanged
- Energy conservation matches baseline
- Performance does not regress
