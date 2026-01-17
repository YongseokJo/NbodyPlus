# Phase 5: SDAR Compatibility — Research

**Researcher:** gsd-phase-researcher (inline)
**Date:** 2026-01-17

## RESEARCH COMPLETE

---

## Executive Summary

SDAR integration in ABYSS operates on `Particle` objects through the `AR::TimeTransformedSymplecticIntegrator` template system. The library accesses particles by direct member access (`p.mass`, `p.position[dim]`, `p.velocity[dim]`) and uses macro aliases to map SDAR naming conventions.

**Key Finding:** SDAR does NOT directly access the main `particles[]` array during integration. Instead, it operates on **copies** stored in `sym_int.particles` — a separate `COMM::ParticleGroup<Particle>`. This means the proxy pattern is not needed at the core SDAR level. Instead, we need sync points before/after SDAR operations.

---

## 1. SDAR Library Architecture

### 1.1 Template Instantiation
```cpp
// group.h:41
AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
```

The SDAR integrator is templated on the `Particle` type. It expects:
- `Particle.mass` — double
- `Particle.position[3]` — double array
- `Particle.velocity[3]` — double array
- `Particle.pid` — int

### 1.2 Macro Aliases (group.h:13-29)
```cpp
#define Mass mass
#define Position position
#define Velocity velocity
#define PID pid
```
These map SDAR's expected capitalized names to ABYSS's lowercase names.

---

## 2. SDAR Entry Points

### 2.1 Where SDAR Particle Copies Are Created

**fb_initialization.cpp:63-76** — `initialIntegrator()`:
```cpp
sym_int.particles.addMemberAndAddress(*members);
```
Here, SDAR receives **copies** of Particle objects, not references.

### 2.2 Where SDAR Integration Runs

**fb_integration.cpp:20-233** — `Group::ARIntegration()`:
```cpp
// Before integration: copy CM data from main array
sym_int.particles.cm.position[dim] = groupCM->position[dim];
sym_int.particles.cm.velocity[dim] = groupCM->velocity[dim];
sym_int.particles.cm.acc_total[dim][j] = groupCM->acc_total[dim][j];

// Integration (SDAR operates on its internal copies)
bin_interrupt = sym_int.integrateToTime(next_time*enzo_time_step);

// After integration: copy results back to main array
particles[members->particle_index].position[dim] = groupCM->position[dim] + members->position[dim];
particles[members->particle_index].velocity[dim] = groupCM->velocity[dim] + members->velocity[dim];
particles[members->particle_index].mass = members->mass;
```

### 2.3 Perturber Access Pattern

**ar_interaction.hpp:218** — Perturbers access main `particles[]` array:
```cpp
pertj = &particles[pert_adr[j]];
```

This is the **only** place SDAR code directly accesses the global `particles[]` array.

---

## 3. Fields Used by SDAR

### 3.1 Directly Accessed in ar_interaction.hpp

| Field | Usage |
|-------|-------|
| `mass` | Force calculations |
| `position[3]` | Distance, acceleration |
| `velocity[3]` | Jerk calculations, prediction |
| `acc_total[3][2]` | Orders 0,1 for prediction (not full [3][4]) |
| `current_time_irr` | Time delta calculation |
| `is_active` | Skip inactive particles |
| `cm_particle_index` | Identify CM particles |
| `pid` | Debug output |
| `num_neighbors` | Neighbor iteration |
| `neighbors_offset` | Neighbor array access |
| `particle_index` | Array indexing |
| `particle_type` | Merger radius logic |
| `radius` | Collision detection |
| `spin_param[3]` | Binary merger physics |
| `binary_state` | Interrupt state |
| `time_check` | Collision candidate timing |

### 3.2 Fields Modified by SDAR

In `ARIntegration()` and related functions:
- `position[3]` — updated from integration
- `velocity[3]` — updated from integration
- `mass` — may change in mergers
- `binary_state` — interrupt state changes
- `current_time_irr` — time advancement
- `is_active` — deactivation on merger
- `spin_param[3]` — remnant spin calculation

---

## 4. Current Sync Pattern Analysis

### 4.1 Before SDAR Integration
```cpp
// Copy from main Particle to SDAR's internal copy
sym_int.particles.cm.position[dim] = groupCM->position[dim];
```

### 4.2 After SDAR Integration
```cpp
// Copy results back to main particles[] array
particles[members->particle_index].position[dim] = ...;
```

**This pattern already exists.** The question is: with SoA, do we need to sync the main `particles[]` array first?

---

## 5. Impact Assessment

### 5.1 Main Particles Array Access

The global `extern Particle *particles` is accessed in two ways:
1. **Group members**: Through `sym_int.particles[]` copies (already synced)
2. **Perturbers**: Direct access via `particles[pert_adr[j]]`

For perturber access, the prediction code reads:
```cpp
pertj->predict_particle_second_order(dt, xp[n_pert_active], vp[...]);
```

This uses `acc_total[dim][0]` and `acc_total[dim][1]` — the first two orders only.

### 5.2 What Needs SoA Sync

Before SDAR operations:
- Perturber particles need `position, velocity, mass, acc_total[dim][0:1]` from SoA

After SDAR operations:
- Modified group member particles need their changes written back to SoA

---

## 6. Recommended Approach

### 6.1 No Proxy Class Needed

SDAR already operates on its own copies. The issue is ensuring `particles[]` is in sync with `ParticleData` (SoA).

### 6.2 Sync Points Needed

1. **Before `NewFBInitialization()`**: Sync group member and potential perturber particles from SoA
2. **Before `Group::ARIntegration()`**: Sync perturber particles from SoA (group members already copied)
3. **After `FBTermination()`**: Sync modified member particles back to SoA

### 6.3 Implementation Strategy

**Option A: Sync-on-entry/exit pattern**
- Add sync calls at FewBody entry/exit points
- Minimal code changes
- May sync more particles than needed

**Option B: Lazy sync with dirty tracking**
- Track which particles are dirty in SoA
- Only sync dirty particles
- More complex but more efficient

**Recommendation:** Option A for Phase 5 (simpler), with Option B as future optimization.

---

## 7. Files to Modify

| File | Changes |
|------|---------|
| `src/FewBody/fb_initialization.cpp` | Add SoA→Particle sync before initialization |
| `src/FewBody/fb_integration.cpp` | Add SoA→Particle sync for perturbers before integration |
| `src/FewBody/fb_termination.cpp` | Add Particle→SoA sync after termination |
| `src/particle_data.h` | (Already has sync methods) |
| `src/particle_data.cpp` | Verify sync methods handle all SDAR fields |

---

## 8. Risk Assessment

### 8.1 Low Risk
- SDAR code doesn't need modification
- Sync methods already exist in ParticleData
- Well-defined entry/exit points

### 8.2 Moderate Risk
- Performance: extra sync overhead
- Correctness: must sync all required fields

### 8.3 Mitigation
- Profile sync overhead after implementation
- Add assertions to verify sync correctness
- Test with existing SDAR test cases

---

## 9. Testing Strategy

1. **Unit test**: Verify sync_from_particle/sync_to_particle handle all SDAR fields
2. **Integration test**: Run binary formation/termination scenarios
3. **Regression test**: Compare energy conservation with/without SoA

---

*Research complete. Proceeding to planning.*
