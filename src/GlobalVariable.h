#ifndef GLOBAL_VARIABLE_H
#define GLOBAL_VARIABLE_H

#include "def.h"

// ============================================================================
// Global state structure (shared across MPI processes)
// ============================================================================
struct GlobalVariable {
    int last_particle_index;    // Index of last active particle
    ull_t next_reg_time_block;  // Next regular time block
};

#endif
