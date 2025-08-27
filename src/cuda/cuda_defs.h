#ifndef CUDA_DEFS_H
#define CUDA_DEFS_H

#define _PROFILE
#define _six 6
#define _two 2
#define _seven 7
#define new_size(A) ((A > 1024) ? int(pow(2,ceil(log(A)/log(2.0)))) : 1024)

#include "../def.h"

struct Iparticle {
    CUDA_REAL posx;
    CUDA_REAL posy;
    CUDA_REAL posz;
    CUDA_REAL r2;

    CUDA_REAL velx;
    CUDA_REAL vely;
    CUDA_REAL velz;
    CUDA_REAL dtr;
};

struct Jparticle {
    CUDA_REAL posx;
    CUDA_REAL posy;
    CUDA_REAL posz;
    CUDA_REAL mass;

    CUDA_REAL velx;
    CUDA_REAL vely;
    CUDA_REAL velz;
    int index;
#ifndef CUDA_FLOAT
    int pad;
#endif
};

#endif
