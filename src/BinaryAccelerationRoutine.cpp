#include <iostream>
#include "global.h"

void BinaryAccelerationRoutine(double next_time) {

    for (Binary* ptclBin: BinaryList) {
        ptclBin->IntegrateBinary(next_time);
    }

}