#ifndef PERFORMANCE_LOG_H
#define PERFORMANCE_LOG_H

#include <string>
#include "performance.h"

void WritePerformanceLogCsv(const std::string& directory_path,
                            const char* output_prefix,
                            int output_num,
                            double sim_time_myr,
                            const Performance& perf,
                            int num_particles);

#endif
