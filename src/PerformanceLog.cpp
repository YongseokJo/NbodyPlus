#include "PerformanceLog.h"

#include <fstream>
#include <iomanip>

static bool fileExists(const std::string& path) {
    std::ifstream input(path);
    return input.good();
}

static std::string safePrefix(const char* prefix) {
    if (prefix == nullptr || prefix[0] == '\0') {
        return "perf";
    }
    return std::string(prefix);
}

void WritePerformanceLogCsv(const std::string& directory_path,
                            const char* output_prefix,
                            int output_num,
                            double sim_time_myr,
                            const Performance& perf,
                            int num_particles) {
    const std::string filename = directory_path + "/" + safePrefix(output_prefix) + "_perf.csv";
    const bool has_file = fileExists(filename);

    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        return;
    }

    if (!has_file) {
        out << "output_num,sim_time_myr,num_particles,elapsed_s,"
            << "irregular_force_s,irregular_update_s,fewbody_termination_s,"
            << "fewbody_initialization_s,regular_send_gpu_s,regular_gpu_s,"
            << "regular_adjust_s,regular_update_s,regular_force_s,"
            << "skiplist_create_s,skiplist_update_s,"
#ifdef MULTIMAP
            << "regular_map_s,"
#else
            << "update_next_reg_time_s,"
#endif
#ifdef SEVN
            << "stellar_evolution_s,"
#ifdef SEVN_BINARY
            << "binary_stellar_evolution_s,"
#endif
#endif
            << "elapsed_ns\n";
    }

    const double ns_to_s = 1e-9;

    out << output_num << ','
        << std::fixed << std::setprecision(6) << sim_time_myr << ','
        << num_particles << ','
        << perf.WholeRoutine * ns_to_s << ','
        << perf.IrregularForce * ns_to_s << ','
        << perf.IrregularUpdate * ns_to_s << ','
        << perf.FewBodyTermination * ns_to_s << ','
        << perf.FewBodyInitialization * ns_to_s << ','
        << perf.RegularSendAllParticlesToGPU * ns_to_s << ','
        << perf.RegularGPU * ns_to_s << ','
        << perf.RegularAdjust * ns_to_s << ','
        << perf.RegularUpdate * ns_to_s << ','
        << perf.RegularForce * ns_to_s << ','
        << perf.SkipListCreate * ns_to_s << ','
        << perf.SkipListUpdate * ns_to_s << ','
#ifdef MULTIMAP
        << perf.RegularMap * ns_to_s << ','
#else
        << perf.UpdateNextRegTime * ns_to_s << ','
#endif
#ifdef SEVN
        << perf.StellarEvolution * ns_to_s << ','
#ifdef SEVN_BINARY
        << perf.BinaryStellarEvolution * ns_to_s << ','
#endif
#endif
        << perf.WholeRoutine
        << '\n';
}
