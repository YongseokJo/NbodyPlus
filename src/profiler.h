#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <climits>

// Timer IDs for fast lookup (avoid string hashing in hot paths)
enum class TimerID : int {
    // Main loop timers
    WholeRoutine = 0,

    // Irregular step timers
    IrregularForce,
    IrregularUpdate,
    IrregularTotal,

    // Regular step timers
    RegularForce,
    RegularUpdate,
    RegularTotal,
    RegularSendToGPU,
    RegularGPU,
    RegularAdjust,

    // Few-body timers
    FewBodyTermination,
    FewBodySearch,
    FewBodyInitialization,
    FewBodyIntegration,  // SDAR integration time

    // Data structure timers
    SkipListCreate,
    SkipListUpdate,
    RegularMap,
    UpdateNextRegTime,

    // MPI communication timers
    MPISend,
    MPIRecv,
    MPIWait,
    MPIBarrier,
    MPIReduce,
    MPIBcast,
    MPIWindowSync,

    // Worker timers
    WorkerCompute,
    WorkerIdle,
    QueueWait,

    // I/O timers
    FileWrite,
    FileRead,

    // Stellar evolution timers
    StellarEvolution,
    BinaryStellarEvolution,

    // Memory operations
    NeighborListUpdate,
    ParticleUpdate,

    // Count for array sizing
    NUM_TIMERS
};

// Statistics for a single timer
struct TimerStats {
    long long total_ns = 0;        // Total time in nanoseconds
    long long count = 0;           // Number of calls
    long long min_ns = LLONG_MAX;  // Minimum single call time
    long long max_ns = 0;          // Maximum single call time
    long long last_ns = 0;         // Last recorded duration

    // For per-output-interval tracking
    long long interval_total_ns = 0;
    long long interval_count = 0;
    long long interval_min_ns = LLONG_MAX;
    long long interval_max_ns = 0;

    void record(long long duration_ns) {
        total_ns += duration_ns;
        count++;
        last_ns = duration_ns;
        min_ns = std::min(min_ns, duration_ns);
        max_ns = std::max(max_ns, duration_ns);

        interval_total_ns += duration_ns;
        interval_count++;
        interval_min_ns = std::min(interval_min_ns, duration_ns);
        interval_max_ns = std::max(interval_max_ns, duration_ns);
    }

    void resetInterval() {
        interval_total_ns = 0;
        interval_count = 0;
        interval_min_ns = LLONG_MAX;
        interval_max_ns = 0;
    }

    double totalSeconds() const { return total_ns * 1e-9; }
    double intervalSeconds() const { return interval_total_ns * 1e-9; }
    double meanNs() const { return count > 0 ? static_cast<double>(total_ns) / count : 0.0; }
    double intervalMeanNs() const {
        return interval_count > 0 ? static_cast<double>(interval_total_ns) / interval_count : 0.0;
    }
    double meanMicros() const { return meanNs() * 1e-3; }
    double intervalMeanMicros() const { return intervalMeanNs() * 1e-3; }
};

// Main profiler class
class Profiler {
public:
    static Profiler& instance() {
        static Profiler inst;
        return inst;
    }

    // Start a timer
    void start(TimerID id) {
        start_times_[static_cast<int>(id)] = std::chrono::high_resolution_clock::now();
    }

    // Stop a timer and record the duration
    void stop(TimerID id) {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - start_times_[static_cast<int>(id)]).count();
        stats_[static_cast<int>(id)].record(duration);
    }

    // Get statistics for a timer
    const TimerStats& getStats(TimerID id) const {
        return stats_[static_cast<int>(id)];
    }

    TimerStats& getStats(TimerID id) {
        return stats_[static_cast<int>(id)];
    }

    // Increment a counter without timing (for tracking operation counts)
    void incrementCount(TimerID id, long long amount = 1) {
        stats_[static_cast<int>(id)].count += amount;
        stats_[static_cast<int>(id)].interval_count += amount;
    }

    // Reset interval statistics (called at each output)
    void resetIntervalStats() {
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            stats_[i].resetInterval();
        }
    }

    // Reset all statistics
    void resetAll() {
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            stats_[i] = TimerStats();
        }
    }

    // Get timer name for output
    static const char* getTimerName(TimerID id) {
        static const char* names[] = {
            "WholeRoutine",
            "IrregularForce",
            "IrregularUpdate",
            "IrregularTotal",
            "RegularForce",
            "RegularUpdate",
            "RegularTotal",
            "RegularSendToGPU",
            "RegularGPU",
            "RegularAdjust",
            "FewBodyTermination",
            "FewBodySearch",
            "FewBodyInitialization",
            "FewBodyIntegration",
            "SkipListCreate",
            "SkipListUpdate",
            "RegularMap",
            "UpdateNextRegTime",
            "MPISend",
            "MPIRecv",
            "MPIWait",
            "MPIBarrier",
            "MPIReduce",
            "MPIBcast",
            "MPIWindowSync",
            "WorkerCompute",
            "WorkerIdle",
            "QueueWait",
            "FileWrite",
            "FileRead",
            "StellarEvolution",
            "BinaryStellarEvolution",
            "NeighborListUpdate",
            "ParticleUpdate"
        };
        return names[static_cast<int>(id)];
    }

    // Print summary to stdout (interval stats)
    void printIntervalSummary(std::ostream& os, double sim_time_myr) const {
        const auto& whole = stats_[static_cast<int>(TimerID::WholeRoutine)];
        if (whole.interval_total_ns == 0) return;

        os << "\n==================== Performance Summary ====================\n";
        os << "Simulation Time: " << std::fixed << std::setprecision(6) << sim_time_myr << " Myr\n";
        os << "Wall-clock Time: " << std::fixed << std::setprecision(2)
           << whole.intervalSeconds() << " s\n";
        os << "=============================================================\n\n";

        os << std::left << std::setw(28) << "Timer"
           << std::right << std::setw(12) << "Time (s)"
           << std::setw(10) << "Percent"
           << std::setw(12) << "Calls"
           << std::setw(14) << "Avg (us)"
           << std::setw(14) << "Max (us)" << "\n";
        os << std::string(90, '-') << "\n";

        // Print in logical groups
        printTimerLine(os, TimerID::IrregularForce, whole.interval_total_ns);
        printTimerLine(os, TimerID::IrregularUpdate, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::RegularForce, whole.interval_total_ns);
        printTimerLine(os, TimerID::RegularSendToGPU, whole.interval_total_ns);
        printTimerLine(os, TimerID::RegularGPU, whole.interval_total_ns);
        printTimerLine(os, TimerID::RegularAdjust, whole.interval_total_ns);
        printTimerLine(os, TimerID::RegularUpdate, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::FewBodyTermination, whole.interval_total_ns);
        printTimerLine(os, TimerID::FewBodyInitialization, whole.interval_total_ns);
        printTimerLine(os, TimerID::FewBodyIntegration, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::SkipListCreate, whole.interval_total_ns);
        printTimerLine(os, TimerID::SkipListUpdate, whole.interval_total_ns);
        printTimerLine(os, TimerID::UpdateNextRegTime, whole.interval_total_ns);
        printTimerLine(os, TimerID::RegularMap, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::MPISend, whole.interval_total_ns);
        printTimerLine(os, TimerID::MPIRecv, whole.interval_total_ns);
        printTimerLine(os, TimerID::MPIWait, whole.interval_total_ns);
        printTimerLine(os, TimerID::MPIBarrier, whole.interval_total_ns);
        printTimerLine(os, TimerID::QueueWait, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::StellarEvolution, whole.interval_total_ns);
        printTimerLine(os, TimerID::BinaryStellarEvolution, whole.interval_total_ns);
        os << "\n";

        printTimerLine(os, TimerID::FileWrite, whole.interval_total_ns);

        os << std::string(90, '=') << "\n";

        // Calculate and print unaccounted time
        long long accounted = 0;
        for (int i = 1; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            accounted += stats_[i].interval_total_ns;
        }
        double unaccounted_pct = 100.0 * (whole.interval_total_ns - accounted) / whole.interval_total_ns;
        if (unaccounted_pct > 1.0) {
            os << "Note: " << std::fixed << std::setprecision(1) << unaccounted_pct
               << "% of time unaccounted (overhead, uninstrumented code)\n";
        }
    }

    // Write detailed stats to CSV file for analysis
    void writeCSV(const std::string& filename, double sim_time, int step) const {
        std::ofstream file(filename, std::ios::app);
        if (!file.is_open()) return;

        // Write header if file is empty
        file.seekp(0, std::ios::end);
        if (file.tellp() == 0) {
            file << "step,sim_time_myr";
            for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
                file << "," << getTimerName(static_cast<TimerID>(i)) << "_ns"
                     << "," << getTimerName(static_cast<TimerID>(i)) << "_count";
            }
            file << "\n";
        }

        file << step << "," << sim_time;
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            file << "," << stats_[i].interval_total_ns
                 << "," << stats_[i].interval_count;
        }
        file << "\n";
    }

    // Write JSON format for programmatic analysis
    void writeJSON(const std::string& filename, double sim_time, int step) const {
        std::ofstream file(filename);
        if (!file.is_open()) return;

        file << "{\n";
        file << "  \"step\": " << step << ",\n";
        file << "  \"sim_time_myr\": " << sim_time << ",\n";
        file << "  \"timers\": {\n";

        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            const auto& s = stats_[i];
            file << "    \"" << getTimerName(static_cast<TimerID>(i)) << "\": {\n";
            file << "      \"total_ns\": " << s.total_ns << ",\n";
            file << "      \"interval_ns\": " << s.interval_total_ns << ",\n";
            file << "      \"count\": " << s.count << ",\n";
            file << "      \"interval_count\": " << s.interval_count << ",\n";
            file << "      \"min_ns\": " << (s.min_ns == LLONG_MAX ? 0 : s.min_ns) << ",\n";
            file << "      \"max_ns\": " << s.max_ns << ",\n";
            file << "      \"mean_ns\": " << s.meanNs() << "\n";
            file << "    }";
            if (i < static_cast<int>(TimerID::NUM_TIMERS) - 1) file << ",";
            file << "\n";
        }

        file << "  }\n";
        file << "}\n";
    }

private:
    Profiler() {
        stats_.resize(static_cast<int>(TimerID::NUM_TIMERS));
        start_times_.resize(static_cast<int>(TimerID::NUM_TIMERS));
    }

    void printTimerLine(std::ostream& os, TimerID id, long long total_ns) const {
        const auto& s = stats_[static_cast<int>(id)];
        if (s.interval_count == 0 && s.interval_total_ns == 0) return;

        double pct = total_ns > 0 ? 100.0 * s.interval_total_ns / total_ns : 0.0;
        double avg_us = s.intervalMeanMicros();
        double max_us = s.interval_max_ns * 1e-3;

        os << std::left << std::setw(28) << getTimerName(id)
           << std::right << std::fixed << std::setprecision(3) << std::setw(12) << s.intervalSeconds()
           << std::setprecision(1) << std::setw(9) << pct << "%"
           << std::setw(12) << s.interval_count
           << std::setprecision(2) << std::setw(14) << avg_us
           << std::setw(14) << max_us << "\n";
    }

    std::vector<TimerStats> stats_;
    std::vector<std::chrono::high_resolution_clock::time_point> start_times_;
};

// RAII-style scoped timer for automatic start/stop
class ScopedTimer {
public:
    explicit ScopedTimer(TimerID id) : id_(id) {
        Profiler::instance().start(id_);
    }

    ~ScopedTimer() {
        Profiler::instance().stop(id_);
    }

    // Disable copy/move
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
    TimerID id_;
};

// Convenience macros for profiling
#ifdef PERFORMANCETRACE

#define PROFILE_SCOPE(timer_id) ScopedTimer _scoped_timer_##__LINE__(timer_id)
#define PROFILE_START(timer_id) Profiler::instance().start(timer_id)
#define PROFILE_STOP(timer_id) Profiler::instance().stop(timer_id)
#define PROFILE_COUNT(timer_id) Profiler::instance().incrementCount(timer_id)
#define PROFILE_COUNT_N(timer_id, n) Profiler::instance().incrementCount(timer_id, n)

#else

#define PROFILE_SCOPE(timer_id) ((void)0)
#define PROFILE_START(timer_id) ((void)0)
#define PROFILE_STOP(timer_id) ((void)0)
#define PROFILE_COUNT(timer_id) ((void)0)
#define PROFILE_COUNT_N(timer_id, n) ((void)0)

#endif

// Global accessor for convenience
inline Profiler& profiler() {
    return Profiler::instance();
}

#endif // PROFILER_H
