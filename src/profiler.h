#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <climits>
#include <memory>
#include <sstream>
#ifdef USE_MPI
#include <mpi.h>
#endif

// Timer IDs for fast lookup (avoid string hashing in hot paths)
enum class TimerID : int {
    // Main loop timers
    WholeRoutine = 0,

    // Irregular step timers
    IrregularForce,
    IrregularUpdate,
    IrregularTotal,

    // Irregular force sub-timers (Phase 8)
    IrregularNeighborLoop,    // Main neighbor loop
    IrregularCMLoop,          // CM particle loop
    IrregularCorrection,      // 4th order correction
    IrregularPredict,         // Self and neighbor predictions
    IrregularPairsEvaluated,  // Count of neighbor pairs computed (work counter)

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

    // Worker-side timers (Phase 8)
    WorkerRecvWait,      // Time waiting in MPI_Recv for next task
    WorkerTaskDispatch,  // Time executing the assigned task
    WorkerSendComplete,  // Time in MPI_Isend + MPI_Wait
    WorkerIdleTime,      // Cumulative idle time (recv wait + sync waits)

    // Queue scheduler timers (Phase 8)
    QueueAssign,      // Time in assignQueueAuto()
    QueueRun,         // Time in runQueueAuto()
    QueueCallback,    // Time in callback()

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

// Aggregated statistics across MPI ranks (Phase 8)
struct AggregatedStats {
    double min_seconds = 0.0;
    double max_seconds = 0.0;
    double avg_seconds = 0.0;
    double total_seconds = 0.0;
    long long min_count = 0;
    long long max_count = 0;
    long long total_count = 0;
    int min_rank = 0;
    int max_rank = 0;
    double load_balance_ratio = 1.0;  // max/avg, >1.2 indicates imbalance
    double throughput = 0.0;          // work_units per second (if applicable)
};

// Histogram class for call-time distributions (Phase 8)
class Histogram {
public:
    static constexpr int NUM_BUCKETS = 20;
    // Buckets cover: <1us, 1-10us, 10-100us, 100us-1ms, 1-10ms, 10-100ms, 100ms-1s, 1-10s, >10s
    // Using log10 scale: each bucket covers one order of magnitude

    void record(long long duration_ns) {
        int idx = getBucketIndex(duration_ns);
        buckets_[idx]++;
        total_count_++;
    }

    void reset() {
        for (int i = 0; i < NUM_BUCKETS; ++i) buckets_[i] = 0;
        total_count_ = 0;
    }

    void merge(const Histogram& other) {
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            buckets_[i] += other.buckets_[i];
        }
        total_count_ += other.total_count_;
    }

    // Output
    void print(std::ostream& os) const {
        static const char* labels[] = {
            "<1us", "1-10us", "10-100us", "100us-1ms",
            "1-10ms", "10-100ms", "100ms-1s", "1-10s", ">10s"
        };
        // Map our 20 buckets to 9 display ranges
        long long display[9] = {0};
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            int display_idx = std::min(i / 2, 8);
            display[display_idx] += buckets_[i];
        }
        for (int i = 0; i < 9; ++i) {
            double pct = total_count_ > 0 ? 100.0 * display[i] / total_count_ : 0.0;
            if (display[i] > 0) {
                os << "  " << std::left << std::setw(12) << labels[i]
                   << std::right << std::setw(8) << display[i]
                   << " (" << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%)\n";
            }
        }
        if (total_count_ > 0) {
            os << "  p50: " << formatTime(getPercentile(0.5))
               << ", p90: " << formatTime(getPercentile(0.9))
               << ", p99: " << formatTime(getPercentile(0.99)) << "\n";
        }
    }

    std::string toJSON() const {
        std::ostringstream oss;
        oss << "{\"buckets\":[";
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            if (i > 0) oss << ",";
            oss << buckets_[i];
        }
        oss << "],\"total\":" << total_count_;
        if (total_count_ > 0) {
            oss << ",\"p50_ns\":" << getPercentile(0.5)
                << ",\"p90_ns\":" << getPercentile(0.9)
                << ",\"p99_ns\":" << getPercentile(0.99);
        }
        oss << "}";
        return oss.str();
    }

    // Statistics
    long long getPercentile(double p) const {
        if (total_count_ == 0) return 0;
        long long target = static_cast<long long>(p * total_count_);
        long long cumulative = 0;
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            cumulative += buckets_[i];
            if (cumulative >= target) {
                // Return approximate value at bucket midpoint
                auto range = getBucketRange(i);
                return (range.first + range.second) / 2;
            }
        }
        return getBucketRange(NUM_BUCKETS - 1).second;
    }

    long long getMedian() const { return getPercentile(0.5); }
    long long getCount() const { return total_count_; }

private:
    long long buckets_[NUM_BUCKETS] = {0};
    long long total_count_ = 0;

    int getBucketIndex(long long ns) const {
        // Log10 scale: bucket i covers [10^(i/2) us, 10^((i+1)/2) us)
        // Convert ns to us first
        if (ns < 1000) return 0;  // <1us
        double us = ns / 1000.0;
        int idx = static_cast<int>(2.0 * std::log10(us));
        return std::min(std::max(idx, 0), NUM_BUCKETS - 1);
    }

    std::pair<long long, long long> getBucketRange(int idx) const {
        // Returns range in nanoseconds
        if (idx == 0) return {0, 1000};  // 0 to 1us
        long long low_us = static_cast<long long>(std::pow(10.0, idx / 2.0));
        long long high_us = static_cast<long long>(std::pow(10.0, (idx + 1) / 2.0));
        return {low_us * 1000, high_us * 1000};
    }

    static std::string formatTime(long long ns) {
        if (ns < 1000) return std::to_string(ns) + "ns";
        if (ns < 1000000) return std::to_string(ns / 1000) + "us";
        if (ns < 1000000000) return std::to_string(ns / 1000000) + "ms";
        return std::to_string(ns / 1000000000) + "s";
    }
};

// Online statistics calculator using Welford's algorithm (Phase 15)
// Computes mean, variance, min, max in single pass without storing values
class OnlineStats {
public:
    void update(long long value) {
        count_++;
        double delta = value - mean_;
        mean_ += delta / count_;
        M2_ += delta * (value - mean_);
        min_val_ = std::min(min_val_, value);
        max_val_ = std::max(max_val_, value);
    }

    void reset() {
        count_ = 0;
        mean_ = 0.0;
        M2_ = 0.0;
        min_val_ = LLONG_MAX;
        max_val_ = 0;
    }

    // Accessors
    long long count() const { return count_; }
    double mean() const { return mean_; }
    double variance() const { return count_ > 1 ? M2_ / (count_ - 1) : 0.0; }
    double stddev() const { return std::sqrt(variance()); }
    long long min_val() const { return count_ > 0 ? min_val_ : 0; }
    long long max_val() const { return max_val_; }

    // Outlier detection (> mean + n*stddev)
    bool isOutlier(long long value, double n_sigma = 2.0) const {
        if (count_ < 2) return false;
        return value > mean_ + n_sigma * stddev();
    }

    // For aggregation across MPI ranks
    void merge(const OnlineStats& other) {
        if (other.count_ == 0) return;
        if (count_ == 0) {
            *this = other;
            return;
        }
        // Combined statistics using parallel algorithm
        long long combined_count = count_ + other.count_;
        double delta = other.mean_ - mean_;
        double combined_mean = mean_ + delta * other.count_ / combined_count;
        double combined_M2 = M2_ + other.M2_ + delta * delta * count_ * other.count_ / combined_count;
        count_ = combined_count;
        mean_ = combined_mean;
        M2_ = combined_M2;
        min_val_ = std::min(min_val_, other.min_val_);
        max_val_ = std::max(max_val_, other.max_val_);
    }

private:
    long long count_ = 0;
    double mean_ = 0.0;
    double M2_ = 0.0;  // Sum of squared differences from mean
    long long min_val_ = LLONG_MAX;
    long long max_val_ = 0;
};

// Online correlation tracker using incremental algorithm (Phase 15)
// Computes Pearson correlation coefficient between two variables
class CorrelationTracker {
public:
    void update(double x, double y) {
        n_++;
        sum_x_ += x;
        sum_y_ += y;
        sum_xy_ += x * y;
        sum_x2_ += x * x;
        sum_y2_ += y * y;
    }

    void reset() {
        n_ = 0;
        sum_x_ = sum_y_ = sum_xy_ = sum_x2_ = sum_y2_ = 0.0;
    }

    double correlation() const {
        if (n_ < 2) return 0.0;
        double num = n_ * sum_xy_ - sum_x_ * sum_y_;
        double den_x = n_ * sum_x2_ - sum_x_ * sum_x_;
        double den_y = n_ * sum_y2_ - sum_y_ * sum_y_;
        double den = std::sqrt(den_x * den_y);
        return den > 0 ? num / den : 0.0;
    }

    long long count() const { return n_; }

    // For aggregation across MPI ranks
    void merge(const CorrelationTracker& other) {
        n_ += other.n_;
        sum_x_ += other.sum_x_;
        sum_y_ += other.sum_y_;
        sum_xy_ += other.sum_xy_;
        sum_x2_ += other.sum_x2_;
        sum_y2_ += other.sum_y2_;
    }

private:
    long long n_ = 0;
    double sum_x_ = 0.0, sum_y_ = 0.0;
    double sum_xy_ = 0.0;
    double sum_x2_ = 0.0, sum_y2_ = 0.0;
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

    // Work unit tracking (Phase 8)
    long long work_units = 0;           // e.g., particles or pairs
    long long interval_work_units = 0;

    // Histogram for call-time distributions (Phase 8)
    std::unique_ptr<Histogram> histogram;

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

        if (histogram) histogram->record(duration_ns);
    }

    void resetInterval() {
        interval_total_ns = 0;
        interval_count = 0;
        interval_min_ns = LLONG_MAX;
        interval_max_ns = 0;
        interval_work_units = 0;
        if (histogram) histogram->reset();
    }

    void recordWork(long long units) {
        work_units += units;
        interval_work_units += units;
    }

    void enableHistogram() {
        if (!histogram) histogram = std::unique_ptr<Histogram>(new Histogram());
    }

    double totalSeconds() const { return total_ns * 1e-9; }
    double intervalSeconds() const { return interval_total_ns * 1e-9; }
    double meanNs() const { return count > 0 ? static_cast<double>(total_ns) / count : 0.0; }
    double intervalMeanNs() const {
        return interval_count > 0 ? static_cast<double>(interval_total_ns) / interval_count : 0.0;
    }
    double meanMicros() const { return meanNs() * 1e-3; }
    double intervalMeanMicros() const { return intervalMeanNs() * 1e-3; }
    double throughputPerSecond() const {
        return interval_total_ns > 0 ?
               interval_work_units * 1e9 / interval_total_ns : 0.0;
    }
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
        // Reset neighbor profiling stats (Phase 15)
        interval_neighbor_count_stats_.reset();
        interval_neighbor_time_correlation_.reset();
        interval_outlier_count_ = 0;
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
            // Irregular force sub-timers (Phase 8)
            "IrregularNeighborLoop",
            "IrregularCMLoop",
            "IrregularCorrection",
            "IrregularPredict",
            "IrregularPairsEvaluated",
            // Regular step timers
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
            // Worker-side timers (Phase 8)
            "WorkerRecvWait",
            "WorkerTaskDispatch",
            "WorkerSendComplete",
            "WorkerIdleTime",
            // Queue scheduler timers (Phase 8)
            "QueueAssign",
            "QueueRun",
            "QueueCallback",
            // I/O and other timers
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

    // Record work units for throughput calculation (Phase 8)
    void recordWork(TimerID id, long long units) {
        stats_[static_cast<int>(id)].recordWork(units);
    }

    // Neighbor count profiling (Phase 15)
    void recordNeighborCount(long long count) {
        neighbor_count_stats_.update(count);
        interval_neighbor_count_stats_.update(count);
    }

    void recordNeighborWithTime(long long neighbor_count, long long compute_time_ns) {
        recordNeighborCount(neighbor_count);
        neighbor_time_correlation_.update(static_cast<double>(neighbor_count),
                                          static_cast<double>(compute_time_ns));
        interval_neighbor_time_correlation_.update(static_cast<double>(neighbor_count),
                                                   static_cast<double>(compute_time_ns));
        // Check for outlier
        if (interval_neighbor_count_stats_.count() > 10 &&
            interval_neighbor_count_stats_.isOutlier(neighbor_count, 2.0)) {
            outlier_count_++;
            interval_outlier_count_++;
        }
    }

    const OnlineStats& getNeighborStats() const { return neighbor_count_stats_; }
    const OnlineStats& getIntervalNeighborStats() const { return interval_neighbor_count_stats_; }
    double getNeighborTimeCorrelation() const { return neighbor_time_correlation_.correlation(); }
    double getIntervalNeighborTimeCorrelation() const {
        return interval_neighbor_time_correlation_.correlation();
    }
    long long getOutlierCount() const { return outlier_count_; }
    long long getIntervalOutlierCount() const { return interval_outlier_count_; }

    // Enable histogram for a specific timer (Phase 8)
    void enableHistogram(TimerID id) {
        stats_[static_cast<int>(id)].enableHistogram();
    }

    // Print histograms for all timers that have them enabled (Phase 8)
    void printHistograms(std::ostream& os) const {
        os << "\n==================== Timing Histograms ====================\n";
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            const auto& s = stats_[i];
            if (s.histogram && s.histogram->getCount() > 0) {
                os << "\n" << getTimerName(static_cast<TimerID>(i)) << " histogram:\n";
                s.histogram->print(os);
            }
        }
        os << "============================================================\n";
    }

#ifdef USE_MPI
    // Aggregate statistics across all MPI ranks (Phase 8)
    void aggregateAcrossRanks(MPI_Comm comm) {
        int rank, num_ranks;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &num_ranks);
        num_ranks_ = num_ranks;
        my_rank_ = rank;

        aggregated_stats_.resize(static_cast<int>(TimerID::NUM_TIMERS));

        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            const auto& local = stats_[i];
            auto& agg = aggregated_stats_[i];

            double local_seconds = local.intervalSeconds();

            // Get sum for average
            double sum_seconds = 0.0;
            MPI_Reduce(&local_seconds, &sum_seconds, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

            // Get min with rank using MPI_MINLOC
            struct { double val; int rank; } local_min = {local_seconds, rank};
            struct { double val; int rank; } global_min;
            MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, comm);

            // Get max with rank using MPI_MAXLOC
            struct { double val; int rank; } local_max = {local_seconds, rank};
            struct { double val; int rank; } global_max;
            MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm);

            // Aggregate counts
            long long local_count = local.interval_count;
            long long sum_count = 0;
            MPI_Reduce(&local_count, &sum_count, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);

            // Aggregate work units
            long long local_work = local.interval_work_units;
            long long sum_work = 0;
            MPI_Reduce(&local_work, &sum_work, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);

            if (rank == 0) {
                agg.min_seconds = global_min.val;
                agg.max_seconds = global_max.val;
                agg.avg_seconds = sum_seconds / num_ranks;
                agg.total_seconds = sum_seconds;
                agg.min_rank = global_min.rank;
                agg.max_rank = global_max.rank;
                agg.total_count = sum_count;
                agg.load_balance_ratio = agg.avg_seconds > 0 ? agg.max_seconds / agg.avg_seconds : 1.0;
                agg.throughput = sum_seconds > 0 ? sum_work / sum_seconds : 0.0;
            }
        }
    }

    // Get load balance ratio for a timer (Phase 8)
    double getLoadBalanceRatio(TimerID id) const {
        if (aggregated_stats_.empty()) return 1.0;
        return aggregated_stats_[static_cast<int>(id)].load_balance_ratio;
    }

    // Get aggregated stats for a timer (Phase 8)
    const AggregatedStats& getAggregatedStats(TimerID id) const {
        static AggregatedStats empty;
        if (aggregated_stats_.empty()) return empty;
        return aggregated_stats_[static_cast<int>(id)];
    }

    // Print aggregated summary across MPI ranks (Phase 8)
    void printAggregatedSummary(std::ostream& os, double sim_time_myr) const {
        if (my_rank_ != 0) return;  // Only root prints

        const auto& whole = aggregated_stats_[static_cast<int>(TimerID::WholeRoutine)];
        if (whole.total_seconds == 0) return;

        os << "\n==================== Aggregated Performance Summary ====================\n";
        os << "Simulation Time: " << std::fixed << std::setprecision(6) << sim_time_myr << " Myr\n";
        os << "MPI Ranks: " << num_ranks_ << "\n";
        os << "Total Wall-clock Time: " << std::fixed << std::setprecision(2)
           << whole.max_seconds << " s (slowest rank)\n";
        os << "========================================================================\n\n";

        os << std::left << std::setw(24) << "Timer"
           << std::right << std::setw(10) << "Min (s)"
           << std::setw(10) << "Avg (s)"
           << std::setw(10) << "Max (s)"
           << std::setw(8) << "LB"
           << std::setw(10) << "MinRank"
           << std::setw(10) << "MaxRank"
           << std::setw(12) << "Throughput" << "\n";
        os << std::string(94, '-') << "\n";

        // Print in logical groups
        printAggregatedLine(os, TimerID::IrregularForce);
        printAggregatedLine(os, TimerID::IrregularNeighborLoop);
        printAggregatedLine(os, TimerID::IrregularCMLoop);
        printAggregatedLine(os, TimerID::IrregularCorrection);
        printAggregatedLine(os, TimerID::IrregularPredict);
        printAggregatedLine(os, TimerID::IrregularUpdate);
        os << "\n";

        printAggregatedLine(os, TimerID::RegularForce);
        printAggregatedLine(os, TimerID::RegularGPU);
        printAggregatedLine(os, TimerID::RegularUpdate);
        os << "\n";

        printAggregatedLine(os, TimerID::FewBodyIntegration);
        os << "\n";

        printAggregatedLine(os, TimerID::WorkerRecvWait);
        printAggregatedLine(os, TimerID::WorkerTaskDispatch);
        printAggregatedLine(os, TimerID::WorkerSendComplete);
        printAggregatedLine(os, TimerID::QueueAssign);
        printAggregatedLine(os, TimerID::QueueRun);
        printAggregatedLine(os, TimerID::QueueWait);
        os << "\n";

        os << std::string(94, '=') << "\n";

        // Print load balance warnings
        os << "\n--- Load Balance Analysis ---\n";
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            const auto& agg = aggregated_stats_[i];
            if (agg.load_balance_ratio > 1.2 && agg.max_seconds > 0.1) {
                os << "WARNING: " << getTimerName(static_cast<TimerID>(i))
                   << " has load imbalance ratio " << std::fixed << std::setprecision(2)
                   << agg.load_balance_ratio << " (max on rank " << agg.max_rank << ")\n";
            }
        }

        // Print throughput metrics
        const auto& pairs = aggregated_stats_[static_cast<int>(TimerID::IrregularPairsEvaluated)];
        if (pairs.throughput > 0) {
            os << "\n--- Throughput Metrics ---\n";
            os << "Neighbor pairs evaluated: " << std::scientific << std::setprecision(2)
               << pairs.throughput << " pairs/s\n";
        }
    }
#endif

private:
    Profiler() {
        stats_.resize(static_cast<int>(TimerID::NUM_TIMERS));
        start_times_.resize(static_cast<int>(TimerID::NUM_TIMERS));

        // Enable histograms for key timers (Phase 8)
        enableHistogramForTimer(TimerID::IrregularForce);
        enableHistogramForTimer(TimerID::RegularForce);
        enableHistogramForTimer(TimerID::FewBodyIntegration);
        enableHistogramForTimer(TimerID::QueueWait);
        enableHistogramForTimer(TimerID::WorkerCompute);
        enableHistogramForTimer(TimerID::IrregularNeighborLoop);
        enableHistogramForTimer(TimerID::WorkerRecvWait);
    }

    void enableHistogramForTimer(TimerID id) {
        stats_[static_cast<int>(id)].enableHistogram();
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

#ifdef USE_MPI
    void printAggregatedLine(std::ostream& os, TimerID id) const {
        const auto& agg = aggregated_stats_[static_cast<int>(id)];
        if (agg.total_count == 0 && agg.total_seconds == 0) return;

        os << std::left << std::setw(24) << getTimerName(id)
           << std::right << std::fixed << std::setprecision(3)
           << std::setw(10) << agg.min_seconds
           << std::setw(10) << agg.avg_seconds
           << std::setw(10) << agg.max_seconds
           << std::setprecision(2) << std::setw(8) << agg.load_balance_ratio
           << std::setw(10) << agg.min_rank
           << std::setw(10) << agg.max_rank;
        if (agg.throughput > 0) {
            os << std::scientific << std::setprecision(1) << std::setw(12) << agg.throughput;
        } else {
            os << std::setw(12) << "-";
        }
        os << "\n";
    }
#endif

    std::vector<TimerStats> stats_;
    std::vector<std::chrono::high_resolution_clock::time_point> start_times_;
#ifdef USE_MPI
    std::vector<AggregatedStats> aggregated_stats_;
    int num_ranks_ = 1;
    int my_rank_ = 0;
#endif

    // Neighbor profiling (Phase 15)
    OnlineStats neighbor_count_stats_;
    OnlineStats interval_neighbor_count_stats_;
    CorrelationTracker neighbor_time_correlation_;
    CorrelationTracker interval_neighbor_time_correlation_;
    long long outlier_count_ = 0;
    long long interval_outlier_count_ = 0;
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
#define PROFILE_WORK(timer_id, units) Profiler::instance().recordWork(timer_id, units)
#define PROFILE_NEIGHBOR(count) Profiler::instance().recordNeighborCount(count)
#define PROFILE_NEIGHBOR_TIME(count, time_ns) Profiler::instance().recordNeighborWithTime(count, time_ns)

#else

#define PROFILE_SCOPE(timer_id) ((void)0)
#define PROFILE_START(timer_id) ((void)0)
#define PROFILE_STOP(timer_id) ((void)0)
#define PROFILE_COUNT(timer_id) ((void)0)
#define PROFILE_COUNT_N(timer_id, n) ((void)0)
#define PROFILE_WORK(timer_id, units) ((void)0)
#define PROFILE_NEIGHBOR(count) ((void)0)
#define PROFILE_NEIGHBOR_TIME(count, time_ns) ((void)0)

#endif

// Global accessor for convenience
inline Profiler& profiler() {
    return Profiler::instance();
}

#endif // PROFILER_H
