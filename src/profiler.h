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

// Phase 19: Linux perf_event for cache counters
#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <cstring>
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

// Queue depth tracker for monitoring pending tasks (Phase 16)
// Samples queue depth to identify dispatch bottlenecks
class QueueDepthTracker {
public:
    void sample(int depth) {
        sample_count_++;
        depth_stats_.update(depth);
        if (depth == 0) {
            empty_count_++;
        }
        // Track if queue was non-empty when sampled (for starvation correlation)
        last_depth_ = depth;
    }

    void reset() {
        sample_count_ = 0;
        empty_count_ = 0;
        depth_stats_.reset();
        last_depth_ = 0;
    }

    // Accessors
    long long sampleCount() const { return sample_count_; }
    long long emptyCount() const { return empty_count_; }
    int lastDepth() const { return last_depth_; }
    double meanDepth() const { return depth_stats_.mean(); }
    double stddevDepth() const { return depth_stats_.stddev(); }
    long long minDepth() const { return depth_stats_.min_val(); }
    long long maxDepth() const { return depth_stats_.max_val(); }
    double emptyRatio() const {
        return sample_count_ > 0 ? static_cast<double>(empty_count_) / sample_count_ : 0.0;
    }

private:
    long long sample_count_ = 0;
    long long empty_count_ = 0;
    int last_depth_ = 0;
    OnlineStats depth_stats_;
};

// Heavy particle info for tracking outliers (Phase 17)
struct HeavyParticleInfo {
    int particle_id;
    int worker_rank;
    long long compute_time_ns;
    long long neighbor_count;
};

// Worker distribution tracker for load balance analysis (Phase 17)
// Tracks per-worker particle counts and compute times
class WorkerDistributionTracker {
public:
    void initialize(int num_workers) {
        num_workers_ = num_workers;
        particles_per_worker_.resize(num_workers + 1, 0);
        compute_time_per_worker_ns_.resize(num_workers + 1, 0);
        heavy_particles_.reserve(100);  // Cap heavy particles per interval
    }

    void recordAssignment(int worker_rank) {
        if (worker_rank > 0 && worker_rank <= num_workers_) {
            particles_per_worker_[worker_rank]++;
        }
    }

    void recordComputeTime(int worker_rank, long long compute_ns) {
        if (worker_rank > 0 && worker_rank <= num_workers_) {
            compute_time_per_worker_ns_[worker_rank] += compute_ns;
            worker_time_stats_.update(compute_ns);
        }
    }

    void recordHeavyParticle(int particle_id, int worker_rank,
                             long long compute_time_ns, long long neighbor_count) {
        if (heavy_particles_.size() < 100) {  // Cap at 100 per interval
            heavy_particles_.push_back({particle_id, worker_rank,
                                       compute_time_ns, neighbor_count});
        }
    }

    void reset() {
        std::fill(particles_per_worker_.begin(), particles_per_worker_.end(), 0);
        std::fill(compute_time_per_worker_ns_.begin(), compute_time_per_worker_ns_.end(), 0);
        heavy_particles_.clear();
        worker_time_stats_.reset();
        particle_count_stats_.reset();
    }

    // Compute statistics after interval completes
    void computeStats() const {
        particle_count_stats_.reset();
        for (int i = 1; i <= num_workers_; i++) {
            if (particles_per_worker_[i] > 0) {
                particle_count_stats_.update(particles_per_worker_[i]);
            }
        }
    }

    // Accessors
    int getNumWorkers() const { return num_workers_; }
    int getParticleCount(int worker_rank) const {
        return (worker_rank > 0 && worker_rank <= num_workers_)
               ? particles_per_worker_[worker_rank] : 0;
    }
    long long getComputeTimeNs(int worker_rank) const {
        return (worker_rank > 0 && worker_rank <= num_workers_)
               ? compute_time_per_worker_ns_[worker_rank] : 0;
    }
    double getComputeTimeSeconds(int worker_rank) const {
        return getComputeTimeNs(worker_rank) * 1e-9;
    }

    // Aggregate statistics
    int getTotalParticles() const {
        int total = 0;
        for (int i = 1; i <= num_workers_; i++) {
            total += particles_per_worker_[i];
        }
        return total;
    }

    const OnlineStats& getParticleCountStats() const { return particle_count_stats_; }

    double getLoadBalanceRatio() const {
        if (num_workers_ == 0) return 1.0;
        long long max_time = 0;
        long long total_time = 0;
        int active_workers = 0;
        for (int i = 1; i <= num_workers_; i++) {
            if (compute_time_per_worker_ns_[i] > 0) {
                max_time = std::max(max_time, compute_time_per_worker_ns_[i]);
                total_time += compute_time_per_worker_ns_[i];
                active_workers++;
            }
        }
        if (active_workers == 0 || total_time == 0) return 1.0;
        double avg_time = static_cast<double>(total_time) / active_workers;
        return max_time / avg_time;
    }

    int getMaxTimeWorker() const {
        int max_worker = 0;
        long long max_time = 0;
        for (int i = 1; i <= num_workers_; i++) {
            if (compute_time_per_worker_ns_[i] > max_time) {
                max_time = compute_time_per_worker_ns_[i];
                max_worker = i;
            }
        }
        return max_worker;
    }

    double getMinComputeTimeSeconds() const {
        long long min_time = LLONG_MAX;
        for (int i = 1; i <= num_workers_; i++) {
            if (compute_time_per_worker_ns_[i] > 0) {
                min_time = std::min(min_time, compute_time_per_worker_ns_[i]);
            }
        }
        return (min_time == LLONG_MAX) ? 0.0 : min_time * 1e-9;
    }

    double getMaxComputeTimeSeconds() const {
        long long max_time = 0;
        for (int i = 1; i <= num_workers_; i++) {
            max_time = std::max(max_time, compute_time_per_worker_ns_[i]);
        }
        return max_time * 1e-9;
    }

    double getMeanComputeTimeSeconds() const {
        long long total_time = 0;
        int active_workers = 0;
        for (int i = 1; i <= num_workers_; i++) {
            if (compute_time_per_worker_ns_[i] > 0) {
                total_time += compute_time_per_worker_ns_[i];
                active_workers++;
            }
        }
        return (active_workers > 0) ? (total_time * 1e-9 / active_workers) : 0.0;
    }

    const std::vector<HeavyParticleInfo>& getHeavyParticles() const {
        return heavy_particles_;
    }

    bool hasSignificantImbalance() const {
        return getLoadBalanceRatio() > 1.5;
    }

    // Phase 20: Get histogram of per-worker compute times (ANLYS-02)
    Histogram getWorkerTimeHistogram() const {
        Histogram hist;
        for (int i = 1; i <= num_workers_; i++) {
            if (compute_time_per_worker_ns_[i] > 0) {
                hist.record(compute_time_per_worker_ns_[i]);
            }
        }
        return hist;
    }

private:
    int num_workers_ = 0;
    std::vector<int> particles_per_worker_;
    std::vector<long long> compute_time_per_worker_ns_;
    std::vector<HeavyParticleInfo> heavy_particles_;
    OnlineStats worker_time_stats_;
    mutable OnlineStats particle_count_stats_;
};

// Particle type breakdown for CM vs regular particles (Phase 18)
struct ParticleTypeStats {
    // Counts
    int cm_particle_count = 0;
    int regular_particle_count = 0;

    // Timing (nanoseconds)
    long long cm_compute_time_ns = 0;
    long long regular_compute_time_ns = 0;

    void reset() {
        cm_particle_count = 0;
        regular_particle_count = 0;
        cm_compute_time_ns = 0;
        regular_compute_time_ns = 0;
    }

    // Accessors
    int getTotalParticles() const {
        return cm_particle_count + regular_particle_count;
    }

    double getCMRatio() const {
        int total = getTotalParticles();
        return (total > 0) ? static_cast<double>(cm_particle_count) / total : 0.0;
    }

    double getCMTimeRatio() const {
        long long total = cm_compute_time_ns + regular_compute_time_ns;
        return (total > 0) ? static_cast<double>(cm_compute_time_ns) / total : 0.0;
    }

    double getCMTimeSeconds() const { return cm_compute_time_ns * 1e-9; }
    double getRegularTimeSeconds() const { return regular_compute_time_ns * 1e-9; }
    double getTotalTimeSeconds() const {
        return (cm_compute_time_ns + regular_compute_time_ns) * 1e-9;
    }
};

// Phase 19: Cache statistics for memory access profiling
struct CacheStats {
    long long l1d_misses = 0;
    long long l1d_accesses = 0;
    long long ll_misses = 0;     // Last-level (L3) cache
    long long ll_accesses = 0;
    long long measurement_count = 0;
    double total_time_seconds = 0.0;

    void reset() {
        l1d_misses = l1d_accesses = 0;
        ll_misses = ll_accesses = 0;
        measurement_count = 0;
        total_time_seconds = 0.0;
    }

    void accumulate(long long l1d_miss, long long l1d_acc,
                    long long ll_miss, long long ll_acc, double time_s) {
        l1d_misses += l1d_miss;
        l1d_accesses += l1d_acc;
        ll_misses += ll_miss;
        ll_accesses += ll_acc;
        measurement_count++;
        total_time_seconds += time_s;
    }

    double getL1DMissRate() const {
        return l1d_accesses > 0 ? static_cast<double>(l1d_misses) / l1d_accesses : 0.0;
    }

    double getLLMissRate() const {
        return ll_accesses > 0 ? static_cast<double>(ll_misses) / ll_accesses : 0.0;
    }

    // Estimate memory bandwidth (assuming 64-byte cache lines)
    double getEstimatedBandwidthGBs() const {
        if (total_time_seconds <= 0) return 0.0;
        // LL misses go to memory, so LL misses * 64 bytes = memory traffic
        return (ll_misses * 64.0) / total_time_seconds / 1e9;
    }

    // Operational intensity estimate (FLOPs/byte)
    // Assumes ~50 FLOPs per particle pair, 64 bytes per LL miss
    double getOperationalIntensity() const {
        if (ll_misses <= 0) return 0.0;
        // Rough estimate: measurement_count ~ particle pairs processed
        double estimated_flops = measurement_count * 50.0;
        double bytes_transferred = ll_misses * 64.0;
        return estimated_flops / bytes_transferred;
    }

    bool isMemoryBound() const {
        // Machine balance for Skylake ~20 FLOPs/byte
        // If OI < 5, definitely memory-bound
        return getOperationalIntensity() < 5.0;
    }
};

// Phase 19: RAII wrapper for perf_event file descriptors
#ifdef __linux__
class PerfEventCounter {
public:
    PerfEventCounter() : fd_(-1), enabled_(false) {}

    ~PerfEventCounter() {
        if (fd_ >= 0) close(fd_);
    }

    // Initialize counter for cache events
    // cache_id: PERF_COUNT_HW_CACHE_L1D or PERF_COUNT_HW_CACHE_LL
    // op_result: PERF_COUNT_HW_CACHE_RESULT_MISS or PERF_COUNT_HW_CACHE_RESULT_ACCESS
    bool init(int cache_id, int op_id, int op_result) {
        struct perf_event_attr pe;
        memset(&pe, 0, sizeof(pe));
        pe.type = PERF_TYPE_HW_CACHE;
        pe.size = sizeof(pe);
        pe.config = cache_id | (op_id << 8) | (op_result << 16);
        pe.disabled = 1;
        pe.exclude_kernel = 1;
        pe.exclude_hv = 1;

        fd_ = syscall(__NR_perf_event_open, &pe, 0, -1, -1, 0);
        return fd_ >= 0;
    }

    void enable() {
        if (fd_ >= 0) {
            ioctl(fd_, PERF_EVENT_IOC_RESET, 0);
            ioctl(fd_, PERF_EVENT_IOC_ENABLE, 0);
            enabled_ = true;
        }
    }

    long long disable_and_read() {
        if (fd_ < 0 || !enabled_) return -1;
        ioctl(fd_, PERF_EVENT_IOC_DISABLE, 0);
        enabled_ = false;
        long long count = 0;
        if (read(fd_, &count, sizeof(count)) != sizeof(count)) {
            return -1;
        }
        return count;
    }

    bool isValid() const { return fd_ >= 0; }

private:
    int fd_;
    bool enabled_;
};
#endif // __linux__

// Histogram for neighbor count distribution (Phase 15)
// Uses linear buckets for small counts, exponential for large
class NeighborHistogram {
public:
    static constexpr int NUM_BUCKETS = 15;
    // Buckets: 0, 1-5, 6-10, 11-20, 21-50, 51-100, 101-200, 201-500,
    //          501-1000, 1001-2000, 2001-5000, 5001-10000, 10001-20000,
    //          20001-50000, 50001+

    void record(long long count) {
        int idx = getBucketIndex(count);
        buckets_[idx]++;
        total_count_++;
    }

    void reset() {
        for (int i = 0; i < NUM_BUCKETS; ++i) buckets_[i] = 0;
        total_count_ = 0;
    }

    void print(std::ostream& os) const {
        static const char* labels[] = {
            "0", "1-5", "6-10", "11-20", "21-50", "51-100", "101-200",
            "201-500", "501-1K", "1K-2K", "2K-5K", "5K-10K", "10K-20K",
            "20K-50K", "50K+"
        };
        os << "Neighbor count distribution:\n";
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            if (buckets_[i] > 0) {
                double pct = 100.0 * buckets_[i] / total_count_;
                os << "  " << std::left << std::setw(10) << labels[i]
                   << std::right << std::setw(8) << buckets_[i]
                   << " (" << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%)\n";
            }
        }
    }

    std::string toJSON() const {
        std::ostringstream oss;
        oss << "{\"buckets\":[";
        for (int i = 0; i < NUM_BUCKETS; ++i) {
            if (i > 0) oss << ",";
            oss << buckets_[i];
        }
        oss << "],\"total\":" << total_count_ << "}";
        return oss.str();
    }

    long long getCount() const { return total_count_; }

    // Get bucket boundaries for analysis
    static std::pair<long long, long long> getBucketRange(int idx) {
        static const long long bounds[] = {
            0, 1, 6, 11, 21, 51, 101, 201, 501, 1001, 2001, 5001, 10001, 20001, 50001
        };
        static const long long upper[] = {
            0, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, LLONG_MAX
        };
        return {bounds[idx], upper[idx]};
    }

private:
    long long buckets_[NUM_BUCKETS] = {0};
    long long total_count_ = 0;

    int getBucketIndex(long long count) const {
        if (count == 0) return 0;
        if (count <= 5) return 1;
        if (count <= 10) return 2;
        if (count <= 20) return 3;
        if (count <= 50) return 4;
        if (count <= 100) return 5;
        if (count <= 200) return 6;
        if (count <= 500) return 7;
        if (count <= 1000) return 8;
        if (count <= 2000) return 9;
        if (count <= 5000) return 10;
        if (count <= 10000) return 11;
        if (count <= 20000) return 12;
        if (count <= 50000) return 13;
        return 14;  // 50001+
    }
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
        interval_neighbor_histogram_.reset();
        // Reset queue dispatch profiling stats (Phase 16)
        interval_queue_depth_tracker_.reset();
        interval_starvation_events_ = 0;
        interval_dispatch_latency_ns_ = 0;
        interval_dispatch_count_ = 0;
        interval_dispatch_latency_stats_.reset();
        // Reset worker distribution tracking (Phase 17)
        interval_worker_distribution_.reset();
        // Reset particle type stats (Phase 18)
        interval_particle_type_stats_.reset();
        // Reset cache stats (Phase 19)
        interval_cache_stats_.reset();
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

        // Phase 15: Neighbor count statistics
        const auto& ns = interval_neighbor_count_stats_;
        if (ns.count() > 0) {
            os << "\n--- Neighbor Count Statistics ---\n";
            os << "Particles processed: " << ns.count() << "\n";
            os << "Neighbor count: min=" << ns.min_val()
               << ", max=" << ns.max_val()
               << ", mean=" << std::fixed << std::setprecision(1) << ns.mean()
               << ", stddev=" << std::setprecision(1) << ns.stddev() << "\n";
            os << "Outliers (>2σ): " << interval_outlier_count_ << "\n";
            os << "Neighbor-Time correlation: " << std::fixed << std::setprecision(3)
               << interval_neighbor_time_correlation_.correlation() << "\n";
            if (interval_neighbor_histogram_.getCount() > 0) {
                interval_neighbor_histogram_.print(os);
            }
        }

        // Phase 16: Queue dispatch statistics
        const auto& qd = interval_queue_depth_tracker_;
        if (qd.sampleCount() > 0) {
            os << "\n--- Queue Dispatch Statistics ---\n";
            os << "Queue depth: min=" << qd.minDepth()
               << ", max=" << qd.maxDepth()
               << ", mean=" << std::fixed << std::setprecision(1) << qd.meanDepth()
               << ", empty=" << std::setprecision(1) << (qd.emptyRatio() * 100) << "%\n";
            os << "Starvation events: " << interval_starvation_events_ << "\n";
            os << "Dispatch count: " << interval_dispatch_count_ << "\n";
            if (interval_dispatch_count_ > 0) {
                os << "Dispatch latency: mean=" << std::fixed << std::setprecision(0)
                   << getIntervalAvgDispatchLatencyNs() / 1000 << "us"
                   << ", stddev=" << std::setprecision(0)
                   << interval_dispatch_latency_stats_.stddev() / 1000 << "us\n";
            }
            os << "Root time breakdown: assign=" << std::fixed << std::setprecision(1)
               << (getIntervalAssignTimeRatio() * 100) << "%"
               << ", wait=" << std::setprecision(1)
               << (getIntervalWaitTimeRatio() * 100) << "%\n";
            if (isDispatchBottleneck()) {
                os << "WARNING: Dispatch appears to be bottleneck (assign > 50%)\n";
            }
            if (interval_starvation_events_ > 0) {
                os << "WARNING: " << interval_starvation_events_
                   << " starvation events (workers waited with non-empty queue)\n";
            }
        }

        // Phase 17: Worker distribution statistics
        const auto& wd = interval_worker_distribution_;
        if (wd.getTotalParticles() > 0) {
            wd.computeStats();  // Compute particle count stats
            const auto& pcs = wd.getParticleCountStats();

            os << "\n--- Worker Distribution Statistics ---\n";
            os << "Total particles: " << wd.getTotalParticles()
               << ", Workers: " << wd.getNumWorkers() << "\n";
            os << "Particles per worker: min=" << pcs.min_val()
               << ", max=" << pcs.max_val()
               << ", mean=" << std::fixed << std::setprecision(1) << pcs.mean()
               << ", stddev=" << std::setprecision(1) << pcs.stddev() << "\n";
            os << "Compute time per worker: min=" << std::setprecision(3)
               << wd.getMinComputeTimeSeconds() << "s"
               << ", max=" << wd.getMaxComputeTimeSeconds() << "s"
               << ", mean=" << wd.getMeanComputeTimeSeconds() << "s\n";
            os << "Load balance ratio: " << std::setprecision(2)
               << wd.getLoadBalanceRatio();
            if (wd.getMaxTimeWorker() > 0) {
                os << " (max on rank " << wd.getMaxTimeWorker() << ")";
            }
            os << "\n";

            const auto& heavy = wd.getHeavyParticles();
            if (!heavy.empty()) {
                os << "Heavy particles (>2σ): " << heavy.size() << "\n";
                int display_count = std::min(static_cast<int>(heavy.size()), 5);
                for (int i = 0; i < display_count; i++) {
                    os << "  PID " << heavy[i].particle_id
                       << " on rank " << heavy[i].worker_rank
                       << ": " << std::setprecision(1)
                       << (heavy[i].compute_time_ns / 1000.0) << "us"
                       << " (neighbors: " << heavy[i].neighbor_count << ")\n";
                }
                if (heavy.size() > 5) {
                    os << "  ... and " << (heavy.size() - 5) << " more\n";
                }
            }

            if (wd.hasSignificantImbalance()) {
                os << "WARNING: Load imbalance detected (ratio > 1.5)\n";
            }

            // Phase 20: Worker compute time histogram (ANLYS-02)
            os << "Worker compute time distribution:\n";
            wd.getWorkerTimeHistogram().print(os);
        }

        // Phase 18: Particle type breakdown
        const auto& pt = interval_particle_type_stats_;
        if (pt.getTotalParticles() > 0) {
            os << "\n--- Particle Type Breakdown ---\n";
            os << "Regular particles: " << pt.regular_particle_count
               << " (" << std::fixed << std::setprecision(1)
               << ((1.0 - pt.getCMRatio()) * 100) << "%)\n";
            os << "CM particles: " << pt.cm_particle_count
               << " (" << std::setprecision(1)
               << (pt.getCMRatio() * 100) << "%)\n";
            os << "Compute time: regular=" << std::setprecision(3)
               << pt.getRegularTimeSeconds() << "s"
               << ", CM=" << pt.getCMTimeSeconds() << "s\n";
            os << "CM time ratio: " << std::setprecision(1)
               << (pt.getCMTimeRatio() * 100) << "%\n";
        }

        // Phase 18: Few-body timers summary
        const auto& fb_search = stats_[static_cast<int>(TimerID::FewBodySearch)];
        const auto& fb_init = stats_[static_cast<int>(TimerID::FewBodyInitialization)];
        const auto& fb_integ = stats_[static_cast<int>(TimerID::FewBodyIntegration)];
        const auto& fb_term = stats_[static_cast<int>(TimerID::FewBodyTermination)];
        long long fb_total = fb_search.interval_total_ns + fb_init.interval_total_ns +
                            fb_integ.interval_total_ns + fb_term.interval_total_ns;
        if (fb_total > 0) {
            os << "\n--- Few-Body Timing Breakdown ---\n";
            os << "Search: " << std::fixed << std::setprecision(3)
               << (fb_search.interval_total_ns * 1e-9) << "s\n";
            os << "Initialization: " << (fb_init.interval_total_ns * 1e-9) << "s\n";
            os << "Integration: " << (fb_integ.interval_total_ns * 1e-9) << "s\n";
            os << "Termination: " << (fb_term.interval_total_ns * 1e-9) << "s\n";
            os << "Total few-body: " << (fb_total * 1e-9) << "s\n";
        }

        // Phase 19: Cache statistics
        const auto& cs = interval_cache_stats_;
        if (cs.measurement_count > 0) {
            os << "\n--- Cache Statistics ---\n";
            if (areCacheCountersAvailable()) {
                os << "L1D miss rate: " << std::fixed << std::setprecision(2)
                   << (cs.getL1DMissRate() * 100) << "%"
                   << " (" << cs.l1d_misses << "/" << cs.l1d_accesses << ")\n";
                os << "LL (L3) miss rate: " << std::setprecision(2)
                   << (cs.getLLMissRate() * 100) << "%"
                   << " (" << cs.ll_misses << "/" << cs.ll_accesses << ")\n";
                os << "Est. memory bandwidth: " << std::setprecision(2)
                   << cs.getEstimatedBandwidthGBs() << " GB/s\n";
                os << "Operational intensity: " << std::setprecision(2)
                   << cs.getOperationalIntensity() << " FLOPs/byte\n";
                os << "Classification: " << (cs.isMemoryBound() ? "MEMORY-BOUND" : "COMPUTE-BOUND") << "\n";
            } else {
                os << "Cache counters not available (perf_event initialization failed)\n";
            }
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
            // Phase 15: Neighbor profiling columns
            file << ",NeighborCount_count,NeighborCount_min,NeighborCount_max"
                 << ",NeighborCount_mean,NeighborCount_stddev"
                 << ",NeighborCount_outliers,NeighborTime_correlation";
            // Phase 16: Queue dispatch profiling columns
            file << ",QueueDepth_samples,QueueDepth_min,QueueDepth_max"
                 << ",QueueDepth_mean,QueueDepth_empty_ratio"
                 << ",Starvation_events,Dispatch_count"
                 << ",DispatchLatency_mean_ns,DispatchLatency_stddev_ns"
                 << ",AssignTime_ratio,WaitTime_ratio";
            // Phase 17: Worker distribution columns
            file << ",WorkerParticles_total,WorkerParticles_min,WorkerParticles_max"
                 << ",WorkerParticles_mean,WorkerParticles_stddev"
                 << ",WorkerComputeTime_min_s,WorkerComputeTime_max_s,WorkerComputeTime_mean_s"
                 << ",LoadBalanceRatio,MaxTimeWorker,HeavyParticleCount";
            // Phase 18: Particle type breakdown columns
            file << ",RegularParticles,CMParticles,CMRatio"
                 << ",RegularTime_s,CMTime_s,CMTimeRatio";
            // Phase 19: Cache statistics columns
            file << ",L1DMisses,L1DAccesses,L1DMissRate"
                 << ",LLMisses,LLAccesses,LLMissRate"
                 << ",EstBandwidth_GBs,OperationalIntensity,MemoryBound";
            file << "\n";
        }

        file << step << "," << sim_time;
        for (int i = 0; i < static_cast<int>(TimerID::NUM_TIMERS); ++i) {
            file << "," << stats_[i].interval_total_ns
                 << "," << stats_[i].interval_count;
        }
        // Phase 15: Neighbor profiling data
        const auto& ns = interval_neighbor_count_stats_;
        file << "," << ns.count()
             << "," << ns.min_val()
             << "," << ns.max_val()
             << "," << std::fixed << std::setprecision(2) << ns.mean()
             << "," << std::fixed << std::setprecision(2) << ns.stddev()
             << "," << interval_outlier_count_
             << "," << std::fixed << std::setprecision(4) << interval_neighbor_time_correlation_.correlation();
        // Phase 16: Queue dispatch profiling data
        const auto& qd = interval_queue_depth_tracker_;
        file << "," << qd.sampleCount()
             << "," << qd.minDepth()
             << "," << qd.maxDepth()
             << "," << std::fixed << std::setprecision(2) << qd.meanDepth()
             << "," << std::fixed << std::setprecision(4) << qd.emptyRatio()
             << "," << interval_starvation_events_
             << "," << interval_dispatch_count_
             << "," << std::fixed << std::setprecision(0) << getIntervalAvgDispatchLatencyNs()
             << "," << std::fixed << std::setprecision(0) << interval_dispatch_latency_stats_.stddev()
             << "," << std::fixed << std::setprecision(4) << getIntervalAssignTimeRatio()
             << "," << std::fixed << std::setprecision(4) << getIntervalWaitTimeRatio();
        // Phase 17: Worker distribution data
        const auto& wd = interval_worker_distribution_;
        wd.computeStats();
        const auto& pcs = wd.getParticleCountStats();
        file << "," << wd.getTotalParticles()
             << "," << pcs.min_val()
             << "," << pcs.max_val()
             << "," << std::fixed << std::setprecision(2) << pcs.mean()
             << "," << std::fixed << std::setprecision(2) << pcs.stddev()
             << "," << std::fixed << std::setprecision(4) << wd.getMinComputeTimeSeconds()
             << "," << std::fixed << std::setprecision(4) << wd.getMaxComputeTimeSeconds()
             << "," << std::fixed << std::setprecision(4) << wd.getMeanComputeTimeSeconds()
             << "," << std::fixed << std::setprecision(3) << wd.getLoadBalanceRatio()
             << "," << wd.getMaxTimeWorker()
             << "," << wd.getHeavyParticles().size();
        // Phase 18: Particle type breakdown data
        const auto& pt = interval_particle_type_stats_;
        file << "," << pt.regular_particle_count
             << "," << pt.cm_particle_count
             << "," << std::fixed << std::setprecision(4) << pt.getCMRatio()
             << "," << std::fixed << std::setprecision(4) << pt.getRegularTimeSeconds()
             << "," << std::fixed << std::setprecision(4) << pt.getCMTimeSeconds()
             << "," << std::fixed << std::setprecision(4) << pt.getCMTimeRatio();
        // Phase 19: Cache statistics data
        const auto& cs = interval_cache_stats_;
        file << "," << cs.l1d_misses
             << "," << cs.l1d_accesses
             << "," << std::fixed << std::setprecision(4) << cs.getL1DMissRate()
             << "," << cs.ll_misses
             << "," << cs.ll_accesses
             << "," << std::fixed << std::setprecision(4) << cs.getLLMissRate()
             << "," << std::fixed << std::setprecision(4) << cs.getEstimatedBandwidthGBs()
             << "," << std::fixed << std::setprecision(4) << cs.getOperationalIntensity()
             << "," << (cs.isMemoryBound() ? 1 : 0);
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

        file << "  },\n";  // Close timers object

        // Phase 15: Neighbor profiling
        file << "  \"neighbor_profiling\": {\n";
        const auto& ns = interval_neighbor_count_stats_;
        file << "    \"count\": " << ns.count() << ",\n";
        file << "    \"min\": " << ns.min_val() << ",\n";
        file << "    \"max\": " << ns.max_val() << ",\n";
        file << "    \"mean\": " << std::fixed << std::setprecision(2) << ns.mean() << ",\n";
        file << "    \"stddev\": " << std::fixed << std::setprecision(2) << ns.stddev() << ",\n";
        file << "    \"outlier_count\": " << interval_outlier_count_ << ",\n";
        file << "    \"time_correlation\": " << std::fixed << std::setprecision(4)
             << interval_neighbor_time_correlation_.correlation() << ",\n";
        file << "    \"histogram\": " << interval_neighbor_histogram_.toJSON() << "\n";
        file << "  },\n";  // Close neighbor_profiling

        // Phase 16: Queue dispatch profiling
        file << "  \"queue_dispatch\": {\n";
        const auto& qd = interval_queue_depth_tracker_;
        file << "    \"depth_samples\": " << qd.sampleCount() << ",\n";
        file << "    \"depth_min\": " << qd.minDepth() << ",\n";
        file << "    \"depth_max\": " << qd.maxDepth() << ",\n";
        file << "    \"depth_mean\": " << std::fixed << std::setprecision(2) << qd.meanDepth() << ",\n";
        file << "    \"depth_empty_ratio\": " << std::fixed << std::setprecision(4) << qd.emptyRatio() << ",\n";
        file << "    \"starvation_events\": " << interval_starvation_events_ << ",\n";
        file << "    \"dispatch_count\": " << interval_dispatch_count_ << ",\n";
        file << "    \"latency_mean_ns\": " << std::fixed << std::setprecision(0) << getIntervalAvgDispatchLatencyNs() << ",\n";
        file << "    \"latency_stddev_ns\": " << std::fixed << std::setprecision(0) << interval_dispatch_latency_stats_.stddev() << ",\n";
        file << "    \"assign_time_ratio\": " << std::fixed << std::setprecision(4) << getIntervalAssignTimeRatio() << ",\n";
        file << "    \"wait_time_ratio\": " << std::fixed << std::setprecision(4) << getIntervalWaitTimeRatio() << ",\n";
        file << "    \"is_dispatch_bottleneck\": " << (isDispatchBottleneck() ? "true" : "false") << "\n";
        file << "  },\n";  // Close queue_dispatch, add comma for worker_distribution

        // Phase 17: Worker distribution
        file << "  \"worker_distribution\": {\n";
        const auto& wd = interval_worker_distribution_;
        wd.computeStats();
        const auto& pcs = wd.getParticleCountStats();
        file << "    \"total_particles\": " << wd.getTotalParticles() << ",\n";
        file << "    \"num_workers\": " << wd.getNumWorkers() << ",\n";
        file << "    \"particles_per_worker\": {\n";
        file << "      \"min\": " << pcs.min_val() << ",\n";
        file << "      \"max\": " << pcs.max_val() << ",\n";
        file << "      \"mean\": " << std::fixed << std::setprecision(2) << pcs.mean() << ",\n";
        file << "      \"stddev\": " << std::fixed << std::setprecision(2) << pcs.stddev() << "\n";
        file << "    },\n";
        file << "    \"compute_time\": {\n";
        file << "      \"min_s\": " << std::fixed << std::setprecision(4) << wd.getMinComputeTimeSeconds() << ",\n";
        file << "      \"max_s\": " << std::fixed << std::setprecision(4) << wd.getMaxComputeTimeSeconds() << ",\n";
        file << "      \"mean_s\": " << std::fixed << std::setprecision(4) << wd.getMeanComputeTimeSeconds() << "\n";
        file << "    },\n";
        file << "    \"load_balance_ratio\": " << std::fixed << std::setprecision(3) << wd.getLoadBalanceRatio() << ",\n";
        file << "    \"max_time_worker\": " << wd.getMaxTimeWorker() << ",\n";
        file << "    \"is_imbalanced\": " << (wd.hasSignificantImbalance() ? "true" : "false") << ",\n";
        const auto& heavy = wd.getHeavyParticles();
        file << "    \"heavy_particle_count\": " << heavy.size() << ",\n";
        file << "    \"heavy_particles\": [\n";
        for (size_t i = 0; i < heavy.size(); i++) {
            file << "      {\"pid\": " << heavy[i].particle_id
                 << ", \"worker\": " << heavy[i].worker_rank
                 << ", \"time_ns\": " << heavy[i].compute_time_ns
                 << ", \"neighbors\": " << heavy[i].neighbor_count << "}";
            if (i < heavy.size() - 1) file << ",";
            file << "\n";
        }
        file << "    ],\n";  // Close heavy_particles, add comma
        // Phase 20: Worker compute time histogram (ANLYS-02)
        file << "    \"worker_time_histogram\": " << wd.getWorkerTimeHistogram().toJSON() << "\n";
        file << "  },\n";  // Close worker_distribution

        // Phase 18: Particle type breakdown
        file << "  \"particle_type_breakdown\": {\n";
        const auto& pt = interval_particle_type_stats_;
        file << "    \"regular_count\": " << pt.regular_particle_count << ",\n";
        file << "    \"cm_count\": " << pt.cm_particle_count << ",\n";
        file << "    \"cm_ratio\": " << std::fixed << std::setprecision(4) << pt.getCMRatio() << ",\n";
        file << "    \"regular_time_s\": " << std::fixed << std::setprecision(4) << pt.getRegularTimeSeconds() << ",\n";
        file << "    \"cm_time_s\": " << std::fixed << std::setprecision(4) << pt.getCMTimeSeconds() << ",\n";
        file << "    \"cm_time_ratio\": " << std::fixed << std::setprecision(4) << pt.getCMTimeRatio() << "\n";
        file << "  },\n";  // Close particle_type_breakdown

        // Phase 18: Few-body timing breakdown
        file << "  \"fewbody_timing\": {\n";
        const auto& fb_search = stats_[static_cast<int>(TimerID::FewBodySearch)];
        const auto& fb_init = stats_[static_cast<int>(TimerID::FewBodyInitialization)];
        const auto& fb_integ = stats_[static_cast<int>(TimerID::FewBodyIntegration)];
        const auto& fb_term = stats_[static_cast<int>(TimerID::FewBodyTermination)];
        file << "    \"search_s\": " << std::fixed << std::setprecision(4) << (fb_search.interval_total_ns * 1e-9) << ",\n";
        file << "    \"initialization_s\": " << (fb_init.interval_total_ns * 1e-9) << ",\n";
        file << "    \"integration_s\": " << (fb_integ.interval_total_ns * 1e-9) << ",\n";
        file << "    \"termination_s\": " << (fb_term.interval_total_ns * 1e-9) << "\n";
        file << "  },\n";  // Close fewbody_timing, add comma

        // Phase 19: Cache statistics
        file << "  \"cache_statistics\": {\n";
        const auto& cs = interval_cache_stats_;
        file << "    \"l1d_misses\": " << cs.l1d_misses << ",\n";
        file << "    \"l1d_accesses\": " << cs.l1d_accesses << ",\n";
        file << "    \"l1d_miss_rate\": " << std::fixed << std::setprecision(4) << cs.getL1DMissRate() << ",\n";
        file << "    \"ll_misses\": " << cs.ll_misses << ",\n";
        file << "    \"ll_accesses\": " << cs.ll_accesses << ",\n";
        file << "    \"ll_miss_rate\": " << std::fixed << std::setprecision(4) << cs.getLLMissRate() << ",\n";
        file << "    \"estimated_bandwidth_gbs\": " << std::fixed << std::setprecision(4) << cs.getEstimatedBandwidthGBs() << ",\n";
        file << "    \"operational_intensity\": " << std::fixed << std::setprecision(4) << cs.getOperationalIntensity() << ",\n";
        file << "    \"memory_bound\": " << (cs.isMemoryBound() ? "true" : "false") << ",\n";
        file << "    \"counters_available\": " << (areCacheCountersAvailable() ? "true" : "false") << "\n";
        file << "  }\n";  // Close cache_statistics
        file << "}\n";    // Close root object
    }

    // Record work units for throughput calculation (Phase 8)
    void recordWork(TimerID id, long long units) {
        stats_[static_cast<int>(id)].recordWork(units);
    }

    // Neighbor count profiling (Phase 15)
    void recordNeighborCount(long long count) {
        neighbor_count_stats_.update(count);
        interval_neighbor_count_stats_.update(count);
        neighbor_histogram_.record(count);
        interval_neighbor_histogram_.record(count);
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

    // Queue dispatch profiling (Phase 16)
    void sampleQueueDepth(int depth) {
        queue_depth_tracker_.sample(depth);
        interval_queue_depth_tracker_.sample(depth);
    }

    void recordStarvationEvent() {
        starvation_events_++;
        interval_starvation_events_++;
    }

    void recordDispatchLatency(long long latency_ns) {
        total_dispatch_latency_ns_ += latency_ns;
        interval_dispatch_latency_ns_ += latency_ns;
        dispatch_count_++;
        interval_dispatch_count_++;
        // Phase 16: Track latency distribution
        dispatch_latency_stats_.update(latency_ns);
        interval_dispatch_latency_stats_.update(latency_ns);
    }

    // Queue dispatch accessors (Phase 16)
    const QueueDepthTracker& getQueueDepthTracker() const { return queue_depth_tracker_; }
    const QueueDepthTracker& getIntervalQueueDepthTracker() const { return interval_queue_depth_tracker_; }
    long long getStarvationEvents() const { return starvation_events_; }
    long long getIntervalStarvationEvents() const { return interval_starvation_events_; }
    double getIntervalAvgDispatchLatencyNs() const {
        return interval_dispatch_count_ > 0 ?
               static_cast<double>(interval_dispatch_latency_ns_) / interval_dispatch_count_ : 0.0;
    }

    // Dispatch latency statistics (Phase 16)
    const OnlineStats& getDispatchLatencyStats() const { return dispatch_latency_stats_; }
    const OnlineStats& getIntervalDispatchLatencyStats() const { return interval_dispatch_latency_stats_; }
    double getDispatchLatencyStddev() const { return dispatch_latency_stats_.stddev(); }
    double getIntervalDispatchLatencyStddev() const { return interval_dispatch_latency_stats_.stddev(); }

    // Root-side dispatch breakdown (Phase 16)
    // Uses existing QueueAssign and QueueWait timers
    double getIntervalAssignTimeRatio() const {
        const auto& assign = stats_[static_cast<int>(TimerID::QueueAssign)];
        const auto& wait = stats_[static_cast<int>(TimerID::QueueWait)];
        double total = assign.interval_total_ns + wait.interval_total_ns;
        return total > 0 ? static_cast<double>(assign.interval_total_ns) / total : 0.0;
    }

    double getIntervalWaitTimeRatio() const {
        return 1.0 - getIntervalAssignTimeRatio();
    }

    bool isDispatchBottleneck() const {
        // If assign time > 50% of (assign + wait), dispatch is bottleneck
        return getIntervalAssignTimeRatio() > 0.5;
    }

    // Worker distribution tracking (Phase 17)
    void initializeWorkerTracking(int num_workers) {
        worker_distribution_.initialize(num_workers);
        interval_worker_distribution_.initialize(num_workers);
    }

    void recordWorkerAssignment(int worker_rank) {
        worker_distribution_.recordAssignment(worker_rank);
        interval_worker_distribution_.recordAssignment(worker_rank);
        current_particle_worker_rank_ = worker_rank;
    }

    void recordWorkerComputeTime(int worker_rank, long long compute_ns) {
        worker_distribution_.recordComputeTime(worker_rank, compute_ns);
        interval_worker_distribution_.recordComputeTime(worker_rank, compute_ns);
    }

    void recordHeavyParticle(int particle_id, int worker_rank,
                             long long compute_time_ns, long long neighbor_count) {
        worker_distribution_.recordHeavyParticle(particle_id, worker_rank,
                                                 compute_time_ns, neighbor_count);
        interval_worker_distribution_.recordHeavyParticle(particle_id, worker_rank,
                                                          compute_time_ns, neighbor_count);
    }

    int getCurrentParticleWorkerRank() const { return current_particle_worker_rank_; }

    const WorkerDistributionTracker& getWorkerDistribution() const {
        return worker_distribution_;
    }
    const WorkerDistributionTracker& getIntervalWorkerDistribution() const {
        return interval_worker_distribution_;
    }

    void computeWorkerStats() {
        worker_distribution_.computeStats();
        interval_worker_distribution_.computeStats();
    }

    // Phase 19: Initialize cache counters (call once at startup)
    void initCacheCounters() {
#ifdef __linux__
        if (cache_counters_initialized_) return;

        bool l1d_ok = l1d_miss_counter_.init(PERF_COUNT_HW_CACHE_L1D,
                                              PERF_COUNT_HW_CACHE_OP_READ,
                                              PERF_COUNT_HW_CACHE_RESULT_MISS);
        l1d_ok = l1d_ok && l1d_access_counter_.init(PERF_COUNT_HW_CACHE_L1D,
                                                     PERF_COUNT_HW_CACHE_OP_READ,
                                                     PERF_COUNT_HW_CACHE_RESULT_ACCESS);
        bool ll_ok = ll_miss_counter_.init(PERF_COUNT_HW_CACHE_LL,
                                            PERF_COUNT_HW_CACHE_OP_READ,
                                            PERF_COUNT_HW_CACHE_RESULT_MISS);
        ll_ok = ll_ok && ll_access_counter_.init(PERF_COUNT_HW_CACHE_LL,
                                                  PERF_COUNT_HW_CACHE_OP_READ,
                                                  PERF_COUNT_HW_CACHE_RESULT_ACCESS);

        cache_counters_initialized_ = l1d_ok && ll_ok;
        if (!cache_counters_initialized_) {
            std::cerr << "Warning: Could not initialize cache counters. "
                      << "Cache profiling will be disabled.\n";
        }
#endif
    }

    void startCacheCounters() {
#ifdef __linux__
        if (!cache_counters_initialized_) return;
        l1d_miss_counter_.enable();
        l1d_access_counter_.enable();
        ll_miss_counter_.enable();
        ll_access_counter_.enable();
#endif
    }

    void stopCacheCounters(double elapsed_seconds) {
#ifdef __linux__
        if (!cache_counters_initialized_) return;
        long long l1d_miss = l1d_miss_counter_.disable_and_read();
        long long l1d_acc = l1d_access_counter_.disable_and_read();
        long long ll_miss = ll_miss_counter_.disable_and_read();
        long long ll_acc = ll_access_counter_.disable_and_read();

        cache_stats_.accumulate(l1d_miss, l1d_acc, ll_miss, ll_acc, elapsed_seconds);
        interval_cache_stats_.accumulate(l1d_miss, l1d_acc, ll_miss, ll_acc, elapsed_seconds);
#else
        (void)elapsed_seconds;
#endif
    }

    bool areCacheCountersAvailable() const {
#ifdef __linux__
        return cache_counters_initialized_;
#else
        return false;
#endif
    }

    const CacheStats& getIntervalCacheStats() const { return interval_cache_stats_; }

    // Particle type tracking (Phase 18)
    void recordParticleType(bool is_cm_particle, long long compute_time_ns) {
        if (is_cm_particle) {
            particle_type_stats_.cm_particle_count++;
            particle_type_stats_.cm_compute_time_ns += compute_time_ns;
            interval_particle_type_stats_.cm_particle_count++;
            interval_particle_type_stats_.cm_compute_time_ns += compute_time_ns;
        } else {
            particle_type_stats_.regular_particle_count++;
            particle_type_stats_.regular_compute_time_ns += compute_time_ns;
            interval_particle_type_stats_.regular_particle_count++;
            interval_particle_type_stats_.regular_compute_time_ns += compute_time_ns;
        }
    }

    const ParticleTypeStats& getParticleTypeStats() const {
        return particle_type_stats_;
    }

    const ParticleTypeStats& getIntervalParticleTypeStats() const {
        return interval_particle_type_stats_;
    }

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

        // Phase 15: Neighbor count statistics (local, not aggregated)
        const auto& ns = interval_neighbor_count_stats_;
        if (ns.count() > 0) {
            os << "\n--- Neighbor Count Statistics (Root rank) ---\n";
            os << "Particles processed: " << ns.count() << "\n";
            os << "Neighbor count: min=" << ns.min_val()
               << ", max=" << ns.max_val()
               << ", mean=" << std::fixed << std::setprecision(1) << ns.mean()
               << ", stddev=" << std::setprecision(1) << ns.stddev() << "\n";
            os << "Outliers (>2σ): " << interval_outlier_count_ << "\n";
            os << "Neighbor-Time correlation: " << std::fixed << std::setprecision(3)
               << interval_neighbor_time_correlation_.correlation() << "\n";
        }

        // Phase 16: Queue dispatch statistics (local, not aggregated)
        const auto& qd = interval_queue_depth_tracker_;
        if (qd.sampleCount() > 0) {
            os << "\n--- Queue Dispatch Statistics (Root rank) ---\n";
            os << "Queue depth: min=" << qd.minDepth()
               << ", max=" << qd.maxDepth()
               << ", mean=" << std::fixed << std::setprecision(1) << qd.meanDepth() << "\n";
            os << "Starvation events: " << interval_starvation_events_ << "\n";
            os << "Root time breakdown: assign=" << std::fixed << std::setprecision(1)
               << (getIntervalAssignTimeRatio() * 100) << "%"
               << ", wait=" << std::setprecision(1)
               << (getIntervalWaitTimeRatio() * 100) << "%\n";
            if (isDispatchBottleneck()) {
                os << "WARNING: Dispatch appears to be bottleneck\n";
            }
        }

        // Phase 17: Worker distribution statistics
        const auto& wd = interval_worker_distribution_;
        if (wd.getTotalParticles() > 0) {
            wd.computeStats();
            const auto& pcs = wd.getParticleCountStats();

            os << "\n--- Worker Distribution Statistics ---\n";
            os << "Total particles: " << wd.getTotalParticles()
               << ", Workers: " << wd.getNumWorkers() << "\n";
            os << "Particles per worker: min=" << pcs.min_val()
               << ", max=" << pcs.max_val()
               << ", mean=" << std::fixed << std::setprecision(1) << pcs.mean() << "\n";
            os << "Compute time per worker: min=" << std::setprecision(3)
               << wd.getMinComputeTimeSeconds() << "s"
               << ", max=" << wd.getMaxComputeTimeSeconds() << "s"
               << ", mean=" << wd.getMeanComputeTimeSeconds() << "s\n";
            os << "Load balance ratio: " << std::setprecision(2)
               << wd.getLoadBalanceRatio();
            if (wd.getMaxTimeWorker() > 0) {
                os << " (max on rank " << wd.getMaxTimeWorker() << ")";
            }
            os << "\n";
            os << "Heavy particles: " << wd.getHeavyParticles().size() << "\n";
            if (wd.hasSignificantImbalance()) {
                os << "WARNING: Load imbalance detected (ratio > 1.5)\n";
            }
        }

        // Phase 18: Particle type breakdown
        const auto& pt = interval_particle_type_stats_;
        if (pt.getTotalParticles() > 0) {
            os << "\n--- Particle Type Breakdown ---\n";
            os << "Regular: " << pt.regular_particle_count
               << ", CM: " << pt.cm_particle_count
               << " (CM ratio: " << std::fixed << std::setprecision(1)
               << (pt.getCMRatio() * 100) << "%)\n";
            os << "Time ratio: CM=" << std::setprecision(1)
               << (pt.getCMTimeRatio() * 100) << "%\n";
        }

        // Phase 19: Cache statistics
        const auto& cs = interval_cache_stats_;
        if (cs.measurement_count > 0 && areCacheCountersAvailable()) {
            os << "\n--- Cache Statistics ---\n";
            os << "L1D miss rate: " << std::fixed << std::setprecision(2)
               << (cs.getL1DMissRate() * 100) << "%\n";
            os << "LL miss rate: " << std::setprecision(2)
               << (cs.getLLMissRate() * 100) << "%\n";
            os << "Est. bandwidth: " << std::setprecision(2)
               << cs.getEstimatedBandwidthGBs() << " GB/s\n";
            os << "Classification: " << (cs.isMemoryBound() ? "MEMORY-BOUND" : "COMPUTE-BOUND") << "\n";
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
    // Neighbor histogram (Phase 15)
    NeighborHistogram neighbor_histogram_;
    NeighborHistogram interval_neighbor_histogram_;

    // Queue dispatch profiling (Phase 16)
    QueueDepthTracker queue_depth_tracker_;
    QueueDepthTracker interval_queue_depth_tracker_;
    long long starvation_events_ = 0;
    long long interval_starvation_events_ = 0;
    long long total_dispatch_latency_ns_ = 0;
    long long interval_dispatch_latency_ns_ = 0;
    long long dispatch_count_ = 0;
    long long interval_dispatch_count_ = 0;
    // Dispatch latency statistics (Phase 16)
    OnlineStats dispatch_latency_stats_;
    OnlineStats interval_dispatch_latency_stats_;

    // Worker distribution tracking (Phase 17)
    WorkerDistributionTracker worker_distribution_;
    WorkerDistributionTracker interval_worker_distribution_;
    int current_particle_worker_rank_ = 0;  // Tracks worker for current particle

    // Phase 18: Particle type breakdown
    ParticleTypeStats particle_type_stats_;
    ParticleTypeStats interval_particle_type_stats_;

    // Phase 19: Cache statistics
    CacheStats cache_stats_;
    CacheStats interval_cache_stats_;
#ifdef __linux__
    PerfEventCounter l1d_miss_counter_;
    PerfEventCounter l1d_access_counter_;
    PerfEventCounter ll_miss_counter_;
    PerfEventCounter ll_access_counter_;
    bool cache_counters_initialized_ = false;
#endif
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
#define PROFILE_QUEUE_DEPTH(depth) Profiler::instance().sampleQueueDepth(depth)
#define PROFILE_STARVATION_EVENT() Profiler::instance().recordStarvationEvent()
#define PROFILE_DISPATCH_LATENCY(latency_ns) Profiler::instance().recordDispatchLatency(latency_ns)
#define PROFILE_WORKER_ASSIGNMENT(worker_rank) Profiler::instance().recordWorkerAssignment(worker_rank)
#define PROFILE_WORKER_COMPUTE(worker_rank, compute_ns) Profiler::instance().recordWorkerComputeTime(worker_rank, compute_ns)
#define PROFILE_HEAVY_PARTICLE(pid, worker, time_ns, neighbors) Profiler::instance().recordHeavyParticle(pid, worker, time_ns, neighbors)
// Phase 18: Particle type tracking
#define PROFILE_PARTICLE_TYPE(is_cm, compute_ns) Profiler::instance().recordParticleType(is_cm, compute_ns)
// Phase 19: Cache profiling
#define PROFILE_CACHE_INIT() Profiler::instance().initCacheCounters()
#define PROFILE_CACHE_START() Profiler::instance().startCacheCounters()
#define PROFILE_CACHE_STOP(elapsed_s) Profiler::instance().stopCacheCounters(elapsed_s)
#define PROFILE_CACHE_AVAILABLE() Profiler::instance().areCacheCountersAvailable()

#else

#define PROFILE_SCOPE(timer_id) ((void)0)
#define PROFILE_START(timer_id) ((void)0)
#define PROFILE_STOP(timer_id) ((void)0)
#define PROFILE_COUNT(timer_id) ((void)0)
#define PROFILE_COUNT_N(timer_id, n) ((void)0)
#define PROFILE_WORK(timer_id, units) ((void)0)
#define PROFILE_NEIGHBOR(count) ((void)0)
#define PROFILE_NEIGHBOR_TIME(count, time_ns) ((void)0)
#define PROFILE_QUEUE_DEPTH(depth) ((void)0)
#define PROFILE_STARVATION_EVENT() ((void)0)
#define PROFILE_DISPATCH_LATENCY(latency_ns) ((void)0)
#define PROFILE_WORKER_ASSIGNMENT(worker_rank) ((void)0)
#define PROFILE_WORKER_COMPUTE(worker_rank, compute_ns) ((void)0)
#define PROFILE_HEAVY_PARTICLE(pid, worker, time_ns, neighbors) ((void)0)
#define PROFILE_PARTICLE_TYPE(is_cm, compute_ns) ((void)0)
// Phase 19: Cache profiling (no-ops)
#define PROFILE_CACHE_INIT() ((void)0)
#define PROFILE_CACHE_START() ((void)0)
#define PROFILE_CACHE_STOP(elapsed_s) ((void)0)
#define PROFILE_CACHE_AVAILABLE() (false)

#endif

// Global accessor for convenience
inline Profiler& profiler() {
    return Profiler::instance();
}

#endif // PROFILER_H
