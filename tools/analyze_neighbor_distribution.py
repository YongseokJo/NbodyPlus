#!/usr/bin/env python3
"""
Analyze neighbor count distribution to evaluate queue system effectiveness.

This script answers: Is dynamic dispatch worth the overhead, or would
static partitioning (brute-force parallelization) work just as well?

Key insight: If neighbor counts (and thus compute times) are uniform,
static partitioning would work without dispatch overhead. If variance
is high, dynamic dispatch helps balance load.
"""

import json
import sys
import numpy as np
from pathlib import Path


def analyze_neighbor_distribution(json_path):
    """Analyze neighbor distribution from profiling JSON."""

    with open(json_path) as f:
        data = json.load(f)

    neighbor = data.get('neighbor_profiling', {})
    histogram = neighbor.get('histogram', {}).get('buckets', [])

    if not histogram:
        print("No neighbor histogram data found!")
        return

    # Histogram parameters
    n_buckets = len(histogram)
    max_neighbors = neighbor.get('max', 513)
    bucket_width = max_neighbors / n_buckets

    print("=" * 70)
    print("NEIGHBOR COUNT DISTRIBUTION ANALYSIS")
    print("=" * 70)
    print(f"\nData from: {json_path}")
    print(f"Step: {data.get('step', 'N/A')}")

    # Basic stats
    print(f"\n--- Basic Statistics ---")
    print(f"Total neighbor evaluations: {neighbor.get('count', 0):,}")
    print(f"Min neighbors: {neighbor.get('min', 0)}")
    print(f"Max neighbors: {neighbor.get('max', 0)}")
    print(f"Mean neighbors: {neighbor.get('mean', 0):.2f}")

    # Histogram analysis
    print(f"\n--- Histogram ({n_buckets} buckets, width={bucket_width:.1f}) ---")
    total_samples = sum(histogram)

    # Calculate statistics from histogram
    bucket_centers = [(i + 0.5) * bucket_width for i in range(n_buckets)]
    weighted_sum = sum(h * c for h, c in zip(histogram, bucket_centers))
    if total_samples > 0:
        hist_mean = weighted_sum / total_samples
        weighted_var = sum(h * (c - hist_mean)**2 for h, c in zip(histogram, bucket_centers)) / total_samples
        hist_std = np.sqrt(weighted_var)
    else:
        hist_mean = hist_std = 0

    print(f"Total samples in histogram: {total_samples:,}")
    print(f"Histogram-derived mean: {hist_mean:.2f}")
    print(f"Histogram-derived stddev: {hist_std:.2f}")
    print(f"Coefficient of variation (CV): {hist_std/hist_mean*100:.1f}%" if hist_mean > 0 else "N/A")

    # Print histogram
    print(f"\n--- Neighbor Count Distribution ---")
    print(f"{'Bucket':<12} {'Range':<16} {'Count':>15} {'Percent':>10} {'Cumulative':>12}")
    print("-" * 70)

    cumulative = 0
    for i, count in enumerate(histogram):
        if count > 0:
            low = i * bucket_width
            high = (i + 1) * bucket_width
            pct = 100.0 * count / total_samples if total_samples > 0 else 0
            cumulative += pct
            bar = '█' * int(pct / 2)
            print(f"{i:<12} {low:>6.0f}-{high:<6.0f}  {count:>15,} {pct:>9.1f}% {cumulative:>10.1f}% {bar}")

    # Queue effectiveness analysis
    print(f"\n" + "=" * 70)
    print("QUEUE EFFECTIVENESS ANALYSIS")
    print("=" * 70)

    worker_dist = data.get('worker_distribution', {})
    queue_dispatch = data.get('queue_dispatch', {})

    # Load balance from actual run
    lb_ratio = worker_dist.get('load_balance_ratio', 0)
    particles_stddev = worker_dist.get('particles_per_worker', {}).get('stddev', 0)
    particles_mean = worker_dist.get('particles_per_worker', {}).get('mean', 0)

    print(f"\n--- Current Queue System Performance ---")
    print(f"Load balance ratio: {lb_ratio:.4f} (1.0 = perfect, >1.2 = problematic)")
    print(f"Particles per worker: mean={particles_mean:.0f}, stddev={particles_stddev:.0f}")
    print(f"Starvation events: {queue_dispatch.get('starvation_events', 0):,}")
    print(f"Dispatch bottleneck: {queue_dispatch.get('is_dispatch_bottleneck', False)}")

    # Estimate static partitioning impact
    print(f"\n--- Static Partitioning Analysis ---")

    # If we used static partitioning, compute time variance would be proportional
    # to neighbor count variance (since compute ~ neighbors)
    cv_neighbors = hist_std / hist_mean * 100 if hist_mean > 0 else 0

    # Estimated load imbalance under static partitioning
    # Assuming ~6000 particles per worker, central limit theorem applies
    num_workers = worker_dist.get('num_workers', 15)
    particles_total = worker_dist.get('total_particles', 10000000)
    particles_per_worker = particles_total / num_workers if num_workers > 0 else 1

    # Standard error of mean decreases with sqrt(n)
    # But max/mean ratio in samples of n increases
    expected_cv_per_worker = cv_neighbors / np.sqrt(particles_per_worker / 1000)

    # Estimate max/mean ratio for static partitioning
    # Rule of thumb: max ≈ mean + 3*stddev for normal distribution
    estimated_static_imbalance = 1 + 3 * expected_cv_per_worker / 100

    print(f"Neighbor count CV: {cv_neighbors:.1f}%")
    print(f"Particles per worker (static): {particles_per_worker:.0f}")
    print(f"Estimated CV per worker (static): {expected_cv_per_worker:.2f}%")
    print(f"Estimated load imbalance (static): {estimated_static_imbalance:.3f}")

    # Compare costs
    print(f"\n--- Cost-Benefit Analysis ---")

    timers = data.get('timers', {})
    whole_time = timers.get('WholeRoutine', {}).get('interval_ns', 0) / 1e9
    mpi_send = timers.get('MPISend', {}).get('interval_ns', 0) / 1e9
    mpi_recv = timers.get('MPIRecv', {}).get('interval_ns', 0) / 1e9
    queue_wait = timers.get('QueueWait', {}).get('interval_ns', 0) / 1e9
    queue_assign = timers.get('QueueAssign', {}).get('interval_ns', 0) / 1e9
    queue_run = timers.get('QueueRun', {}).get('interval_ns', 0) / 1e9
    queue_callback = timers.get('QueueCallback', {}).get('interval_ns', 0) / 1e9

    dispatch_overhead = mpi_send + mpi_recv + queue_wait + queue_assign + queue_run + queue_callback
    dispatch_pct = 100 * dispatch_overhead / whole_time if whole_time > 0 else 0

    print(f"Current wall time: {whole_time:.2f}s")
    print(f"Dispatch overhead: {dispatch_overhead:.2f}s ({dispatch_pct:.1f}%)")
    print(f"  - MPISend: {mpi_send:.2f}s")
    print(f"  - MPIRecv: {mpi_recv:.2f}s")
    print(f"  - QueueWait: {queue_wait:.2f}s")
    print(f"  - QueueAssign: {queue_assign:.2f}s")
    print(f"  - QueueRun: {queue_run:.2f}s")
    print(f"  - QueueCallback: {queue_callback:.2f}s")

    # Static partitioning estimate
    compute_time = whole_time - dispatch_overhead
    static_time_estimate = compute_time * estimated_static_imbalance

    print(f"\nCompute time (excluding dispatch): {compute_time:.2f}s")
    print(f"Estimated static partitioning time: {static_time_estimate:.2f}s")

    # Verdict
    print(f"\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)

    if static_time_estimate < whole_time:
        savings = whole_time - static_time_estimate
        savings_pct = 100 * savings / whole_time
        print(f"\n⚠️  STATIC PARTITIONING MAY BE FASTER!")
        print(f"   Estimated savings: {savings:.2f}s ({savings_pct:.1f}%)")
        print(f"   Reason: Dispatch overhead ({dispatch_pct:.1f}%) > load imbalance cost")
    else:
        cost = static_time_estimate - whole_time
        cost_pct = 100 * cost / whole_time
        print(f"\n✓  DYNAMIC DISPATCH IS BENEFICIAL")
        print(f"   Avoiding imbalance saves: {cost:.2f}s ({cost_pct:.1f}%)")
        print(f"   Reason: Load imbalance cost > dispatch overhead")

    print(f"\n--- Recommendation ---")
    if dispatch_pct > 25 and cv_neighbors < 50:
        print("→ Consider MPI batching (Phase 24) to reduce dispatch overhead")
        print("  With batching, you get dynamic dispatch benefits with lower overhead")
    elif cv_neighbors > 100:
        print("→ High neighbor variance justifies dynamic dispatch")
        print("→ Focus on optimizing per-particle compute time instead")
    else:
        print("→ Current approach reasonable; batching will improve it further")

    return {
        'cv_neighbors': cv_neighbors,
        'lb_ratio': lb_ratio,
        'dispatch_overhead_pct': dispatch_pct,
        'estimated_static_imbalance': estimated_static_imbalance
    }


def main():
    if len(sys.argv) < 2:
        # Default to latest profiling file
        run_dirs = sorted(Path('workflow/runs').glob('*/work/output/profiling_99.json'))
        if run_dirs:
            json_path = str(run_dirs[-1])
            print(f"Using: {json_path}")
        else:
            print("Usage: python analyze_neighbor_distribution.py <profiling.json>")
            sys.exit(1)
    else:
        json_path = sys.argv[1]

    analyze_neighbor_distribution(json_path)


if __name__ == '__main__':
    main()
