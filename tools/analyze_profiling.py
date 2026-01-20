#!/usr/bin/env python3
"""
ABYSS Performance Profiling Analysis Tool

Analyzes profiling data output by ABYSS to identify performance bottlenecks
and optimization opportunities.

Usage:
    python analyze_profiling.py <output_directory>
    python analyze_profiling.py profiling.csv
    python analyze_profiling.py profiling_*.json
"""

import argparse
import json
import os
import sys
from pathlib import Path

try:
    import pandas as pd
    import numpy as np
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("Warning: pandas not available. Some analysis features disabled.")


def load_csv(filepath):
    """Load profiling CSV data."""
    if not HAS_PANDAS:
        print("Error: pandas required for CSV analysis")
        return None
    return pd.read_csv(filepath)


def load_json(filepath):
    """Load profiling JSON data."""
    with open(filepath) as f:
        return json.load(f)


def analyze_csv(df):
    """Analyze profiling CSV data and print summary."""
    print("\n" + "=" * 70)
    print("ABYSS Performance Analysis Summary")
    print("=" * 70)

    # Get timer columns (those ending in _ns)
    timer_cols = [c for c in df.columns if c.endswith('_ns')]

    # Calculate total time
    if 'WholeRoutine_ns' in df.columns:
        total_time = df['WholeRoutine_ns'].sum()
        print(f"\nTotal simulation wall-clock time: {total_time * 1e-9:.2f} seconds")
        print(f"Simulation time range: {df['sim_time_myr'].min():.6f} - {df['sim_time_myr'].max():.6f} Myr")
        print(f"Number of output intervals: {len(df)}")

    print("\n" + "-" * 70)
    print("Time Breakdown by Component")
    print("-" * 70)

    # Calculate totals and percentages
    results = []
    for col in timer_cols:
        total = df[col].sum()
        if total > 0:
            results.append({
                'Timer': col.replace('_ns', ''),
                'Total (s)': total * 1e-9,
                'Percent': 100.0 * total / (df['WholeRoutine_ns'].sum() if 'WholeRoutine_ns' in df.columns else total),
                'Mean (ms)': total * 1e-6 / len(df),
                'Max (ms)': df[col].max() * 1e-6
            })

    # Sort by total time
    results.sort(key=lambda x: x['Total (s)'], reverse=True)

    print(f"\n{'Timer':<30} {'Total (s)':>12} {'Percent':>10} {'Mean (ms)':>12} {'Max (ms)':>12}")
    print("-" * 78)
    for r in results[:15]:  # Top 15
        print(f"{r['Timer']:<30} {r['Total (s)']:>12.3f} {r['Percent']:>9.1f}% {r['Mean (ms)']:>12.2f} {r['Max (ms)']:>12.2f}")

    # Identify bottlenecks
    print("\n" + "-" * 70)
    print("Optimization Recommendations")
    print("-" * 70)

    recommendations = []

    if results:
        top_timer = results[0]
        if top_timer['Percent'] > 50:
            recommendations.append(f"- {top_timer['Timer']} dominates ({top_timer['Percent']:.1f}% of time). Focus optimization here.")

    # Check for MPI overhead
    mpi_timers = [r for r in results if r['Timer'].startswith('MPI')]
    mpi_total = sum(r['Total (s)'] for r in mpi_timers)
    mpi_pct = sum(r['Percent'] for r in mpi_timers)
    if mpi_pct > 10:
        recommendations.append(f"- MPI communication overhead is {mpi_pct:.1f}% of total time. Consider:")
        recommendations.append("  * Reducing message frequency (batching)")
        recommendations.append("  * Using asynchronous communication")
        recommendations.append("  * Better load balancing")

    # Check for queue wait time
    queue_wait = next((r for r in results if r['Timer'] == 'QueueWait'), None)
    if queue_wait and queue_wait['Percent'] > 5:
        recommendations.append(f"- QueueWait is {queue_wait['Percent']:.1f}% of time. Workers may be idle. Consider:")
        recommendations.append("  * Better task distribution")
        recommendations.append("  * Overlapping computation and communication")

    # Check for data structure overhead
    skiplist_timers = [r for r in results if 'SkipList' in r['Timer']]
    skiplist_pct = sum(r['Percent'] for r in skiplist_timers)
    if skiplist_pct > 5:
        recommendations.append(f"- SkipList operations are {skiplist_pct:.1f}% of time. Consider:")
        recommendations.append("  * Alternative data structures")
        recommendations.append("  * Lazy updates")

    if not recommendations:
        recommendations.append("- Performance looks balanced. Consider profiling at finer granularity.")

    for rec in recommendations:
        print(rec)

    return results


def analyze_load_balance(df):
    """Analyze Phase 15-19 load balance metrics (Phase 20: ANLYS-04)."""
    print("\n" + "-" * 70)
    print("Load Balance Analysis (Phases 15-19)")
    print("-" * 70)

    findings = []

    # Phase 15: Neighbor count analysis
    if 'NeighborCount_stddev' in df.columns:
        mean_stddev = df['NeighborCount_stddev'].mean()
        mean_mean = df['NeighborCount_mean'].mean()
        if mean_mean > 0:
            cv = mean_stddev / mean_mean  # Coefficient of variation
            if cv > 0.5:
                findings.append({
                    'source': 'Neighbor Count Variance',
                    'severity': 'HIGH' if cv > 1.0 else 'MEDIUM',
                    'metric': f'CV={cv:.2f} (stddev/mean)',
                    'recommendation': 'Consider work stealing or neighbor-aware scheduling'
                })
            print(f"  Neighbor count: mean={mean_mean:.0f}, stddev={mean_stddev:.0f}, CV={cv:.2f}")

    # Phase 16: Queue dispatch analysis
    if 'Starvation_events' in df.columns:
        total_starvation = df['Starvation_events'].sum()
        if total_starvation > 0:
            findings.append({
                'source': 'Worker Starvation',
                'severity': 'HIGH' if total_starvation > 100 else 'MEDIUM',
                'metric': f'{total_starvation} starvation events',
                'recommendation': 'Consider MPI batching or larger work chunks'
            })
        print(f"  Starvation events: {total_starvation}")

    if 'AssignTime_ratio' in df.columns:
        mean_assign_ratio = df['AssignTime_ratio'].mean()
        if mean_assign_ratio > 0.5:
            findings.append({
                'source': 'Queue Dispatch Overhead',
                'severity': 'HIGH' if mean_assign_ratio > 0.7 else 'MEDIUM',
                'metric': f'Assign ratio={mean_assign_ratio:.1%}',
                'recommendation': 'Root is bottleneck - batch task assignments'
            })
        print(f"  Queue assign time ratio: {mean_assign_ratio:.1%}")

    # Phase 17: Worker distribution analysis
    if 'LoadBalanceRatio' in df.columns:
        mean_lb_ratio = df['LoadBalanceRatio'].mean()
        max_lb_ratio = df['LoadBalanceRatio'].max()
        if mean_lb_ratio > 1.5:
            findings.append({
                'source': 'Worker Load Imbalance',
                'severity': 'HIGH' if mean_lb_ratio > 2.0 else 'MEDIUM',
                'metric': f'Mean LB ratio={mean_lb_ratio:.2f}, Max={max_lb_ratio:.2f}',
                'recommendation': 'Consider dynamic scheduling or work stealing'
            })
        print(f"  Load balance ratio: mean={mean_lb_ratio:.2f}, max={max_lb_ratio:.2f}")

    if 'HeavyParticleCount' in df.columns:
        total_heavy = df['HeavyParticleCount'].sum()
        if total_heavy > 0:
            print(f"  Heavy particles (>2σ): {total_heavy} total")

    # Phase 18: Particle type analysis
    if 'CMTimeRatio' in df.columns:
        mean_cm_time_ratio = df['CMTimeRatio'].mean()
        mean_cm_ratio = df['CMRatio'].mean() if 'CMRatio' in df.columns else 0
        if mean_cm_time_ratio > 0.3:
            findings.append({
                'source': 'CM Particle Overhead',
                'severity': 'MEDIUM',
                'metric': f'CM time ratio={mean_cm_time_ratio:.1%}, CM count ratio={mean_cm_ratio:.1%}',
                'recommendation': 'CM particles are expensive - consider CM-specific optimization'
            })
        print(f"  CM particle time ratio: {mean_cm_time_ratio:.1%}")

    # Phase 19: Memory analysis
    if 'MemoryBound' in df.columns:
        memory_bound_pct = df['MemoryBound'].mean() * 100
        if memory_bound_pct > 50:
            findings.append({
                'source': 'Memory-Bound Execution',
                'severity': 'HIGH' if memory_bound_pct > 80 else 'MEDIUM',
                'metric': f'{memory_bound_pct:.0f}% of intervals memory-bound',
                'recommendation': 'Focus on cache optimization and data locality'
            })
        print(f"  Memory-bound intervals: {memory_bound_pct:.0f}%")

    if 'L1DMissRate' in df.columns:
        mean_l1d_miss = df['L1DMissRate'].mean()
        print(f"  L1D cache miss rate: {mean_l1d_miss:.2%}")

    # Rank findings by severity
    severity_order = {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}
    findings.sort(key=lambda x: severity_order.get(x['severity'], 2))

    return findings


def print_findings(findings):
    """Print ranked findings and recommendations (Phase 20: ANLYS-04)."""
    print("\n" + "-" * 70)
    print("Top Imbalance Sources (Ranked)")
    print("-" * 70)

    if not findings:
        print("  No significant imbalance sources detected.")
        return

    for i, f in enumerate(findings[:3], 1):
        print(f"\n  {i}. [{f['severity']}] {f['source']}")
        print(f"     Metric: {f['metric']}")
        print(f"     Recommendation: {f['recommendation']}")


def generate_report(df, findings, output_path):
    """Generate markdown analysis report (Phase 20: ANLYS-04)."""
    import datetime

    with open(output_path, 'w') as f:
        f.write("# ABYSS Load Balance Analysis Report\n\n")
        f.write(f"**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Summary section
        f.write("## Summary\n\n")
        if 'WholeRoutine_ns' in df.columns:
            total_time = df['WholeRoutine_ns'].sum() * 1e-9
            f.write(f"- Total wall-clock time: {total_time:.2f} seconds\n")
        f.write(f"- Output intervals analyzed: {len(df)}\n")
        if 'sim_time_myr' in df.columns:
            f.write(f"- Simulation time range: {df['sim_time_myr'].min():.6f} - {df['sim_time_myr'].max():.6f} Myr\n")
        f.write("\n")

        # Findings section
        f.write("## Top Imbalance Sources\n\n")
        if not findings:
            f.write("No significant imbalance sources detected.\n\n")
        else:
            for i, finding in enumerate(findings[:3], 1):
                f.write(f"### {i}. {finding['source']} ({finding['severity']})\n\n")
                f.write(f"**Metric:** {finding['metric']}\n\n")
                f.write(f"**Recommendation:** {finding['recommendation']}\n\n")

        # v2.3 Recommendations section
        f.write("## Recommended Optimization Approach for v2.3\n\n")
        if findings:
            top = findings[0]
            if 'Neighbor' in top['source'] or 'Worker Load' in top['source']:
                f.write("Based on the analysis, the primary imbalance source is **work distribution**.\n\n")
                f.write("Recommended approach:\n")
                f.write("1. Implement work stealing between workers\n")
                f.write("2. Consider neighbor-count-aware task scheduling\n")
                f.write("3. Profile specific heavy particles for targeted optimization\n")
            elif 'Queue' in top['source'] or 'Starvation' in top['source']:
                f.write("Based on the analysis, the primary bottleneck is **queue dispatch overhead**.\n\n")
                f.write("Recommended approach:\n")
                f.write("1. Implement MPI message batching\n")
                f.write("2. Consider larger work chunks per assignment\n")
                f.write("3. Evaluate OpenMP hybrid parallelism\n")
            elif 'Memory' in top['source']:
                f.write("Based on the analysis, execution is **memory-bound**.\n\n")
                f.write("Recommended approach:\n")
                f.write("1. Improve data locality in force calculation\n")
                f.write("2. Consider cache-blocking techniques\n")
                f.write("3. Evaluate GPU offload for irregular forces\n")
            else:
                f.write("Based on the analysis, consider the specific recommendations above.\n")
        else:
            f.write("Performance appears balanced. Consider:\n")
            f.write("- Running longer simulations for more data\n")
            f.write("- Profiling at finer granularity\n")

        f.write("\n---\n")
        f.write("*Report generated by ABYSS analyze_profiling.py*\n")

    print(f"\nReport written to: {output_path}")


def analyze_json(data):
    """Analyze a single JSON profiling snapshot.

    Supports both old format (no schema_version) and new format (schema_version 2.3+).
    """
    print("\n" + "=" * 70)
    print(f"Profiling Snapshot - Step {data.get('step', 'N/A')}")
    print(f"Simulation time: {data.get('sim_time_myr', 'N/A')} Myr")

    # Check for schema version (Phase 22)
    schema_version = data.get('schema_version', '1.0')
    print(f"Schema version: {schema_version}")
    print("=" * 70)

    # Display summary section if available (schema 2.3+)
    summary = data.get('summary', {})
    if summary:
        print("\n--- Summary ---")
        print(f"  Wall time: {summary.get('wall_time_s', 'N/A')} s")
        print(f"  Throughput: {summary.get('throughput_particles_per_s', 'N/A'):.0f} particles/s")
        print(f"  Load balance ratio: {summary.get('load_balance_ratio', 'N/A')}")
        print(f"  Primary bottleneck: {summary.get('primary_bottleneck', 'N/A')}")

    timers = data.get('timers', {})

    # Handle both old format (interval_ns in each timer) and new format (just interval_ns)
    # New format filters out zero timers automatically
    whole_time = timers.get('WholeRoutine', {}).get('interval_ns', 1)

    results = []
    for name, stats in timers.items():
        interval_ns = stats.get('interval_ns', 0)
        count = stats.get('count', stats.get('interval_count', 0))

        if interval_ns > 0 or count > 0:
            results.append({
                'Timer': name,
                'Time (s)': interval_ns * 1e-9,
                'Percent': 100.0 * interval_ns / whole_time if whole_time > 0 else 0,
                'Calls': count,
                'Avg (us)': stats.get('mean_ns', interval_ns / count if count > 0 else 0) * 1e-3
            })

    results.sort(key=lambda x: x['Time (s)'], reverse=True)

    print(f"\n{'Timer':<28} {'Time (s)':>10} {'Percent':>8} {'Calls':>10} {'Avg (us)':>12}")
    print("-" * 70)
    for r in results[:20]:
        print(f"{r['Timer']:<28} {r['Time (s)']:>10.4f} {r['Percent']:>7.1f}% {r['Calls']:>10} {r['Avg (us)']:>12.2f}")


def plot_timeline(df, output_path=None):
    """Plot time breakdown over simulation progress."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available for plotting")
        return

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Plot 1: Stacked area chart of time breakdown
    timer_cols = [c for c in df.columns if c.endswith('_ns') and c != 'WholeRoutine_ns']
    data = df[timer_cols].copy()
    data.columns = [c.replace('_ns', '') for c in data.columns]

    # Normalize to percentages
    total = data.sum(axis=1)
    data_pct = data.div(total, axis=0) * 100

    data_pct.plot.area(ax=axes[0], stacked=True, alpha=0.7)
    axes[0].set_xlabel('Output Interval')
    axes[0].set_ylabel('Time (%)')
    axes[0].set_title('Time Distribution Over Simulation')
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Plot 2: Total time per interval
    axes[1].plot(df['sim_time_myr'], df['WholeRoutine_ns'] * 1e-9, 'b-', linewidth=2)
    axes[1].set_xlabel('Simulation Time (Myr)')
    axes[1].set_ylabel('Wall-clock Time (s)')
    axes[1].set_title('Wall-clock Time per Output Interval')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Analyze ABYSS profiling data')
    parser.add_argument('input', nargs='+', help='Profiling CSV file, JSON file(s), or output directory')
    parser.add_argument('--plot', action='store_true', help='Generate plots (requires matplotlib)')
    parser.add_argument('--output', '-o', help='Output file for plot')
    parser.add_argument('--report', '-r', help='Generate markdown report to specified path (Phase 20)')
    args = parser.parse_args()

    for input_path in args.input:
        path = Path(input_path)

        if path.is_dir():
            # Look for profiling.csv in directory
            csv_path = path / 'profiling.csv'
            if csv_path.exists():
                df = load_csv(csv_path)
                if df is not None:
                    results = analyze_csv(df)
                    # Phase 20: Load balance analysis
                    findings = analyze_load_balance(df)
                    print_findings(findings)
                    if args.report:
                        generate_report(df, findings, args.report)
                    if args.plot:
                        plot_path = args.output or str(path / 'profiling_plot.png')
                        plot_timeline(df, plot_path)
            else:
                # Look for JSON files
                json_files = sorted(path.glob('profiling_*.json'))
                if json_files:
                    # Analyze the latest one
                    data = load_json(json_files[-1])
                    analyze_json(data)
                else:
                    print(f"No profiling data found in {path}")

        elif path.suffix == '.csv':
            df = load_csv(path)
            if df is not None:
                results = analyze_csv(df)
                # Phase 20: Load balance analysis
                findings = analyze_load_balance(df)
                print_findings(findings)
                if args.report:
                    generate_report(df, findings, args.report)
                if args.plot:
                    plot_path = args.output or str(path.with_suffix('.png'))
                    plot_timeline(df, plot_path)

        elif path.suffix == '.json':
            data = load_json(path)
            analyze_json(data)

        else:
            print(f"Unknown file type: {path}")


if __name__ == '__main__':
    main()
