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


def analyze_json(data):
    """Analyze a single JSON profiling snapshot."""
    print("\n" + "=" * 70)
    print(f"Profiling Snapshot - Step {data.get('step', 'N/A')}")
    print(f"Simulation time: {data.get('sim_time_myr', 'N/A')} Myr")
    print("=" * 70)

    timers = data.get('timers', {})
    whole_time = timers.get('WholeRoutine', {}).get('interval_ns', 1)

    results = []
    for name, stats in timers.items():
        if stats.get('interval_ns', 0) > 0 or stats.get('interval_count', 0) > 0:
            results.append({
                'Timer': name,
                'Time (s)': stats['interval_ns'] * 1e-9,
                'Percent': 100.0 * stats['interval_ns'] / whole_time if whole_time > 0 else 0,
                'Calls': stats['interval_count'],
                'Avg (us)': stats['interval_ns'] * 1e-3 / stats['interval_count'] if stats['interval_count'] > 0 else 0
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
    args = parser.parse_args()

    for input_path in args.input:
        path = Path(input_path)

        if path.is_dir():
            # Look for profiling.csv in directory
            csv_path = path / 'profiling.csv'
            if csv_path.exists():
                df = load_csv(csv_path)
                if df is not None:
                    analyze_csv(df)
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
                analyze_csv(df)
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
