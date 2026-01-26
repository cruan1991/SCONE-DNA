#!/usr/bin/env python
"""
SCONE Camera-Ready Experiments

Run full experiments for SCONE and generate publishable metrics.
"""

import json
import csv
import time
import random
import os
import math
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, asdict
from collections import Counter
import statistics

from fsm_constraint import FSMConstraint, INDEX_TO_BASE, NUM_BASES
from scone_fsm_arith import encode_fsm, decode_fsm


# ==============================================================================
# Configuration
# ==============================================================================

@dataclass
class ExperimentConfig:
    """Experiment configuration"""
    num_sequences: int = 5000
    sequence_length: int = 100
    base_probs: List[float] = None
    gc_window: int = 20
    gc_low: float = 0.45
    gc_high: float = 0.55
    max_homopolymer: int = 3
    random_seed: int = 42
    output_dir: str = "experiment_results"
    
    def __post_init__(self):
        if self.base_probs is None:
            self.base_probs = [0.25, 0.25, 0.25, 0.25]


@dataclass
class SequenceMetrics:
    """Metrics for a single sequence"""
    sequence_id: int
    latent_length: int
    dna_length: int
    bitstream_length: int
    gc_ratio: float
    max_homopolymer: int
    coding_bpn: float      # len(bitstream) / len(dna) - actual coding rate
    capacity_bpn: float    # mean log2(allowed_count) - channel capacity under FSM
    encode_time_ms: float
    decode_time_ms: float
    roundtrip_success: bool


@dataclass
class ExperimentResults:
    """Aggregated experiment results"""
    config: dict
    num_sequences: int
    gc_mean: float
    gc_std: float
    gc_min: float
    gc_max: float
    homopolymer_mean: float
    homopolymer_max: int
    homopolymer_p95: int
    coding_bpn_mean: float
    coding_bpn_std: float
    capacity_bpn_mean: float
    capacity_bpn_std: float
    encode_time_mean_ms: float
    decode_time_mean_ms: float
    total_time_s: float
    roundtrip_success_rate: float


# ==============================================================================
# Measurement Functions
# ==============================================================================

def measure_constraints(dna: str) -> Dict:
    """
    Compute constraint metrics for a DNA sequence.
    
    Args:
        dna: DNA sequence string
    
    Returns:
        Dict containing GC ratio and max homopolymer length
    """
    if not dna:
        return {
            'gc_ratio': 0.0,
            'max_homopolymer': 0,
            'length': 0
        }
    
    # Calculate GC ratio
    gc_count = sum(1 for base in dna if base in 'GC')
    gc_ratio = gc_count / len(dna)
    
    # Calculate max homopolymer length
    max_run = 1
    current_run = 1
    
    for i in range(1, len(dna)):
        if dna[i] == dna[i-1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    
    return {
        'gc_ratio': gc_ratio,
        'max_homopolymer': max_run,
        'length': len(dna)
    }


def evaluate_single(
    latent: List[int],
    base_probs: List[float],
    config: ExperimentConfig,
    sequence_id: int = 0
) -> SequenceMetrics:
    """
    Evaluate a single latent sequence.
    
    Args:
        latent: Latent symbol sequence
        base_probs: Base probability distribution
        config: Experiment configuration
        sequence_id: Sequence ID
    
    Returns:
        SequenceMetrics object
    """
    # Encode
    encode_start = time.perf_counter()
    bitstream, dna = encode_fsm(
        latent,
        base_probs,
        gc_window=config.gc_window,
        gc_low=config.gc_low,
        gc_high=config.gc_high,
        max_homopolymer=config.max_homopolymer
    )
    encode_time = (time.perf_counter() - encode_start) * 1000  # ms
    
    # Measure constraints
    constraints = measure_constraints(dna)
    
    # Calculate rate metrics
    dna_length = len(dna)
    bitstream_length = len(bitstream)
    
    # (1) coding_bpn = len(bitstream) / len(dna) - actual coding rate
    coding_bpn = bitstream_length / dna_length if dna_length > 0 else 0
    
    # (2) capacity_bpn = mean_t log2(allowed_count_t) - channel capacity under FSM
    # Simulate FSM along the latent sequence to compute capacity
    fsm_sim = FSMConstraint(
        gc_window=config.gc_window,
        gc_low=config.gc_low,
        gc_high=config.gc_high,
        max_homopolymer=config.max_homopolymer
    )
    accumulated_log2 = 0.0
    for s in latent:
        mask = fsm_sim.get_mask()
        allowed_count = sum(mask)
        accumulated_log2 += math.log2(allowed_count) if allowed_count > 0 else 0
        fsm_sim.update(s)
    capacity_bpn = accumulated_log2 / len(latent) if len(latent) > 0 else 0
    
    # Decode
    decode_start = time.perf_counter()
    decoded, decoded_dna = decode_fsm(
        bitstream,
        base_probs,
        gc_window=config.gc_window,
        gc_low=config.gc_low,
        gc_high=config.gc_high,
        max_homopolymer=config.max_homopolymer,
        max_symbols=len(latent) + 100
    )
    decode_time = (time.perf_counter() - decode_start) * 1000  # ms
    
    # Verify roundtrip
    roundtrip_success = (decoded == latent)
    
    if not roundtrip_success:
        print(f"WARNING: Sequence {sequence_id} roundtrip failed!")
        print(f"  Input length: {len(latent)}, Decoded length: {len(decoded)}")
    
    return SequenceMetrics(
        sequence_id=sequence_id,
        latent_length=len(latent),
        dna_length=dna_length,
        bitstream_length=bitstream_length,
        gc_ratio=constraints['gc_ratio'],
        max_homopolymer=constraints['max_homopolymer'],
        coding_bpn=coding_bpn,
        capacity_bpn=capacity_bpn,
        encode_time_ms=encode_time,
        decode_time_ms=decode_time,
        roundtrip_success=roundtrip_success
    )


def generate_valid_latent_sequence(
    length: int,
    config: ExperimentConfig,
    rng: random.Random
) -> List[int]:
    """
    Generate a valid latent sequence (respecting FSM constraints).
    
    Args:
        length: Sequence length
        config: Experiment configuration
        rng: Random number generator
    
    Returns:
        Valid latent symbol sequence
    """
    fsm = FSMConstraint(
        gc_window=config.gc_window,
        gc_low=config.gc_low,
        gc_high=config.gc_high,
        max_homopolymer=config.max_homopolymer
    )
    
    sequence = []
    for _ in range(length):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = rng.choice(allowed)
        sequence.append(symbol)
        fsm.update(symbol)
    
    return sequence


# ==============================================================================
# Large-Scale Experiments
# ==============================================================================

def run_large_scale_experiment(config: ExperimentConfig) -> Tuple[ExperimentResults, List[SequenceMetrics]]:
    """
    Run large-scale experiment.
    
    Args:
        config: Experiment configuration
    
    Returns:
        (ExperimentResults, List[SequenceMetrics])
    """
    print("=" * 70)
    print("SCONE Large-Scale Experiment")
    print("=" * 70)
    print(f"Configuration:")
    print(f"  Num sequences: {config.num_sequences}")
    print(f"  Sequence length: {config.sequence_length}")
    print(f"  GC window: {config.gc_window}")
    print(f"  GC range: [{config.gc_low:.2%}, {config.gc_high:.2%}]")
    print(f"  Max homopolymer: {config.max_homopolymer}")
    print(f"  Random seed: {config.random_seed}")
    print("=" * 70)
    
    # Initialize random number generator
    rng = random.Random(config.random_seed)
    
    # Collect metrics
    all_metrics: List[SequenceMetrics] = []
    
    total_start = time.time()
    
    # Progress report interval
    report_interval = max(1, config.num_sequences // 10)
    
    for i in range(config.num_sequences):
        # Generate random latent sequence
        latent = generate_valid_latent_sequence(
            config.sequence_length,
            config,
            rng
        )
        
        # Evaluate
        metrics = evaluate_single(latent, config.base_probs, config, sequence_id=i)
        all_metrics.append(metrics)
        
        # Progress report
        if (i + 1) % report_interval == 0:
            progress = (i + 1) / config.num_sequences * 100
            elapsed = time.time() - total_start
            eta = elapsed / (i + 1) * (config.num_sequences - i - 1)
            print(f"  Progress: {progress:.0f}% ({i+1}/{config.num_sequences}), "
                  f"Elapsed: {elapsed:.1f}s, ETA: {eta:.1f}s")
    
    total_time = time.time() - total_start
    
    print(f"\nExperiment complete. Total time: {total_time:.2f}s")
    
    # Calculate aggregate statistics
    gc_ratios = [m.gc_ratio for m in all_metrics]
    homopolymers = [m.max_homopolymer for m in all_metrics]
    coding_bpns = [m.coding_bpn for m in all_metrics]
    capacity_bpns = [m.capacity_bpn for m in all_metrics]
    encode_times = [m.encode_time_ms for m in all_metrics]
    decode_times = [m.decode_time_ms for m in all_metrics]
    roundtrip_successes = [m.roundtrip_success for m in all_metrics]
    
    # Calculate percentiles
    sorted_homopolymers = sorted(homopolymers)
    p95_idx = int(0.95 * len(sorted_homopolymers))
    
    results = ExperimentResults(
        config=asdict(config),
        num_sequences=config.num_sequences,
        gc_mean=statistics.mean(gc_ratios),
        gc_std=statistics.stdev(gc_ratios) if len(gc_ratios) > 1 else 0,
        gc_min=min(gc_ratios),
        gc_max=max(gc_ratios),
        homopolymer_mean=statistics.mean(homopolymers),
        homopolymer_max=max(homopolymers),
        homopolymer_p95=sorted_homopolymers[p95_idx] if p95_idx < len(sorted_homopolymers) else sorted_homopolymers[-1],
        coding_bpn_mean=statistics.mean(coding_bpns),
        coding_bpn_std=statistics.stdev(coding_bpns) if len(coding_bpns) > 1 else 0,
        capacity_bpn_mean=statistics.mean(capacity_bpns),
        capacity_bpn_std=statistics.stdev(capacity_bpns) if len(capacity_bpns) > 1 else 0,
        encode_time_mean_ms=statistics.mean(encode_times),
        decode_time_mean_ms=statistics.mean(decode_times),
        total_time_s=total_time,
        roundtrip_success_rate=sum(roundtrip_successes) / len(roundtrip_successes)
    )
    
    return results, all_metrics


def compute_histogram(values: List[float], bins: int = 20) -> Dict:
    """Compute histogram of values"""
    if not values:
        return {'bins': [], 'counts': []}
    
    min_val = min(values)
    max_val = max(values)
    
    if min_val == max_val:
        return {'bins': [min_val], 'counts': [len(values)]}
    
    bin_width = (max_val - min_val) / bins
    bin_edges = [min_val + i * bin_width for i in range(bins + 1)]
    
    counts = [0] * bins
    for v in values:
        bin_idx = min(int((v - min_val) / bin_width), bins - 1)
        counts[bin_idx] += 1
    
    return {
        'bins': bin_edges,
        'counts': counts
    }


# ==============================================================================
# Save Results
# ==============================================================================

def save_results(
    results: ExperimentResults,
    all_metrics: List[SequenceMetrics],
    output_dir: str
):
    """
    Save experiment results.
    
    Args:
        results: Aggregated results
        all_metrics: All sequence metrics
        output_dir: Output directory
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save metrics.json
    results_dict = asdict(results)
    
    # Add histogram data
    gc_ratios = [m.gc_ratio for m in all_metrics]
    homopolymers = [m.max_homopolymer for m in all_metrics]
    coding_bpns = [m.coding_bpn for m in all_metrics]
    capacity_bpns = [m.capacity_bpn for m in all_metrics]
    
    results_dict['histograms'] = {
        'gc_ratio': compute_histogram(gc_ratios),
        'max_homopolymer': compute_histogram([float(h) for h in homopolymers], bins=results.homopolymer_max + 1),
        'coding_bpn': compute_histogram(coding_bpns),
        'capacity_bpn': compute_histogram(capacity_bpns)
    }
    
    json_path = os.path.join(output_dir, 'metrics.json')
    with open(json_path, 'w') as f:
        json.dump(results_dict, f, indent=2)
    print(f"Saved metrics to: {json_path}")
    
    # Save CSV
    csv_path = os.path.join(output_dir, 'sequence_metrics.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=asdict(all_metrics[0]).keys())
        writer.writeheader()
        for m in all_metrics:
            writer.writerow(asdict(m))
    print(f"Saved CSV to: {csv_path}")
    
    # Save summary CSV (for plotting)
    summary_csv_path = os.path.join(output_dir, 'summary.csv')
    with open(summary_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Metric', 'Value'])
        writer.writerow(['GC_mean', f'{results.gc_mean:.4f}'])
        writer.writerow(['GC_std', f'{results.gc_std:.4f}'])
        writer.writerow(['GC_min', f'{results.gc_min:.4f}'])
        writer.writerow(['GC_max', f'{results.gc_max:.4f}'])
        writer.writerow(['Homopolymer_mean', f'{results.homopolymer_mean:.2f}'])
        writer.writerow(['Homopolymer_max', f'{results.homopolymer_max}'])
        writer.writerow(['Homopolymer_p95', f'{results.homopolymer_p95}'])
        writer.writerow(['Coding_bpn_mean', f'{results.coding_bpn_mean:.4f}'])
        writer.writerow(['Coding_bpn_std', f'{results.coding_bpn_std:.4f}'])
        writer.writerow(['Capacity_bpn_mean', f'{results.capacity_bpn_mean:.4f}'])
        writer.writerow(['Capacity_bpn_std', f'{results.capacity_bpn_std:.4f}'])
        writer.writerow(['Encode_time_mean_ms', f'{results.encode_time_mean_ms:.4f}'])
        writer.writerow(['Decode_time_mean_ms', f'{results.decode_time_mean_ms:.4f}'])
        writer.writerow(['Roundtrip_success_rate', f'{results.roundtrip_success_rate:.4f}'])
    print(f"Saved summary CSV to: {summary_csv_path}")


def print_summary_table(results: ExperimentResults):
    """
    Print summary table.
    
    Args:
        results: Experiment results
    """
    print("\n" + "=" * 70)
    print("SCONE EXPERIMENT RESULTS SUMMARY")
    print("=" * 70)
    
    print("\n┌─────────────────────────────────────────────────────────────────┐")
    print("│                     SCONE Performance Metrics                    │")
    print("├──────────────────────────────┬──────────────────────────────────┤")
    print(f"│ Number of Sequences          │ {results.num_sequences:>32} │")
    print("├──────────────────────────────┼──────────────────────────────────┤")
    print(f"│ GC Ratio (mean ± std)        │ {results.gc_mean:.4f} ± {results.gc_std:.4f}{' ':>17} │")
    print(f"│ GC Ratio (min, max)          │ ({results.gc_min:.4f}, {results.gc_max:.4f}){' ':>15} │")
    print("├──────────────────────────────┼──────────────────────────────────┤")
    print(f"│ Homopolymer (mean)           │ {results.homopolymer_mean:>32.2f} │")
    print(f"│ Homopolymer (max)            │ {results.homopolymer_max:>32} │")
    print(f"│ Homopolymer (95th pctl)      │ {results.homopolymer_p95:>32} │")
    print("├──────────────────────────────┼──────────────────────────────────┤")
    print(f"│ Coding bpn (mean ± std)      │ {results.coding_bpn_mean:.4f} ± {results.coding_bpn_std:.4f}{' ':>17} │")
    print(f"│ Capacity bpn (mean ± std)    │ {results.capacity_bpn_mean:.4f} ± {results.capacity_bpn_std:.4f}{' ':>17} │")
    print("├──────────────────────────────┼──────────────────────────────────┤")
    print(f"│ Encode Time (mean)           │ {results.encode_time_mean_ms:>29.4f} ms │")
    print(f"│ Decode Time (mean)           │ {results.decode_time_mean_ms:>29.4f} ms │")
    print("├──────────────────────────────┼──────────────────────────────────┤")
    print(f"│ Roundtrip Success Rate       │ {results.roundtrip_success_rate*100:>31.2f}% │")
    print(f"│ Total Experiment Time        │ {results.total_time_s:>30.2f} s │")
    print("└──────────────────────────────┴──────────────────────────────────┘")
    
    # LaTeX output (for paper)
    print("\n" + "-" * 70)
    print("LaTeX Table Row (for paper):")
    print("-" * 70)
    print(f"SCONE & {results.gc_mean:.3f} ± {results.gc_std:.3f} & "
          f"{results.homopolymer_p95} & "
          f"{results.coding_bpn_mean:.3f} & "
          f"{results.capacity_bpn_mean:.3f} & "
          f"{results.encode_time_mean_ms:.2f} \\\\")
    
    # Markdown output
    print("\n" + "-" * 70)
    print("Markdown Table:")
    print("-" * 70)
    print("| Metric | Value |")
    print("|--------|-------|")
    print(f"| GC Ratio | {results.gc_mean:.3f} ± {results.gc_std:.3f} |")
    print(f"| Homopolymer (p95) | {results.homopolymer_p95} |")
    print(f"| Coding bpn | {results.coding_bpn_mean:.3f} ± {results.coding_bpn_std:.3f} |")
    print(f"| Capacity bpn | {results.capacity_bpn_mean:.3f} ± {results.capacity_bpn_std:.3f} |")
    print(f"| Encode Time | {results.encode_time_mean_ms:.2f} ms |")
    print(f"| Decode Time | {results.decode_time_mean_ms:.2f} ms |")
    print(f"| Success Rate | {results.roundtrip_success_rate*100:.1f}% |")


# ==============================================================================
# Main Functions
# ==============================================================================

def run_quick_test():
    """Run quick test (for verification)"""
    print("Running quick test...")
    
    config = ExperimentConfig(
        num_sequences=100,
        sequence_length=50,
        random_seed=42
    )
    
    results, all_metrics = run_large_scale_experiment(config)
    print_summary_table(results)
    
    # Verify all roundtrips succeeded
    assert results.roundtrip_success_rate == 1.0, "Roundtrip verification failed!"
    print("\n✓ Quick test passed")
    
    return results, all_metrics


def run_full_experiment():
    """Run full experiment"""
    print("Running full experiment...")
    
    config = ExperimentConfig(
        num_sequences=5000,
        sequence_length=100,
        random_seed=42,
        output_dir="experiment_results"
    )
    
    results, all_metrics = run_large_scale_experiment(config)
    print_summary_table(results)
    
    # Save results
    save_results(results, all_metrics, config.output_dir)
    
    # Verify all roundtrips succeeded
    assert results.roundtrip_success_rate == 1.0, "Roundtrip verification failed!"
    print("\n✓ Full experiment complete")
    
    return results, all_metrics


def run_parameter_sweep():
    """Run parameter sweep experiment"""
    print("Running parameter sweep experiment...")
    
    # Different GC constraint configurations
    gc_configs = [
        {'gc_low': 0.40, 'gc_high': 0.60, 'name': 'relaxed'},
        {'gc_low': 0.45, 'gc_high': 0.55, 'name': 'standard'},
        {'gc_low': 0.48, 'gc_high': 0.52, 'name': 'strict'},
    ]
    
    # Different homopolymer limits
    homopolymer_limits = [2, 3, 4]
    
    all_results = []
    
    for gc_config in gc_configs:
        for max_hp in homopolymer_limits:
            config = ExperimentConfig(
                num_sequences=1000,
                sequence_length=100,
                gc_low=gc_config['gc_low'],
                gc_high=gc_config['gc_high'],
                max_homopolymer=max_hp,
                random_seed=42,
                output_dir=f"experiment_results/{gc_config['name']}_hp{max_hp}"
            )
            
            print(f"\n--- Config: {gc_config['name']}, max_homopolymer={max_hp} ---")
            results, all_metrics = run_large_scale_experiment(config)
            
            all_results.append({
                'gc_config': gc_config['name'],
                'max_homopolymer': max_hp,
                'gc_mean': results.gc_mean,
                'gc_std': results.gc_std,
                'homopolymer_p95': results.homopolymer_p95,
                'coding_bpn': results.coding_bpn_mean,
                'capacity_bpn': results.capacity_bpn_mean,
                'encode_time_mean_ms': results.encode_time_mean_ms
            })
    
    # Print parameter sweep results table
    print("\n" + "=" * 100)
    print("Parameter Sweep Results")
    print("=" * 100)
    print(f"{'Config':<12} {'HP Limit':<10} {'GC Mean':<10} {'GC Std':<10} {'HP p95':<8} {'Coding':<10} {'Capacity':<10} {'Enc(ms)':<10}")
    print("-" * 100)
    for r in all_results:
        print(f"{r['gc_config']:<12} {r['max_homopolymer']:<10} {r['gc_mean']:.4f}{' ':>4} {r['gc_std']:.4f}{' ':>4} {r['homopolymer_p95']:<8} {r['coding_bpn']:.4f}{' ':>4} {r['capacity_bpn']:.4f}{' ':>4} {r['encode_time_mean_ms']:.4f}")
    
    return all_results


def run_ablation_study():
    """Run ablation study for ISCAS paper"""
    print("=" * 90)
    print("SCONE Ablation Study")
    print("=" * 90)
    
    num_sequences = 1000
    sequence_length = 100
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    results = []
    
    # 1. Full SCONE (standard constraints)
    print("\n--- Full SCONE (GC: 0.45-0.55, HP: 3) ---")
    config = ExperimentConfig(
        num_sequences=num_sequences,
        sequence_length=sequence_length,
        gc_low=0.45,
        gc_high=0.55,
        max_homopolymer=3,
        random_seed=42
    )
    full_results, _ = run_large_scale_experiment(config)
    results.append({
        'config': 'Full SCONE',
        'gc_mean': full_results.gc_mean,
        'gc_std': full_results.gc_std,
        'hp_max': full_results.homopolymer_max,
        'coding_bpn': full_results.coding_bpn_mean,
        'capacity_bpn': full_results.capacity_bpn_mean
    })
    
    # 2. No FSM (unconstrained) - just random sequences
    print("\n--- No FSM (Unconstrained Baseline) ---")
    rng = random.Random(42)
    unconstrained_metrics = []
    
    for _ in range(num_sequences):
        # Generate completely random sequence (no FSM)
        symbols = [rng.randint(0, 3) for _ in range(sequence_length)]
        dna = ''.join(['ATGC'[s] for s in symbols])
        
        constraints = measure_constraints(dna)
        unconstrained_metrics.append({
            'gc_ratio': constraints['gc_ratio'],
            'max_homopolymer': constraints['max_homopolymer']
        })
    
    gc_unconstrained = [m['gc_ratio'] for m in unconstrained_metrics]
    hp_unconstrained = [m['max_homopolymer'] for m in unconstrained_metrics]
    
    # No FSM = full capacity of 2.0 bits/nt (log2(4) = 2)
    results.append({
        'config': 'No FSM',
        'gc_mean': statistics.mean(gc_unconstrained),
        'gc_std': statistics.stdev(gc_unconstrained),
        'hp_max': max(hp_unconstrained),
        'coding_bpn': 2.00,    # No compression overhead
        'capacity_bpn': 2.00   # Full capacity: log2(4) = 2
    })
    
    print(f"  GC Mean: {statistics.mean(gc_unconstrained):.4f} ± {statistics.stdev(gc_unconstrained):.4f}")
    print(f"  HP Max: {max(hp_unconstrained)}")
    print(f"  Capacity: 2.00 bpn (unconstrained)")
    
    # 3. GC only (no HP limit)
    print("\n--- GC Only (no HP constraint) ---")
    config_gc_only = ExperimentConfig(
        num_sequences=num_sequences,
        sequence_length=sequence_length,
        gc_low=0.45,
        gc_high=0.55,
        max_homopolymer=100,  # Effectively no limit
        random_seed=42
    )
    gc_only_results, _ = run_large_scale_experiment(config_gc_only)
    results.append({
        'config': 'GC Only',
        'gc_mean': gc_only_results.gc_mean,
        'gc_std': gc_only_results.gc_std,
        'hp_max': gc_only_results.homopolymer_max,
        'coding_bpn': gc_only_results.coding_bpn_mean,
        'capacity_bpn': gc_only_results.capacity_bpn_mean
    })
    
    # 4. HP only (no GC limit)
    print("\n--- HP Only (no GC constraint) ---")
    config_hp_only = ExperimentConfig(
        num_sequences=num_sequences,
        sequence_length=sequence_length,
        gc_low=0.0,
        gc_high=1.0,  # No GC constraint
        max_homopolymer=3,
        random_seed=42
    )
    hp_only_results, _ = run_large_scale_experiment(config_hp_only)
    results.append({
        'config': 'HP Only',
        'gc_mean': hp_only_results.gc_mean,
        'gc_std': hp_only_results.gc_std,
        'hp_max': hp_only_results.homopolymer_max,
        'coding_bpn': hp_only_results.coding_bpn_mean,
        'capacity_bpn': hp_only_results.capacity_bpn_mean
    })
    
    # Print ablation table
    print("\n" + "=" * 90)
    print("ABLATION STUDY RESULTS")
    print("=" * 90)
    print(f"{'Config':<15} {'GC Mean':<10} {'GC Std':<10} {'HP Max':<10} {'Coding':<10} {'Capacity':<10}")
    print("-" * 90)
    for r in results:
        hp_str = f"{r['hp_max']}+" if r['hp_max'] > 5 else str(r['hp_max'])
        print(f"{r['config']:<15} {r['gc_mean']:.4f}{' ':>4} {r['gc_std']:.4f}{' ':>4} {hp_str:<10} {r['coding_bpn']:.4f}{' ':>4} {r['capacity_bpn']:.4f}")
    
    # LaTeX output for paper
    print("\n" + "-" * 90)
    print("LaTeX Table Rows (for ISCAS paper):")
    print("-" * 90)
    print("\\begin{tabular}{lccccc}")
    print("\\toprule")
    print("\\textbf{Config} & \\textbf{GC Mean} & \\textbf{GC Std} & \\textbf{HP Max} & \\textbf{Coding} & \\textbf{Capacity} \\\\")
    print("\\midrule")
    for r in results:
        hp_str = f"{r['hp_max']}+" if r['hp_max'] > 5 else str(r['hp_max'])
        print(f"{r['config']} & {r['gc_mean']:.3f} & {r['gc_std']:.3f} & {hp_str} & {r['coding_bpn']:.3f} & {r['capacity_bpn']:.3f} \\\\")
    print("\\bottomrule")
    print("\\end{tabular}")
    
    # Markdown table
    print("\n" + "-" * 90)
    print("Markdown Table:")
    print("-" * 90)
    print("| Configuration | GC Mean | GC Std | HP Max | Coding bpn | Capacity bpn |")
    print("|---------------|---------|--------|--------|------------|--------------|")
    for r in results:
        hp_str = f"{r['hp_max']}+" if r['hp_max'] > 5 else str(r['hp_max'])
        print(f"| {r['config']} | {r['gc_mean']:.3f} | {r['gc_std']:.3f} | {hp_str} | {r['coding_bpn']:.3f} | {r['capacity_bpn']:.3f} |")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='SCONE Experiments')
    parser.add_argument('--mode', type=str, default='full',
                       choices=['quick', 'full', 'sweep', 'ablation'],
                       help='Experiment mode: quick, full, sweep, or ablation')
    parser.add_argument('--num', type=int, default=None,
                       help='Override number of sequences')
    parser.add_argument('--length', type=int, default=None,
                       help='Override sequence length')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed')
    parser.add_argument('--output', type=str, default='experiment_results',
                       help='Output directory')
    
    args = parser.parse_args()
    
    if args.mode == 'quick':
        run_quick_test()
    elif args.mode == 'full':
        if args.num is not None or args.length is not None:
            config = ExperimentConfig(
                num_sequences=args.num or 5000,
                sequence_length=args.length or 100,
                random_seed=args.seed,
                output_dir=args.output
            )
            results, all_metrics = run_large_scale_experiment(config)
            print_summary_table(results)
            save_results(results, all_metrics, config.output_dir)
        else:
            run_full_experiment()
    elif args.mode == 'sweep':
        run_parameter_sweep()
    elif args.mode == 'ablation':
        run_ablation_study()
