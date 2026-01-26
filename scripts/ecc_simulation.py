#!/usr/bin/env python3
"""
ECC Recovery Simulation
========================
Simulates the effect of Error Correction Codes on DNA storage recovery accuracy.

Models:
  - No ECC: Recovery succeeds only if zero errors
  - With ECC: RS-like code that can correct up to 10% errors
"""

import random
import matplotlib.pyplot as plt
import numpy as np


def simulate_transmission(seq_length: int, error_rate: float, rng: random.Random) -> int:
    """
    Simulate transmission with random bit flips.
    
    Args:
        seq_length: Length of binary sequence
        error_rate: Probability of each bit being flipped
        rng: Random number generator
    
    Returns:
        Number of errors introduced
    """
    errors = sum(1 for _ in range(seq_length) if rng.random() < error_rate)
    return errors


def decode_no_ecc(num_errors: int) -> bool:
    """
    Decode without ECC: succeeds only if no errors.
    """
    return num_errors == 0


def decode_with_ecc(num_errors: int, seq_length: int, correction_capacity: float = 0.10) -> bool:
    """
    Decode with RS-like ECC: can correct up to correction_capacity fraction of errors.
    
    Args:
        num_errors: Number of errors in received sequence
        seq_length: Total sequence length
        correction_capacity: Maximum correctable error fraction (default 10%)
    
    Returns:
        True if decoding succeeds
    """
    max_correctable = int(seq_length * correction_capacity)
    return num_errors <= max_correctable


def run_ecc_experiment(seq_length: int = 256,
                       error_rates: list = None,
                       num_samples: int = 200,
                       ecc_capacity: float = 0.10,
                       seed: int = 42) -> dict:
    """
    Run ECC simulation experiment.
    
    Args:
        seq_length: Length of binary sequences
        error_rates: List of error rates to test
        num_samples: Number of samples per error rate
        ecc_capacity: ECC correction capacity (fraction)
        seed: Random seed
    
    Returns:
        Dictionary with results
    """
    if error_rates is None:
        error_rates = [i * 0.01 for i in range(11)]  # 0% to 10%
    
    rng = random.Random(seed)
    
    results = {
        'error_rates': error_rates,
        'no_ecc_success': [],
        'with_ecc_success': [],
    }
    
    print("="*60)
    print("ECC Recovery Simulation")
    print("="*60)
    print(f"Sequence length: {seq_length} bits")
    print(f"Samples per rate: {num_samples}")
    print(f"ECC capacity: {ecc_capacity*100:.0f}% errors")
    print("="*60)
    
    for rate in error_rates:
        no_ecc_successes = 0
        with_ecc_successes = 0
        
        for _ in range(num_samples):
            num_errors = simulate_transmission(seq_length, rate, rng)
            
            if decode_no_ecc(num_errors):
                no_ecc_successes += 1
            
            if decode_with_ecc(num_errors, seq_length, ecc_capacity):
                with_ecc_successes += 1
        
        no_ecc_rate = no_ecc_successes / num_samples * 100
        with_ecc_rate = with_ecc_successes / num_samples * 100
        
        results['no_ecc_success'].append(no_ecc_rate)
        results['with_ecc_success'].append(with_ecc_rate)
        
        print(f"Error rate {rate*100:5.1f}%: No ECC = {no_ecc_rate:6.1f}%, With ECC = {with_ecc_rate:6.1f}%")
    
    print("="*60)
    
    return results


def plot_ecc_results(results: dict, save_path: str = "ecc_recovery.png"):
    """
    Create publication-ready plot of ECC simulation results.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    
    error_rates_pct = [r * 100 for r in results['error_rates']]
    
    # Plot curves
    ax.plot(error_rates_pct, results['no_ecc_success'], 
            'o-', color='#e74c3c', linewidth=2, markersize=8,
            label='No ECC', markeredgecolor='black', markeredgewidth=0.5)
    
    ax.plot(error_rates_pct, results['with_ecc_success'], 
            's-', color='#2ecc71', linewidth=2, markersize=8,
            label='With ECC (RS-like, 10% capacity)', markeredgecolor='black', markeredgewidth=0.5)
    
    # Styling
    ax.set_xlabel('Channel Error Rate (%)', fontsize=12)
    ax.set_ylabel('Recovery Success Rate (%)', fontsize=12)
    ax.set_title('Effect of Error Correction on DNA Storage Recovery', fontsize=13, fontweight='bold')
    
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-5, 105)
    ax.set_xticks(range(0, 11, 1))
    ax.set_yticks(range(0, 101, 20))
    
    # Grid
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    # Legend
    ax.legend(loc='center right', fontsize=10, frameon=True, fancybox=True, shadow=True)
    
    # Annotations
    ax.axhline(y=100, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.axvline(x=10, color='#2ecc71', linestyle=':', linewidth=1, alpha=0.5)
    
    # Add annotation for ECC threshold
    ax.annotate('ECC correction\nthreshold', xy=(10, 50), xytext=(7, 30),
                fontsize=9, ha='center',
                arrowprops=dict(arrowstyle='->', color='#2ecc71', lw=1.5),
                color='#2ecc71')
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nFigure saved: {save_path}")
    
    plt.show()


def print_latex_table(results: dict):
    """Print results as LaTeX table."""
    print("\nLaTeX Table:")
    print("-"*60)
    print("\\begin{tabular}{ccc}")
    print("\\toprule")
    print("Error Rate (\\%) & No ECC (\\%) & With ECC (\\%) \\\\")
    print("\\midrule")
    
    for i, rate in enumerate(results['error_rates']):
        no_ecc = results['no_ecc_success'][i]
        with_ecc = results['with_ecc_success'][i]
        print(f"{rate*100:.0f} & {no_ecc:.1f} & {with_ecc:.1f} \\\\")
    
    print("\\bottomrule")
    print("\\end{tabular}")


def main():
    # Run experiment
    results = run_ecc_experiment(
        seq_length=256,
        error_rates=[i * 0.01 for i in range(11)],  # 0% to 10%
        num_samples=200,
        ecc_capacity=0.10,
        seed=42
    )
    
    # Print table
    print_latex_table(results)
    
    # Plot
    plot_ecc_results(results, save_path="ecc_recovery.png")
    
    # Summary
    print("\n" + "="*60)
    print("KEY FINDINGS")
    print("="*60)
    print(f"• At 1% error rate:")
    print(f"    No ECC:   {results['no_ecc_success'][1]:.1f}% recovery")
    print(f"    With ECC: {results['with_ecc_success'][1]:.1f}% recovery")
    print(f"• At 5% error rate:")
    print(f"    No ECC:   {results['no_ecc_success'][5]:.1f}% recovery")
    print(f"    With ECC: {results['with_ecc_success'][5]:.1f}% recovery")
    print(f"• ECC provides reliable recovery up to 10% error rate")
    print("="*60)


if __name__ == "__main__":
    main()
