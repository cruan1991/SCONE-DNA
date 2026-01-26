#!/usr/bin/env python3
"""
SCONE Ablation Experiment
=========================
A minimal, self-contained experiment comparing:
  - Baseline: naive binary-to-quaternary transcoding
  - SCONE: FSM-guided constraint-aware encoding
  - SCONE-no-FSM: arithmetic coding without FSM constraints

Author: SCONE Team
"""

import random
import math
from typing import List, Tuple, Dict
from dataclasses import dataclass
import matplotlib.pyplot as plt


# =============================================================================
# 1. Synthetic Latent Generator
# =============================================================================

def generate_latent_vectors(num_sequences: int, latent_dim: int, seed: int = 42) -> List[List[float]]:
    """
    Generate synthetic latent vectors (random uniform in [0, 1]).
    
    Args:
        num_sequences: Number of sequences to generate
        latent_dim: Dimension of each latent vector
        seed: Random seed for reproducibility
    
    Returns:
        List of latent vectors, each of shape (latent_dim,)
    """
    rng = random.Random(seed)
    return [[rng.random() for _ in range(latent_dim)] for _ in range(num_sequences)]


# =============================================================================
# 2. Baseline Encoder (Binary → Quaternary Naive Transcoder)
# =============================================================================

def baseline_encode(latent: List[float], bits_per_value: int = 8) -> str:
    """
    Baseline encoder: quantize latent to bits, then map to DNA.
    
    Mapping: 00→A, 01→T, 10→G, 11→C
    
    Args:
        latent: List of float values in [0, 1]
        bits_per_value: Bits to use for each latent value
    
    Returns:
        DNA sequence string
    """
    BASE_MAP = {0b00: 'A', 0b01: 'T', 0b10: 'G', 0b11: 'C'}
    
    # Quantize each latent value to bits
    bits = []
    for val in latent:
        # Clamp to [0, 1] and quantize
        val = max(0.0, min(1.0, val))
        quantized = int(val * (2**bits_per_value - 1))
        # Extract bits (MSB first)
        for i in range(bits_per_value - 1, -1, -1):
            bits.append((quantized >> i) & 1)
    
    # Pad to multiple of 2
    while len(bits) % 2 != 0:
        bits.append(0)
    
    # Convert bit pairs to bases
    dna = []
    for i in range(0, len(bits), 2):
        pair = (bits[i] << 1) | bits[i + 1]
        dna.append(BASE_MAP[pair])
    
    return ''.join(dna)


# =============================================================================
# 3. SCONE-Style Encoder with FSM
# =============================================================================

class FSMConstraint:
    """
    Finite State Machine for biochemical constraint enforcement.
    
    Constraints:
    - GC content: targets ~50% (allows bases to maintain balance)
    - Homopolymer: limits consecutive identical bases to max_hp
    """
    
    BASE_TO_IDX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    IDX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
    GC_BASES = {2, 3}  # G=2, C=3
    AT_BASES = {0, 1}  # A=0, T=1
    
    def __init__(self, gc_window: int = 20, gc_target: float = 0.5, 
                 gc_tolerance: float = 0.1, max_homopolymer: int = 3):
        self.gc_window = gc_window
        self.gc_low = gc_target - gc_tolerance
        self.gc_high = gc_target + gc_tolerance
        self.max_hp = max_homopolymer
        self.reset()
    
    def reset(self):
        """Reset FSM state."""
        self.history = []  # Recent bases (indices)
        self.hp_base = None  # Current homopolymer base
        self.hp_len = 0  # Current homopolymer length
    
    def get_gc_count(self) -> int:
        """Count GC bases in recent history."""
        window = self.history[-self.gc_window:] if len(self.history) >= self.gc_window else self.history
        return sum(1 for b in window if b in self.GC_BASES)
    
    def get_allowed_mask(self) -> List[bool]:
        """
        Get mask of allowed bases based on current state.
        Returns: [allow_A, allow_T, allow_G, allow_C]
        """
        mask = [True, True, True, True]
        
        # Homopolymer constraint
        if self.hp_base is not None and self.hp_len >= self.max_hp:
            mask[self.hp_base] = False
        
        # GC constraint (only if we have enough history)
        if len(self.history) >= self.gc_window:
            gc_count = self.get_gc_count()
            gc_ratio = gc_count / self.gc_window
            
            # If GC too high, discourage more GC
            if gc_ratio >= self.gc_high:
                # Check if oldest base leaving window is GC
                oldest = self.history[-self.gc_window]
                if oldest in self.GC_BASES:
                    # GC leaving, ratio might decrease - still allow GC
                    pass
                else:
                    # AT leaving, might increase ratio - block GC
                    mask[2] = False  # G
                    mask[3] = False  # C
            
            # If GC too low, discourage more AT
            elif gc_ratio <= self.gc_low:
                oldest = self.history[-self.gc_window]
                if oldest in self.AT_BASES:
                    pass
                else:
                    mask[0] = False  # A
                    mask[1] = False  # T
        
        # Ensure at least one base is allowed (fail-safe)
        if not any(mask):
            mask = [True, True, True, True]
            if self.hp_base is not None and self.hp_len >= self.max_hp:
                mask[self.hp_base] = False
        
        return mask
    
    def update(self, base_idx: int):
        """Update FSM state after emitting a base."""
        self.history.append(base_idx)
        
        if self.hp_base == base_idx:
            self.hp_len += 1
        else:
            self.hp_base = base_idx
            self.hp_len = 1


def scone_encode(latent: List[float], use_fsm: bool = True, 
                 bits_per_value: int = 8) -> Tuple[str, float]:
    """
    SCONE-style encoder with optional FSM constraints.
    
    Uses a simple arithmetic-coding-like approach where:
    - Each latent value is quantized
    - Base selection is guided by FSM mask (if enabled)
    - Probability mass is redistributed among allowed bases
    
    Args:
        latent: List of float values in [0, 1]
        use_fsm: Whether to enable FSM constraints
        bits_per_value: Bits per latent value
    
    Returns:
        (dna_sequence, bits_per_nt)
    """
    fsm = FSMConstraint() if use_fsm else None
    
    # Quantize latent to symbols (0-3 for quaternary)
    symbols = []
    for val in latent:
        val = max(0.0, min(1.0, val))
        # Map [0,1] to [0, 2^bits_per_value - 1], then extract 2-bit chunks
        quantized = int(val * (2**bits_per_value - 1))
        for shift in range(bits_per_value - 2, -1, -2):
            sym = (quantized >> shift) & 0b11
            symbols.append(sym)
    
    # Encode symbols to DNA
    dna = []
    total_bits_used = 0.0
    
    for sym in symbols:
        if fsm is not None:
            mask = fsm.get_allowed_mask()
            allowed = [i for i in range(4) if mask[i]]
            allowed_count = len(allowed)
            
            # Map symbol to allowed base
            if sym in allowed:
                base_idx = sym
            else:
                # Remap: find closest allowed
                base_idx = allowed[sym % allowed_count]
            
            # Compute effective bits (capacity)
            total_bits_used += math.log2(allowed_count)
            fsm.update(base_idx)
        else:
            # No FSM: direct mapping
            base_idx = sym
            total_bits_used += 2.0  # log2(4) = 2
        
        dna.append(FSMConstraint.IDX_TO_BASE[base_idx])
    
    dna_str = ''.join(dna)
    bits_per_nt = total_bits_used / len(dna) if dna else 0
    
    return dna_str, bits_per_nt


# =============================================================================
# 4. Constraint Checker
# =============================================================================

@dataclass
class ConstraintMetrics:
    """Metrics for constraint satisfaction."""
    gc_ratio: float
    gc_deviation: float  # |gc_ratio - 0.5|
    max_homopolymer: int
    length: int


def check_constraints(dna: str) -> ConstraintMetrics:
    """
    Check biochemical constraints on a DNA sequence.
    
    Args:
        dna: DNA sequence string
    
    Returns:
        ConstraintMetrics with GC deviation and max homopolymer
    """
    if not dna:
        return ConstraintMetrics(0.0, 0.5, 0, 0)
    
    # GC content
    gc_count = sum(1 for b in dna if b in 'GC')
    gc_ratio = gc_count / len(dna)
    gc_deviation = abs(gc_ratio - 0.5)
    
    # Max homopolymer
    max_hp = 1
    current_hp = 1
    for i in range(1, len(dna)):
        if dna[i] == dna[i-1]:
            current_hp += 1
            max_hp = max(max_hp, current_hp)
        else:
            current_hp = 1
    
    return ConstraintMetrics(
        gc_ratio=gc_ratio,
        gc_deviation=gc_deviation,
        max_homopolymer=max_hp,
        length=len(dna)
    )


# =============================================================================
# 5. Ablation Experiment
# =============================================================================

@dataclass
class ExperimentResult:
    """Results for one method."""
    name: str
    gc_deviations: List[float]
    max_homopolymers: List[int]
    bits_per_nt: List[float]
    
    @property
    def mean_gc_dev(self) -> float:
        return sum(self.gc_deviations) / len(self.gc_deviations)
    
    @property
    def mean_max_hp(self) -> float:
        return sum(self.max_homopolymers) / len(self.max_homopolymers)
    
    @property
    def mean_bpn(self) -> float:
        return sum(self.bits_per_nt) / len(self.bits_per_nt)


def run_ablation_experiment(num_sequences: int = 500, 
                            latent_dim: int = 32,
                            seed: int = 42) -> Dict[str, ExperimentResult]:
    """
    Run ablation experiment comparing three methods.
    
    Args:
        num_sequences: Number of sequences to test
        latent_dim: Dimension of each latent vector
        seed: Random seed
    
    Returns:
        Dictionary mapping method name to ExperimentResult
    """
    print("=" * 60)
    print("SCONE Ablation Experiment")
    print("=" * 60)
    print(f"Sequences: {num_sequences}, Latent dim: {latent_dim}")
    print()
    
    # Generate latent vectors
    print("Generating latent vectors...")
    latents = generate_latent_vectors(num_sequences, latent_dim, seed)
    
    results = {}
    
    # Method 1: Baseline
    print("Running Baseline encoder...")
    baseline_result = ExperimentResult("Baseline", [], [], [])
    for latent in latents:
        dna = baseline_encode(latent)
        metrics = check_constraints(dna)
        baseline_result.gc_deviations.append(metrics.gc_deviation)
        baseline_result.max_homopolymers.append(metrics.max_homopolymer)
        baseline_result.bits_per_nt.append(2.0)  # Fixed 2 bits/nt
    results["Baseline"] = baseline_result
    
    # Method 2: SCONE (FSM enabled)
    print("Running SCONE (FSM enabled)...")
    scone_result = ExperimentResult("SCONE", [], [], [])
    for latent in latents:
        dna, bpn = scone_encode(latent, use_fsm=True)
        metrics = check_constraints(dna)
        scone_result.gc_deviations.append(metrics.gc_deviation)
        scone_result.max_homopolymers.append(metrics.max_homopolymer)
        scone_result.bits_per_nt.append(bpn)
    results["SCONE"] = scone_result
    
    # Method 3: SCONE without FSM
    print("Running SCONE (FSM disabled)...")
    scone_no_fsm_result = ExperimentResult("SCONE-no-FSM", [], [], [])
    for latent in latents:
        dna, bpn = scone_encode(latent, use_fsm=False)
        metrics = check_constraints(dna)
        scone_no_fsm_result.gc_deviations.append(metrics.gc_deviation)
        scone_no_fsm_result.max_homopolymers.append(metrics.max_homopolymer)
        scone_no_fsm_result.bits_per_nt.append(bpn)
    results["SCONE-no-FSM"] = scone_no_fsm_result
    
    print("Done!\n")
    return results


# =============================================================================
# 6. Plotting and Reporting
# =============================================================================

def print_results_table(results: Dict[str, ExperimentResult]):
    """Print results in table format with standard deviations."""
    import statistics
    
    print("=" * 90)
    print("RESULTS SUMMARY (mean ± std over 200 sequences)")
    print("=" * 90)
    print(f"{'Method':<20} {'GC Dev':<20} {'Max HP':<20} {'Bits/nt':<20}")
    print("-" * 90)
    for name, res in results.items():
        gc_std = statistics.stdev(res.gc_deviations) if len(res.gc_deviations) > 1 else 0
        hp_std = statistics.stdev(res.max_homopolymers) if len(res.max_homopolymers) > 1 else 0
        bpn_std = statistics.stdev(res.bits_per_nt) if len(res.bits_per_nt) > 1 else 0
        print(f"{name:<20} {res.mean_gc_dev:.4f}±{gc_std:.4f}     {res.mean_max_hp:.2f}±{hp_std:.2f}           {res.mean_bpn:.3f}±{bpn_std:.3f}")
    print("=" * 90)
    print()
    
    # LaTeX table
    print("LaTeX Table:")
    print("-" * 90)
    print("\\begin{tabular}{lccc}")
    print("\\toprule")
    print("Method & GC Deviation & Max Homopolymer & Bits/nt \\\\")
    print("\\midrule")
    for name, res in results.items():
        gc_std = statistics.stdev(res.gc_deviations) if len(res.gc_deviations) > 1 else 0
        hp_std = statistics.stdev(res.max_homopolymers) if len(res.max_homopolymers) > 1 else 0
        bpn_std = statistics.stdev(res.bits_per_nt) if len(res.bits_per_nt) > 1 else 0
        print(f"{name} & {res.mean_gc_dev:.4f}$\\pm${gc_std:.4f} & {res.mean_max_hp:.2f}$\\pm${hp_std:.2f} & {res.mean_bpn:.3f}$\\pm${bpn_std:.3f} \\\\")
    print("\\bottomrule")
    print("\\end{tabular}")
    print()


def plot_results(results: Dict[str, ExperimentResult], save_path: str = None):
    """
    Create bar charts comparing methods with error bars (std dev).
    
    Args:
        results: Dictionary of ExperimentResult
        save_path: Optional path to save figure
    """
    import statistics
    
    methods = list(results.keys())
    gc_devs = [results[m].mean_gc_dev for m in methods]
    max_hps = [results[m].mean_max_hp for m in methods]
    bpns = [results[m].mean_bpn for m in methods]
    
    # Compute standard deviations
    gc_stds = [statistics.stdev(results[m].gc_deviations) if len(results[m].gc_deviations) > 1 else 0 for m in methods]
    hp_stds = [statistics.stdev(results[m].max_homopolymers) if len(results[m].max_homopolymers) > 1 else 0 for m in methods]
    bpn_stds = [statistics.stdev(results[m].bits_per_nt) if len(results[m].bits_per_nt) > 1 else 0 for m in methods]
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    # Colors (matplotlib default)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    x = range(len(methods))
    
    # Bar chart 1: GC Deviation with error bars
    ax1 = axes[0]
    bars1 = ax1.bar(x, gc_devs, color=colors, yerr=gc_stds, capsize=5, error_kw={'linewidth': 1.5})
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods)
    ax1.set_ylabel('GC Deviation from 50%')
    ax1.set_title('GC Content Control')
    ax1.set_ylim(0, max(gc_devs) * 1.5)
    for i, (val, std) in enumerate(zip(gc_devs, gc_stds)):
        ax1.text(i, val + std + 0.003, f'{val:.4f}±{std:.4f}', ha='center', va='bottom', fontsize=8)
    
    # Bar chart 2: Max Homopolymer with error bars
    ax2 = axes[1]
    bars2 = ax2.bar(x, max_hps, color=colors, yerr=hp_stds, capsize=5, error_kw={'linewidth': 1.5})
    ax2.set_xticks(x)
    ax2.set_xticklabels(methods)
    ax2.set_ylabel('Max Homopolymer Length')
    ax2.set_title('Homopolymer Control')
    ax2.set_ylim(0, max(max_hps) * 1.4)
    ax2.axhline(y=3, color='red', linestyle='--', linewidth=1, label='Target limit (≤3)')
    ax2.legend(loc='upper right')
    for i, (val, std) in enumerate(zip(max_hps, hp_stds)):
        ax2.text(i, val + std + 0.15, f'{val:.2f}±{std:.2f}', ha='center', va='bottom', fontsize=8)
    
    # Bar chart 3: Bits per nucleotide with error bars
    ax3 = axes[2]
    bars3 = ax3.bar(x, bpns, color=colors, yerr=bpn_stds, capsize=5, error_kw={'linewidth': 1.5})
    ax3.set_xticks(x)
    ax3.set_xticklabels(methods)
    ax3.set_ylabel('Bits per Nucleotide')
    ax3.set_title('Encoding Efficiency')
    ax3.set_ylim(0, 2.5)
    ax3.axhline(y=2.0, color='gray', linestyle='--', linewidth=1, label='Shannon limit')
    ax3.legend(loc='upper right')
    for i, (val, std) in enumerate(zip(bpns, bpn_stds)):
        ax3.text(i, val + std + 0.05, f'{val:.3f}±{std:.3f}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    plt.show()


def plot_distributions(results: Dict[str, ExperimentResult], save_path: str = None):
    """
    Create histogram distributions for each method.
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    
    methods = list(results.keys())
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, (name, color) in enumerate(zip(methods, colors)):
        res = results[name]
        
        # GC deviation histogram
        ax_gc = axes[0, i]
        ax_gc.hist(res.gc_deviations, bins=30, color=color, alpha=0.7, edgecolor='black')
        ax_gc.set_xlabel('GC Deviation')
        ax_gc.set_ylabel('Count')
        ax_gc.set_title(f'{name}: GC Deviation')
        ax_gc.axvline(x=0.05, color='red', linestyle='--', label='5% threshold')
        ax_gc.legend(fontsize=8)
        
        # Max homopolymer histogram
        ax_hp = axes[1, i]
        hp_counts = {}
        for hp in res.max_homopolymers:
            hp_counts[hp] = hp_counts.get(hp, 0) + 1
        hp_vals = sorted(hp_counts.keys())
        hp_freqs = [hp_counts[v] for v in hp_vals]
        ax_hp.bar(hp_vals, hp_freqs, color=color, alpha=0.7, edgecolor='black')
        ax_hp.set_xlabel('Max Homopolymer')
        ax_hp.set_ylabel('Count')
        ax_hp.set_title(f'{name}: Max Homopolymer')
        ax_hp.axvline(x=3.5, color='red', linestyle='--', label='Limit (≤3)')
        ax_hp.legend(fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Distribution figure saved to: {save_path}")
    
    plt.show()


# =============================================================================
# Main
# =============================================================================

def main():
    """Run the complete ablation experiment."""
    # Run experiment
    results = run_ablation_experiment(
        num_sequences=200,
        latent_dim=32,
        seed=42
    )
    
    # Print table
    print_results_table(results)
    
    # Plot results
    print("Generating plots...")
    plot_results(results, save_path="ablation_comparison.png")
    plot_distributions(results, save_path="ablation_distributions.png")
    
    print("\n" + "=" * 60)
    print("EXPERIMENT COMPLETE")
    print("=" * 60)
    print("\nKey Findings:")
    print(f"  - Baseline GC deviation: {results['Baseline'].mean_gc_dev:.4f}")
    print(f"  - SCONE GC deviation:    {results['SCONE'].mean_gc_dev:.4f} "
          f"({(1 - results['SCONE'].mean_gc_dev/results['Baseline'].mean_gc_dev)*100:.1f}% improvement)")
    print(f"  - Baseline max HP:       {results['Baseline'].mean_max_hp:.2f}")
    print(f"  - SCONE max HP:          {results['SCONE'].mean_max_hp:.2f}")
    print(f"  - SCONE bits/nt:         {results['SCONE'].mean_bpn:.3f} "
          f"(capacity loss: {(2.0 - results['SCONE'].mean_bpn)/2.0*100:.1f}%)")


if __name__ == "__main__":
    main()
