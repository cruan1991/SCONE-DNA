#!/usr/bin/env python3
"""
SCONE Metrics Visualization
============================
Publication-ready bar charts comparing Baseline, SCONE, and SCONE-no-FSM.

Figures:
  - Fig 2(a): GC deviation + Max homopolymer
  - Fig 2(b): Bits per nucleotide
"""

import matplotlib.pyplot as plt
import numpy as np


def run_experiment(num_sequences: int = 200, latent_dim: int = 32, seed: int = 42):
    """
    Run experiment and collect metrics for all three methods.
    Returns dict with mean and std for each metric.
    """
    import random
    import math
    import statistics
    
    # FSM class (inline for self-contained script)
    class FSM:
        GC = {2, 3}
        AT = {0, 1}
        
        def __init__(self, gc_window=20, gc_tol=0.1, max_hp=3):
            self.gc_window = gc_window
            self.gc_low, self.gc_high = 0.5 - gc_tol, 0.5 + gc_tol
            self.max_hp = max_hp
            self.reset()
        
        def reset(self):
            self.history = []
            self.hp_base = None
            self.hp_len = 0
        
        def get_mask(self):
            mask = [True] * 4
            if self.hp_base is not None and self.hp_len >= self.max_hp:
                mask[self.hp_base] = False
            if len(self.history) >= self.gc_window:
                gc = sum(1 for b in self.history[-self.gc_window:] if b in self.GC)
                ratio = gc / self.gc_window
                if ratio >= self.gc_high:
                    mask[2] = mask[3] = False
                elif ratio <= self.gc_low:
                    mask[0] = mask[1] = False
            if not any(mask):
                mask = [True] * 4
                if self.hp_base is not None and self.hp_len >= self.max_hp:
                    mask[self.hp_base] = False
            return mask
        
        def update(self, b):
            self.history.append(b)
            if self.hp_base == b:
                self.hp_len += 1
            else:
                self.hp_base, self.hp_len = b, 1
    
    def analyze(seq):
        gc = sum(1 for b in seq if b in {2, 3}) / len(seq) if seq else 0
        max_hp = 1
        hp = 1
        for i in range(1, len(seq)):
            if seq[i] == seq[i-1]:
                hp += 1
                max_hp = max(max_hp, hp)
            else:
                hp = 1
        return abs(gc - 0.5), max_hp
    
    def baseline_encode(latent, bits=8):
        bits_list = []
        for v in latent:
            q = int(max(0, min(1, v)) * (2**bits - 1))
            for i in range(bits - 1, -1, -1):
                bits_list.append((q >> i) & 1)
        while len(bits_list) % 2:
            bits_list.append(0)
        return [(bits_list[i] << 1) | bits_list[i+1] for i in range(0, len(bits_list), 2)]
    
    def scone_encode(latent, use_fsm, bits=8):
        fsm = FSM() if use_fsm else None
        symbols = []
        for v in latent:
            q = int(max(0, min(1, v)) * (2**bits - 1))
            for shift in range(bits - 2, -1, -2):
                symbols.append((q >> shift) & 0b11)
        
        seq = []
        total_bits = 0.0
        for sym in symbols:
            if fsm:
                mask = fsm.get_mask()
                allowed = [i for i in range(4) if mask[i]]
                base = sym if sym in allowed else allowed[sym % len(allowed)]
                total_bits += math.log2(len(allowed))
                fsm.update(base)
            else:
                base = sym
                total_bits += 2.0
            seq.append(base)
        
        bpn = total_bits / len(seq) if seq else 0
        return seq, bpn
    
    # Generate latents
    rng = random.Random(seed)
    latents = [[rng.random() for _ in range(latent_dim)] for _ in range(num_sequences)]
    
    results = {
        'Baseline': {'gc_dev': [], 'max_hp': [], 'bpn': []},
        'SCONE': {'gc_dev': [], 'max_hp': [], 'bpn': []},
        'SCONE-no-FSM': {'gc_dev': [], 'max_hp': [], 'bpn': []},
    }
    
    for latent in latents:
        # Baseline
        seq = baseline_encode(latent)
        gc_dev, max_hp = analyze(seq)
        results['Baseline']['gc_dev'].append(gc_dev)
        results['Baseline']['max_hp'].append(max_hp)
        results['Baseline']['bpn'].append(2.0)
        
        # SCONE (FSM enabled)
        seq, bpn = scone_encode(latent, use_fsm=True)
        gc_dev, max_hp = analyze(seq)
        results['SCONE']['gc_dev'].append(gc_dev)
        results['SCONE']['max_hp'].append(max_hp)
        results['SCONE']['bpn'].append(bpn)
        
        # SCONE-no-FSM
        seq, bpn = scone_encode(latent, use_fsm=False)
        gc_dev, max_hp = analyze(seq)
        results['SCONE-no-FSM']['gc_dev'].append(gc_dev)
        results['SCONE-no-FSM']['max_hp'].append(max_hp)
        results['SCONE-no-FSM']['bpn'].append(bpn)
    
    # Compute stats
    stats = {}
    for method, data in results.items():
        stats[method] = {
            'gc_dev_mean': statistics.mean(data['gc_dev']),
            'gc_dev_std': statistics.stdev(data['gc_dev']),
            'max_hp_mean': statistics.mean(data['max_hp']),
            'max_hp_std': statistics.stdev(data['max_hp']),
            'bpn_mean': statistics.mean(data['bpn']),
            'bpn_std': statistics.stdev(data['bpn']) if len(set(data['bpn'])) > 1 else 0,
        }
    
    return stats


def plot_figure_2(stats: dict, save_prefix: str = "fig2"):
    """
    Create publication-ready figures:
      Fig 2(a): GC deviation + Max homopolymer (grouped bar)
      Fig 2(b): Bits per nucleotide
    """
    methods = ['Baseline', 'SCONE', 'SCONE-no-FSM']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    # =========================================================================
    # Fig 2(a): GC Deviation + Max Homopolymer (side-by-side)
    # =========================================================================
    fig_a, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    x = np.arange(len(methods))
    width = 0.6
    
    # Left panel: GC Deviation
    gc_means = [stats[m]['gc_dev_mean'] for m in methods]
    gc_stds = [stats[m]['gc_dev_std'] for m in methods]
    
    bars1 = ax1.bar(x, gc_means, width, yerr=gc_stds, capsize=5, 
                    color=colors, edgecolor='black', linewidth=0.8,
                    error_kw={'linewidth': 1.5, 'capthick': 1.5})
    
    ax1.set_ylabel('GC Deviation from 50%', fontsize=11)
    ax1.set_xlabel('Method', fontsize=11)
    ax1.set_title('(a) GC Content Control', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods, fontsize=10)
    ax1.set_ylim(0, max(gc_means) * 1.6)
    ax1.axhline(y=0.05, color='red', linestyle='--', linewidth=1, alpha=0.7, label='5% threshold')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Add value labels
    for i, (mean, std) in enumerate(zip(gc_means, gc_stds)):
        ax1.text(i, mean + std + 0.003, f'{mean:.3f}', ha='center', va='bottom', fontsize=9)
    
    # Right panel: Max Homopolymer
    hp_means = [stats[m]['max_hp_mean'] for m in methods]
    hp_stds = [stats[m]['max_hp_std'] for m in methods]
    
    bars2 = ax2.bar(x, hp_means, width, yerr=hp_stds, capsize=5,
                    color=colors, edgecolor='black', linewidth=0.8,
                    error_kw={'linewidth': 1.5, 'capthick': 1.5})
    
    ax2.set_ylabel('Max Homopolymer Length', fontsize=11)
    ax2.set_xlabel('Method', fontsize=11)
    ax2.set_title('(b) Homopolymer Control', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(methods, fontsize=10)
    ax2.set_ylim(0, max(hp_means) * 1.5)
    ax2.axhline(y=3, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Target limit (≤3)')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    for i, (mean, std) in enumerate(zip(hp_means, hp_stds)):
        ax2.text(i, mean + std + 0.15, f'{mean:.2f}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    fig_a.savefig(f'{save_prefix}_a_constraints.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {save_prefix}_a_constraints.png")
    
    # =========================================================================
    # Fig 2(b): Bits per Nucleotide
    # =========================================================================
    fig_b, ax3 = plt.subplots(figsize=(6, 4))
    
    bpn_means = [stats[m]['bpn_mean'] for m in methods]
    bpn_stds = [stats[m]['bpn_std'] for m in methods]
    
    bars3 = ax3.bar(x, bpn_means, width, yerr=bpn_stds, capsize=5,
                    color=colors, edgecolor='black', linewidth=0.8,
                    error_kw={'linewidth': 1.5, 'capthick': 1.5})
    
    ax3.set_ylabel('Bits per Nucleotide', fontsize=11)
    ax3.set_xlabel('Method', fontsize=11)
    ax3.set_title('Encoding Efficiency', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(methods, fontsize=10)
    ax3.set_ylim(0, 2.5)
    ax3.axhline(y=2.0, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Shannon limit (2 bits/nt)')
    ax3.legend(loc='upper right', fontsize=9)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    for i, (mean, std) in enumerate(zip(bpn_means, bpn_stds)):
        label = f'{mean:.3f}' if std == 0 else f'{mean:.3f}±{std:.3f}'
        ax3.text(i, mean + std + 0.05, label, ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    fig_b.savefig(f'{save_prefix}_b_efficiency.png', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {save_prefix}_b_efficiency.png")
    
    plt.show()
    
    return fig_a, fig_b


def print_table(stats: dict):
    """Print results as formatted table."""
    methods = ['Baseline', 'SCONE', 'SCONE-no-FSM']
    
    print("\n" + "="*80)
    print("EXPERIMENTAL RESULTS (mean ± std, n=200)")
    print("="*80)
    print(f"{'Method':<15} {'GC Deviation':<20} {'Max Homopolymer':<20} {'Bits/nt':<15}")
    print("-"*80)
    
    for m in methods:
        s = stats[m]
        gc = f"{s['gc_dev_mean']:.4f}±{s['gc_dev_std']:.4f}"
        hp = f"{s['max_hp_mean']:.2f}±{s['max_hp_std']:.2f}"
        bpn = f"{s['bpn_mean']:.3f}±{s['bpn_std']:.3f}" if s['bpn_std'] > 0 else f"{s['bpn_mean']:.3f}"
        print(f"{m:<15} {gc:<20} {hp:<20} {bpn:<15}")
    
    print("="*80)
    
    # LaTeX
    print("\nLaTeX Table:")
    print("-"*80)
    print("\\begin{tabular}{lccc}")
    print("\\toprule")
    print("Method & GC Deviation & Max HP & Bits/nt \\\\")
    print("\\midrule")
    for m in methods:
        s = stats[m]
        print(f"{m} & {s['gc_dev_mean']:.4f}$\\pm${s['gc_dev_std']:.4f} & "
              f"{s['max_hp_mean']:.2f}$\\pm${s['max_hp_std']:.2f} & "
              f"{s['bpn_mean']:.3f} \\\\")
    print("\\bottomrule")
    print("\\end{tabular}")


def main():
    print("Running experiment...")
    stats = run_experiment(num_sequences=200, latent_dim=32, seed=42)
    
    print_table(stats)
    
    print("\nGenerating figures...")
    plot_figure_2(stats, save_prefix="fig2")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
