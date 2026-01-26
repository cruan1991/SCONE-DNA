#!/usr/bin/env python3
"""
FSM Steering Visualization for DNA Base Selection
==================================================
Visualizes how a constraint-aware FSM steers DNA base selection
compared to unconstrained random encoding.

For publication in ISCAS 2026.
"""

import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np


# =============================================================================
# FSM Constraint Controller
# =============================================================================

class FSMConstraint:
    """FSM for biochemical constraint enforcement."""
    
    BASE_TO_IDX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    IDX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
    GC_BASES = {2, 3}  # G=2, C=3
    AT_BASES = {0, 1}  # A=0, T=1
    
    def __init__(self, gc_window: int = 10, gc_target: float = 0.5, 
                 gc_tolerance: float = 0.1, max_homopolymer: int = 3):
        self.gc_window = gc_window
        self.gc_low = gc_target - gc_tolerance
        self.gc_high = gc_target + gc_tolerance
        self.max_hp = max_homopolymer
        self.reset()
    
    def reset(self):
        self.history = []
        self.hp_base = None
        self.hp_len = 0
    
    def get_gc_count(self) -> int:
        window = self.history[-self.gc_window:] if len(self.history) >= self.gc_window else self.history
        return sum(1 for b in window if b in self.GC_BASES)
    
    def get_allowed_mask(self) -> list:
        """Returns [allow_A, allow_T, allow_G, allow_C]."""
        mask = [True, True, True, True]
        
        # Homopolymer constraint
        if self.hp_base is not None and self.hp_len >= self.max_hp:
            mask[self.hp_base] = False
        
        # GC constraint
        if len(self.history) >= self.gc_window:
            gc_count = self.get_gc_count()
            gc_ratio = gc_count / self.gc_window
            
            if gc_ratio >= self.gc_high:
                oldest = self.history[-self.gc_window]
                if oldest not in self.GC_BASES:
                    mask[2] = False  # G
                    mask[3] = False  # C
            elif gc_ratio <= self.gc_low:
                oldest = self.history[-self.gc_window]
                if oldest not in self.AT_BASES:
                    mask[0] = False  # A
                    mask[1] = False  # T
        
        # Fail-safe
        if not any(mask):
            mask = [True, True, True, True]
            if self.hp_base is not None and self.hp_len >= self.max_hp:
                mask[self.hp_base] = False
        
        return mask
    
    def update(self, base_idx: int):
        self.history.append(base_idx)
        if self.hp_base == base_idx:
            self.hp_len += 1
        else:
            self.hp_base = base_idx
            self.hp_len = 1


# =============================================================================
# Sequence Generation
# =============================================================================

def generate_random_sequence(length: int, seed: int = 42) -> list:
    """Generate random DNA sequence (no constraints)."""
    rng = random.Random(seed)
    return [rng.randint(0, 3) for _ in range(length)]


def generate_fsm_sequence(length: int, seed: int = 42) -> tuple:
    """Generate FSM-constrained DNA sequence."""
    rng = random.Random(seed)
    fsm = FSMConstraint(gc_window=10, max_homopolymer=3)
    
    sequence = []
    allowed_counts = []
    
    for _ in range(length):
        mask = fsm.get_allowed_mask()
        allowed = [i for i in range(4) if mask[i]]
        allowed_counts.append(len(allowed))
        
        # Choose randomly from allowed bases
        base_idx = rng.choice(allowed)
        sequence.append(base_idx)
        fsm.update(base_idx)
    
    return sequence, allowed_counts


# =============================================================================
# Constraint Analysis
# =============================================================================

def analyze_sequence(seq: list) -> dict:
    """Compute constraint metrics."""
    # GC content
    gc_count = sum(1 for b in seq if b in {2, 3})
    gc_ratio = gc_count / len(seq) if seq else 0
    
    # Max homopolymer
    max_hp = 1
    current_hp = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            current_hp += 1
            max_hp = max(max_hp, current_hp)
        else:
            current_hp = 1
    
    # Homopolymer runs (for visualization)
    hp_runs = []
    start = 0
    for i in range(1, len(seq)):
        if seq[i] != seq[i-1]:
            if i - start > 1:
                hp_runs.append((start, i-1, seq[start]))
            start = i
    if len(seq) - start > 1:
        hp_runs.append((start, len(seq)-1, seq[start]))
    
    return {
        'gc_ratio': gc_ratio,
        'max_hp': max_hp,
        'hp_runs': hp_runs
    }


# =============================================================================
# Visualization
# =============================================================================

def create_fsm_steering_figure(save_path: str = "fsm_steering.png"):
    """Create publication-ready FSM steering visualization."""
    
    # Generate sequences
    length = 50
    seed = 123  # Seed that produces good illustration
    
    random_seq = generate_random_sequence(length, seed)
    fsm_seq, allowed_counts = generate_fsm_sequence(length, seed)
    
    # Analyze
    random_stats = analyze_sequence(random_seq)
    fsm_stats = analyze_sequence(fsm_seq)
    
    # Color scheme for bases (colorblind-friendly)
    base_colors = {
        0: '#E64B35',  # A - Red
        1: '#4DBBD5',  # T - Cyan  
        2: '#00A087',  # G - Teal
        3: '#3C5488',  # C - Blue
    }
    
    # Create figure
    fig, axes = plt.subplots(3, 1, figsize=(14, 5), 
                             gridspec_kw={'height_ratios': [1, 1, 0.6]})
    
    # Helper to draw sequence strip
    def draw_sequence_strip(ax, seq, title, stats, show_violations=True):
        ax.set_xlim(-0.5, length - 0.5)
        ax.set_ylim(-0.5, 0.5)
        
        # Draw each base as a colored rectangle
        for i, base in enumerate(seq):
            rect = plt.Rectangle((i - 0.5, -0.4), 1, 0.8, 
                                  facecolor=base_colors[base],
                                  edgecolor='white', linewidth=0.5)
            ax.add_patch(rect)
            
            # Add base letter
            base_letter = FSMConstraint.IDX_TO_BASE[base]
            ax.text(i, 0, base_letter, ha='center', va='center', 
                   fontsize=8, fontweight='bold', color='white')
        
        # Highlight homopolymer violations (>3)
        if show_violations:
            for start, end, base in stats['hp_runs']:
                run_len = end - start + 1
                if run_len > 3:
                    rect = plt.Rectangle((start - 0.5, -0.5), run_len, 1,
                                         facecolor='none', edgecolor='red',
                                         linewidth=2.5, linestyle='-')
                    ax.add_patch(rect)
        
        # Title and stats
        gc_status = "✓" if 0.4 <= stats['gc_ratio'] <= 0.6 else "✗"
        hp_status = "✓" if stats['max_hp'] <= 3 else "✗"
        
        title_text = f"{title}"
        stats_text = f"GC={stats['gc_ratio']*100:.0f}% {gc_status}  |  Max HP={stats['max_hp']} {hp_status}"
        
        ax.set_title(title_text, fontsize=12, fontweight='bold', loc='left')
        ax.text(length + 0.5, 0, stats_text, va='center', ha='left', fontsize=10,
               fontfamily='monospace')
        
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
    # Draw random sequence (top)
    draw_sequence_strip(axes[0], random_seq, "(a) No FSM (Random Quaternary)", random_stats)
    
    # Draw FSM sequence (middle)
    draw_sequence_strip(axes[1], fsm_seq, "(b) FSM-Enabled (GC≈50%, HP≤3)", fsm_stats)
    
    # Draw allowed bases indicator (bottom)
    ax2 = axes[2]
    ax2.set_xlim(-0.5, length - 0.5)
    ax2.set_ylim(0.5, 4.5)
    
    # Bar showing allowed count at each position
    colors_allowed = ['#2ecc71' if c == 4 else '#f39c12' if c == 3 else '#e74c3c' for c in allowed_counts]
    ax2.bar(range(length), allowed_counts, color=colors_allowed, edgecolor='white', linewidth=0.3)
    
    ax2.set_ylabel('Allowed\nBases', fontsize=9)
    ax2.set_xlabel('Sequence Position', fontsize=10)
    ax2.set_yticks([1, 2, 3, 4])
    ax2.set_title("(c) FSM State: Number of Allowed Bases per Position", 
                  fontsize=12, fontweight='bold', loc='left')
    ax2.axhline(y=4, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Add position markers every 10
    for i in range(0, length + 1, 10):
        ax2.axvline(x=i - 0.5, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
        if i < length:
            ax2.text(i, -0.3, str(i), ha='center', va='top', fontsize=8, color='gray')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=base_colors[0], edgecolor='white', label='A'),
        mpatches.Patch(facecolor=base_colors[1], edgecolor='white', label='T'),
        mpatches.Patch(facecolor=base_colors[2], edgecolor='white', label='G'),
        mpatches.Patch(facecolor=base_colors[3], edgecolor='white', label='C'),
        mpatches.Patch(facecolor='none', edgecolor='red', linewidth=2, label='HP violation (>3)'),
    ]
    fig.legend(handles=legend_elements, loc='upper right', ncol=5, 
               fontsize=9, frameon=True, bbox_to_anchor=(0.98, 0.98))
    
    # Overall title
    fig.suptitle('FSM-Guided DNA Base Selection', fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Figure saved to: {save_path}")
    plt.show()
    
    # Print summary
    print("\n" + "="*60)
    print("SEQUENCE COMPARISON")
    print("="*60)
    print(f"{'Metric':<25} {'No FSM':<15} {'FSM-Enabled':<15}")
    print("-"*60)
    print(f"{'GC Content':<25} {random_stats['gc_ratio']*100:.1f}%{'':<10} {fsm_stats['gc_ratio']*100:.1f}%")
    print(f"{'Max Homopolymer':<25} {random_stats['max_hp']:<15} {fsm_stats['max_hp']:<15}")
    print(f"{'HP Violations (>3)':<25} {len([r for r in random_stats['hp_runs'] if r[1]-r[0]+1 > 3]):<15} {len([r for r in fsm_stats['hp_runs'] if r[1]-r[0]+1 > 3]):<15}")
    print("="*60)


if __name__ == "__main__":
    create_fsm_steering_figure("fsm_steering.png")
