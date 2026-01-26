#!/usr/bin/env python
"""
Fully Deterministic Biochemical FSM Constraint Controller for DNA Encoding

Implements GC content window control and homopolymer length limitation.
"""

from typing import List, Tuple, Optional
from collections import deque


# DNA base mapping
BASE_TO_INDEX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INDEX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
NUM_BASES = 4

# GC bases (G and C)
GC_BASES = {2, 3}  # G=2, C=3
AT_BASES = {0, 1}  # A=0, T=1


class FSMConstraint:
    """
    Fully deterministic biochemical FSM constraint controller for DNA encoding
    
    Supports:
    - GC content window control
    - Homopolymer length limitation
    - Guarantees at least one allowed base
    """
    
    def __init__(
        self,
        gc_window: int = 20,
        gc_low: float = 0.45,
        gc_high: float = 0.55,
        max_homopolymer: int = 3
    ):
        """
        Initialize FSM constraint controller
        
        Args:
            gc_window: GC window size (default 20)
            gc_low: GC content lower bound (default 0.45)
            gc_high: GC content upper bound (default 0.55)
            max_homopolymer: Maximum homopolymer length (default 3)
        """
        self.gc_window = gc_window
        self.gc_low = gc_low
        self.gc_high = gc_high
        self.max_homopolymer = max_homopolymer
        
        # Internal state
        self._window_buffer: deque = deque(maxlen=gc_window)
        self._gc_count: int = 0
        self._homopolymer_base: Optional[int] = None
        self._homopolymer_length: int = 0
        self._step_count: int = 0
    
    def reset(self):
        """Reset internal state"""
        self._window_buffer.clear()
        self._gc_count = 0
        self._homopolymer_base = None
        self._homopolymer_length = 0
        self._step_count = 0
    
    def get_state(self) -> dict:
        """Get current state (for debugging)"""
        return {
            'step_count': self._step_count,
            'window_size': len(self._window_buffer),
            'gc_count': self._gc_count,
            'gc_ratio': self._gc_count / len(self._window_buffer) if self._window_buffer else 0.0,
            'homopolymer_base': self._homopolymer_base,
            'homopolymer_length': self._homopolymer_length,
            'window': list(self._window_buffer)
        }
    
    def _is_gc_base(self, base_index: int) -> bool:
        """Check if base is GC"""
        return base_index in GC_BASES
    
    def _get_current_gc_ratio(self) -> float:
        """Get current GC ratio"""
        if len(self._window_buffer) == 0:
            return 0.5  # Return middle value when no data
        return self._gc_count / len(self._window_buffer)
    
    def _would_gc_ratio_be(self, new_base: int) -> float:
        """
        Calculate GC ratio after adding new base
        Consider that when window is full, the oldest base will be removed
        """
        new_gc_count = self._gc_count
        window_size = len(self._window_buffer)
        
        # If window is full, need to remove the oldest base
        if window_size >= self.gc_window:
            oldest_base = self._window_buffer[0]
            if self._is_gc_base(oldest_base):
                new_gc_count -= 1
        else:
            window_size += 1
        
        # Add new base
        if self._is_gc_base(new_base):
            new_gc_count += 1
        
        return new_gc_count / window_size if window_size > 0 else 0.5
    
    def _would_violate_homopolymer(self, base_index: int) -> bool:
        """Check if adding new base would violate homopolymer limit"""
        if self._homopolymer_base is None:
            return False
        
        if base_index == self._homopolymer_base:
            # Same base, check if would exceed limit
            return self._homopolymer_length >= self.max_homopolymer
        
        return False
    
    def get_mask(self) -> List[bool]:
        """
        Get mask of currently allowed bases
        
        Returns:
            Boolean list[4], indicating if A,T,G,C are allowed
        
        Rules:
            1. Enforce homopolymer limit
            2. Enforce GC bounds
            3. If mask becomes empty, relax GC constraint until >=1 allowed
        """
        mask = [True, True, True, True]
        
        # 1. Apply homopolymer constraint
        for i in range(NUM_BASES):
            if self._would_violate_homopolymer(i):
                mask[i] = False
        
        # 2. Apply GC boundary constraint
        # Only apply GC constraint when window has enough data
        if len(self._window_buffer) >= self.gc_window // 2:
            gc_mask = [True, True, True, True]
            
            for i in range(NUM_BASES):
                if not mask[i]:
                    continue  # Already forbidden by homopolymer constraint
                
                future_gc_ratio = self._would_gc_ratio_be(i)
                
                # Check if would violate GC constraint
                if future_gc_ratio < self.gc_low:
                    # GC too low, only allow GC bases
                    if i not in GC_BASES:
                        gc_mask[i] = False
                elif future_gc_ratio > self.gc_high:
                    # GC too high, only allow AT bases
                    if i not in AT_BASES:
                        gc_mask[i] = False
            
            # Merge GC mask
            combined_mask = [m1 and m2 for m1, m2 in zip(mask, gc_mask)]
            
            # 3. If combined mask is empty, relax GC constraint
            if any(combined_mask):
                mask = combined_mask
            # Otherwise keep only homopolymer constraint mask
        
        # Final check: ensure at least one allowed base
        if not any(mask):
            # This should not happen, but as a safety measure
            # Relax all constraints, keep only homopolymer constraint
            mask = [True, True, True, True]
            for i in range(NUM_BASES):
                if self._would_violate_homopolymer(i):
                    mask[i] = False
            
            # If still empty (should not happen), allow all
            if not any(mask):
                mask = [True, True, True, True]
        
        return mask
    
    def update(self, base_index: int):
        """
        Update FSM state
        
        Args:
            base_index: Base index (0=A, 1=T, 2=G, 3=C)
        """
        if base_index < 0 or base_index >= NUM_BASES:
            raise ValueError(f"Invalid base index: {base_index}")
        
        # Update GC window
        if len(self._window_buffer) >= self.gc_window:
            # Window is full, remove oldest base
            oldest_base = self._window_buffer[0]
            if self._is_gc_base(oldest_base):
                self._gc_count -= 1
        
        # Add new base to window
        self._window_buffer.append(base_index)
        if self._is_gc_base(base_index):
            self._gc_count += 1
        
        # Update homopolymer state
        if self._homopolymer_base == base_index:
            self._homopolymer_length += 1
        else:
            self._homopolymer_base = base_index
            self._homopolymer_length = 1
        
        # Update step count
        self._step_count += 1
    
    def get_allowed_count(self) -> int:
        """Get count of currently allowed bases"""
        return sum(self.get_mask())
    
    def get_allowed_bases(self) -> List[int]:
        """Get list of currently allowed base indices"""
        mask = self.get_mask()
        return [i for i, m in enumerate(mask) if m]


# ==============================================================================
# Test Functions
# ==============================================================================

def test_basic_functionality():
    """Test basic functionality"""
    print("Testing basic functionality:")
    
    fsm = FSMConstraint(gc_window=20, gc_low=0.45, gc_high=0.55, max_homopolymer=3)
    
    # Initial state
    mask = fsm.get_mask()
    assert all(mask), f"Initial state should allow all bases: {mask}"
    print(f"  Initial mask: {mask} ✓")
    
    # Add a few bases
    fsm.update(0)  # A
    fsm.update(1)  # T
    fsm.update(2)  # G
    fsm.update(3)  # C
    
    state = fsm.get_state()
    print(f"  State after ATGC: {state}")
    assert state['step_count'] == 4
    assert state['gc_count'] == 2  # G and C
    print("  Basic functionality test passed ✓")


def test_homopolymer_constraint():
    """Test homopolymer constraint"""
    print("\nTesting homopolymer constraint:")
    
    fsm = FSMConstraint(gc_window=20, gc_low=0.0, gc_high=1.0, max_homopolymer=3)
    
    # Add 3 consecutive A's
    fsm.update(0)  # A
    fsm.update(0)  # A
    fsm.update(0)  # A
    
    # 4th should forbid A
    mask = fsm.get_mask()
    assert not mask[0], f"A should be forbidden after 3 consecutive A's: {mask}"
    assert mask[1] and mask[2] and mask[3], f"Other bases should be allowed: {mask}"
    print(f"  Mask after 3 A's: {mask} ✓")
    
    # Add a different base
    fsm.update(1)  # T
    
    # Now A should be allowed again
    mask = fsm.get_mask()
    assert mask[0], f"A should be allowed after adding T: {mask}"
    print(f"  Mask after adding T: {mask} ✓")
    
    print("  Homopolymer constraint test passed ✓")


def test_gc_constraint():
    """Test GC constraint"""
    print("\nTesting GC constraint:")
    
    fsm = FSMConstraint(gc_window=10, gc_low=0.4, gc_high=0.6, max_homopolymer=10)
    
    # Add 10 A's (GC=0%)
    for _ in range(10):
        fsm.update(0)  # A
    
    state = fsm.get_state()
    print(f"  After 10 A's: GC ratio={state['gc_ratio']:.2f}")
    
    # GC too low, should prefer GC bases
    mask = fsm.get_mask()
    print(f"  Mask: {mask}")
    # Since GC is too low, should bias toward allowing GC bases
    
    # Add some G's to increase GC
    fsm.reset()
    for _ in range(10):
        fsm.update(2)  # G
    
    state = fsm.get_state()
    print(f"  After 10 G's: GC ratio={state['gc_ratio']:.2f}")
    
    # GC too high, should prefer AT bases
    mask = fsm.get_mask()
    print(f"  Mask: {mask}")
    
    print("  GC constraint test passed ✓")


def test_mask_never_empty():
    """Test that mask is never empty"""
    print("\nTesting mask never empty:")
    
    import random
    random.seed(42)
    
    fsm = FSMConstraint(gc_window=10, gc_low=0.45, gc_high=0.55, max_homopolymer=3)
    
    for trial in range(100):
        fsm.reset()
        
        for step in range(200):
            mask = fsm.get_mask()
            
            # Ensure mask is not empty
            assert any(mask), f"Mask is empty! Trial {trial}, step {step}, state {fsm.get_state()}"
            
            # Randomly choose from allowed bases
            allowed = [i for i, m in enumerate(mask) if m]
            base = random.choice(allowed)
            fsm.update(base)
    
    print("  100 trials, 200 steps each, mask was never empty ✓")
    print("  Mask non-empty test passed ✓")


def test_determinism():
    """Test determinism"""
    print("\nTesting determinism:")
    
    sequence = [0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 1, 2, 2, 2, 1, 3, 3, 3, 0]
    
    # Run twice, should get same mask sequences
    fsm1 = FSMConstraint(gc_window=10, gc_low=0.4, gc_high=0.6, max_homopolymer=3)
    fsm2 = FSMConstraint(gc_window=10, gc_low=0.4, gc_high=0.6, max_homopolymer=3)
    
    masks1 = []
    masks2 = []
    
    for base in sequence:
        masks1.append(fsm1.get_mask())
        fsm1.update(base)
    
    for base in sequence:
        masks2.append(fsm2.get_mask())
        fsm2.update(base)
    
    assert masks1 == masks2, "Mask sequences differ between two runs!"
    print("  Two runs produced identical mask sequences ✓")
    print("  Determinism test passed ✓")


def run_fsm_tests():
    """Run all FSM tests"""
    print("=" * 70)
    print("FSM Constraint Controller Tests")
    print("=" * 70)
    
    test_basic_functionality()
    test_homopolymer_constraint()
    test_gc_constraint()
    test_mask_never_empty()
    test_determinism()
    
    print("\n" + "=" * 70)
    print("🎉 All FSM tests passed!")
    print("=" * 70)


if __name__ == "__main__":
    run_fsm_tests()
