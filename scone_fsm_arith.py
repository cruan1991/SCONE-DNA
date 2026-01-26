#!/usr/bin/env python
"""
SCONE FSM Arithmetic Encoder

Combines FSM constraint controller and masked arithmetic encoder for DNA encoding.
"""

from typing import List, Tuple, Optional
from fsm_constraint import FSMConstraint, INDEX_TO_BASE, BASE_TO_INDEX, NUM_BASES
from masked_arithmetic_codec import (
    apply_mask_and_renormalize,
    probs_to_freqs_with_mask,
    freqs_to_cumfreq,
    get_allowed_indices,
    symbol_to_masked_index,
    masked_index_to_symbol,
    MaskedArithmeticEncoder,
    MaskedArithmeticDecoder,
    TOT
)


# ==============================================================================
# Main Encode/Decode Functions
# ==============================================================================

def encode_fsm(
    latent_symbols: List[int],
    base_probs: List[float],
    gc_window: int = 20,
    gc_low: float = 0.45,
    gc_high: float = 0.55,
    max_homopolymer: int = 3
) -> Tuple[List[int], str]:
    """
    Encode latent symbol sequence using FSM constraints
    
    Args:
        latent_symbols: Latent symbol list (each symbol ∈ {0,1,2,3} representing base index to encode)
        base_probs: Base probability distribution [p_A, p_T, p_G, p_C]
        gc_window: GC window size
        gc_low: GC content lower bound
        gc_high: GC content upper bound
        max_homopolymer: Maximum homopolymer length
    
    Returns:
        (bitstream, dna_string): Bitstream and DNA sequence string
    """
    # Create FSM constraint controller
    fsm = FSMConstraint(
        gc_window=gc_window,
        gc_low=gc_low,
        gc_high=gc_high,
        max_homopolymer=max_homopolymer
    )
    
    # Create arithmetic encoder
    encoder = MaskedArithmeticEncoder()
    
    # DNA sequence
    dna_bases = []
    
    # Encode each symbol
    for symbol in latent_symbols:
        if symbol < 0 or symbol >= NUM_BASES:
            raise ValueError(f"Invalid symbol: {symbol}, must be in range [0, {NUM_BASES-1}]")
        
        # Get current mask
        mask = fsm.get_mask()
        
        # Check if symbol is allowed
        if not mask[symbol]:
            raise ValueError(f"Symbol {symbol} is forbidden in current FSM state")
        
        # Always use arithmetic coding (including allowed_count==1 case, for symmetry with decoder)
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # Convert original symbol index to masked index
        masked_idx = symbol_to_masked_index(symbol, mask)
        
        # Encode
        encoder.encode_symbol(masked_idx, cumfreq, TOT)
        
        # Update FSM state
        fsm.update(symbol)
        
        # Add to DNA sequence
        dna_bases.append(INDEX_TO_BASE[symbol])
    
    # Encode EOS
    # Use current FSM state's mask
    # Note: Always encode EOS, even if allowed_count==1, because decoder needs to know when to stop
    mask = fsm.get_mask()
    masked_probs = apply_mask_and_renormalize(base_probs, mask)
    freqs = probs_to_freqs_with_mask(masked_probs, mask)
    cumfreq = freqs_to_cumfreq(freqs)
    eos_idx = len(freqs) - 1
    encoder.encode_symbol(eos_idx, cumfreq, TOT)
    
    # Finish encoding
    encoder.finish()
    
    # Return bitstream and DNA sequence
    bitstream = encoder.get_bitstream()
    dna_string = ''.join(dna_bases)
    
    return bitstream, dna_string


def decode_fsm(
    bits: List[int],
    base_probs: List[float],
    gc_window: int = 20,
    gc_low: float = 0.45,
    gc_high: float = 0.55,
    max_homopolymer: int = 3,
    max_symbols: Optional[int] = None
) -> Tuple[List[int], str]:
    """
    Decode bitstream using FSM constraints
    
    Args:
        bits: Bitstream
        base_probs: Base probability distribution
        gc_window: GC window size
        gc_low: GC content lower bound
        gc_high: GC content upper bound
        max_homopolymer: Maximum homopolymer length
        max_symbols: Maximum number of symbols to decode (safety limit)
    
    Returns:
        (symbols, dna_string): Decoded symbol list and DNA sequence string
    """
    if max_symbols is None:
        max_symbols = 100000  # Default maximum limit
    
    # Create FSM constraint controller
    fsm = FSMConstraint(
        gc_window=gc_window,
        gc_low=gc_low,
        gc_high=gc_high,
        max_homopolymer=max_homopolymer
    )
    
    # Create arithmetic decoder
    decoder = MaskedArithmeticDecoder()
    decoder.initialize(bits)
    
    # Decode results
    decoded_symbols = []
    dna_bases = []
    
    # Decode loop
    for _ in range(max_symbols):
        # Get current mask
        mask = fsm.get_mask()
        
        # Always use arithmetic decoding (even if allowed_count==1, because need to check EOS)
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # Decode
        masked_idx = decoder.decode_symbol(cumfreq, TOT)
        
        # EOS check
        eos_idx = len(freqs) - 1
        if masked_idx == eos_idx:
            break
        
        # Convert masked index to original symbol index
        base = masked_index_to_symbol(masked_idx, mask)
        if base == -1:
            break  # EOS
        
        # Update FSM state
        fsm.update(base)
        
        # Add to results
        decoded_symbols.append(base)
        dna_bases.append(INDEX_TO_BASE[base])
    
    # Return decoded results
    dna_string = ''.join(dna_bases)
    
    return decoded_symbols, dna_string


# ==============================================================================
# Helper Functions
# ==============================================================================

def calculate_gc_content(dna_string: str) -> float:
    """Calculate GC content of DNA sequence"""
    if not dna_string:
        return 0.0
    gc_count = sum(1 for base in dna_string if base in 'GC')
    return gc_count / len(dna_string)


def calculate_max_homopolymer(dna_string: str) -> int:
    """Calculate maximum homopolymer length in DNA sequence"""
    if not dna_string:
        return 0
    
    max_run = 1
    current_run = 1
    
    for i in range(1, len(dna_string)):
        if dna_string[i] == dna_string[i-1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    
    return max_run


def calculate_windowed_gc(dna_string: str, window_size: int) -> Tuple[float, float]:
    """
    Calculate windowed GC content range of DNA sequence
    
    Returns:
        (min_gc, max_gc): Minimum and maximum window GC content
    """
    if len(dna_string) < window_size:
        gc = calculate_gc_content(dna_string)
        return gc, gc
    
    min_gc = 1.0
    max_gc = 0.0
    
    for i in range(len(dna_string) - window_size + 1):
        window = dna_string[i:i+window_size]
        gc = calculate_gc_content(window)
        min_gc = min(min_gc, gc)
        max_gc = max(max_gc, gc)
    
    return min_gc, max_gc


# ==============================================================================
# Test Functions
# ==============================================================================

def test_basic_encode_decode():
    """Test basic encode/decode"""
    print("Testing basic encode/decode:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Test simple sequence
    symbols = [0, 1, 2, 3, 0, 1, 2, 3]  # ATGCATGC
    
    bits, dna = encode_fsm(symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=len(symbols)+10)
    
    print(f"  Input symbols: {symbols}")
    print(f"  Encoded DNA: {dna}")
    print(f"  Bitstream length: {len(bits)}")
    print(f"  Decoded symbols: {decoded}")
    print(f"  Decoded DNA: {decoded_dna}")
    
    assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
    assert decoded_dna == dna, f"DNA mismatch: {decoded_dna} != {dna}"
    print("  Basic encode/decode test passed ✓")


def test_gc_control():
    """Test GC control"""
    print("\nTesting GC control:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    gc_window = 20
    gc_low = 0.40
    gc_high = 0.60
    
    # Generate random sequence (choosing from allowed bases)
    fsm = FSMConstraint(gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10)
    
    symbols = []
    for _ in range(100):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # Encode
    bits, dna = encode_fsm(symbols, base_probs, gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10)
    
    # Check GC content
    min_gc, max_gc = calculate_windowed_gc(dna, gc_window)
    overall_gc = calculate_gc_content(dna)
    
    print(f"  Sequence length: {len(dna)}")
    print(f"  Overall GC content: {overall_gc:.2%}")
    print(f"  Window GC range: [{min_gc:.2%}, {max_gc:.2%}]")
    print(f"  Target range: [{gc_low:.2%}, {gc_high:.2%}]")
    
    # For long enough sequences, window GC should be within range
    if len(dna) >= gc_window:
        # Allow some boundary tolerance
        assert min_gc >= gc_low - 0.1, f"Min GC too low: {min_gc:.2%}"
        assert max_gc <= gc_high + 0.1, f"Max GC too high: {max_gc:.2%}"
    
    # Decode
    decoded, decoded_dna = decode_fsm(bits, base_probs, gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"Decode failed"
    print("  GC control test passed ✓")


def test_homopolymer_control():
    """Test homopolymer control"""
    print("\nTesting homopolymer control:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    max_homopolymer = 3
    
    # Generate random sequence (choosing from allowed bases)
    fsm = FSMConstraint(gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer)
    
    symbols = []
    for _ in range(200):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # Encode
    bits, dna = encode_fsm(symbols, base_probs, gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer)
    
    # Check homopolymer length
    max_run = calculate_max_homopolymer(dna)
    
    print(f"  Sequence length: {len(dna)}")
    print(f"  Max homopolymer length: {max_run}")
    print(f"  Target limit: {max_homopolymer}")
    
    assert max_run <= max_homopolymer, f"Homopolymer exceeds limit: {max_run} > {max_homopolymer}"
    
    # Decode
    decoded, decoded_dna = decode_fsm(bits, base_probs, gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"Decode failed"
    print("  Homopolymer control test passed ✓")


def test_reversibility():
    """Test reversibility"""
    print("\nTesting reversibility:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    passed = 0
    failed = 0
    
    # Test different parameter combinations
    test_configs = [
        {'gc_window': 20, 'gc_low': 0.45, 'gc_high': 0.55, 'max_homopolymer': 3},
        {'gc_window': 10, 'gc_low': 0.40, 'gc_high': 0.60, 'max_homopolymer': 2},
        {'gc_window': 30, 'gc_low': 0.35, 'gc_high': 0.65, 'max_homopolymer': 4},
    ]
    
    for config in test_configs:
        for trial in range(10):
            # Generate random sequence
            fsm = FSMConstraint(**config)
            
            seq_len = random.randint(10, 100)
            symbols = []
            
            for _ in range(seq_len):
                mask = fsm.get_mask()
                allowed = [i for i, m in enumerate(mask) if m]
                symbol = random.choice(allowed)
                symbols.append(symbol)
                fsm.update(symbol)
            
            # Encode
            bits, dna = encode_fsm(symbols, base_probs, **config)
            
            # Decode
            decoded, decoded_dna = decode_fsm(bits, base_probs, **config, max_symbols=len(symbols)+10)
            
            if decoded == symbols:
                passed += 1
            else:
                failed += 1
                if failed <= 3:
                    print(f"  Failed case:")
                    print(f"    Config: {config}")
                    print(f"    Input length: {len(symbols)}")
                    print(f"    Decoded length: {len(decoded)}")
    
    print(f"  Passed: {passed}, Failed: {failed}")
    assert failed == 0, f"Reversibility test failed: {failed} cases"
    print("  Reversibility test passed ✓")


def test_random_1000_steps():
    """Test random 1000-step sequence"""
    print("\nTesting random 1000-step sequence:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    config = {'gc_window': 20, 'gc_low': 0.45, 'gc_high': 0.55, 'max_homopolymer': 3}
    
    # Generate 1000-step sequence
    fsm = FSMConstraint(**config)
    
    symbols = []
    for _ in range(1000):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # Encode
    bits, dna = encode_fsm(symbols, base_probs, **config)
    
    print(f"  Sequence length: {len(dna)}")
    print(f"  Bitstream length: {len(bits)}")
    print(f"  Bits per base: {len(bits)/len(dna):.3f}")
    
    # Check constraints
    gc = calculate_gc_content(dna)
    max_run = calculate_max_homopolymer(dna)
    min_gc, max_gc = calculate_windowed_gc(dna, config['gc_window'])
    
    print(f"  Overall GC content: {gc:.2%}")
    print(f"  Window GC range: [{min_gc:.2%}, {max_gc:.2%}]")
    print(f"  Max homopolymer: {max_run}")
    
    assert max_run <= config['max_homopolymer'], f"Homopolymer exceeds limit: {max_run}"
    
    # Decode
    decoded, decoded_dna = decode_fsm(bits, base_probs, **config, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"Decode failed: length {len(decoded)} vs {len(symbols)}"
    assert decoded_dna == dna, f"DNA mismatch"
    
    print("  Random 1000-step test passed ✓")


def test_edge_cases():
    """Test edge cases"""
    print("\nTesting edge cases:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Test empty sequence
    symbols = []
    bits, dna = encode_fsm(symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=10)
    assert decoded == symbols, f"Empty sequence decode failed"
    print("  Empty sequence ✓")
    
    # Test single symbol
    for s in range(4):
        symbols = [s]
        bits, dna = encode_fsm(symbols, base_probs)
        decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=10)
        assert decoded == symbols, f"Single symbol {s} decode failed"
    print("  Single symbol ✓")
    
    # Test uneven probabilities
    base_probs = [0.5, 0.3, 0.15, 0.05]
    symbols = [0, 0, 1, 2, 3, 0, 1, 0]
    
    # Generate valid sequence
    fsm = FSMConstraint()
    valid_symbols = []
    for s in symbols:
        mask = fsm.get_mask()
        if mask[s]:
            valid_symbols.append(s)
            fsm.update(s)
    
    bits, dna = encode_fsm(valid_symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=len(valid_symbols)+10)
    assert decoded == valid_symbols, f"Uneven probabilities decode failed"
    print("  Uneven probabilities ✓")
    
    print("  Edge cases test passed ✓")


def run_all_tests():
    """Run all tests"""
    print("=" * 70)
    print("SCONE FSM Arithmetic Encoder Tests")
    print("=" * 70)
    
    test_basic_encode_decode()
    test_gc_control()
    test_homopolymer_control()
    test_reversibility()
    test_random_1000_steps()
    test_edge_cases()
    
    print("\n" + "=" * 70)
    print("🎉 All tests passed!")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
