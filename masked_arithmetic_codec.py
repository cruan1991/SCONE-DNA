#!/usr/bin/env python
"""
DNA Base Arithmetic Encoder with Static Masking Layer

Adds masking support on top of the standard arithmetic encoder for FSM-constrained DNA encoding.
Supports dynamic masking where allowed set size can be 4, 3, 2, or 1.
"""

from typing import List, Tuple, Optional


# ==============================================================================
# Constants
# ==============================================================================

# DNA base mapping
BASE_TO_INDEX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INDEX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
NUM_BASES = 4

# Arithmetic coding constants
HALF = 0x80000000  # 2^31
QUARTER = 0x40000000  # 2^30
THREE_QUARTER = 0xC0000000  # 3 * 2^30
LOW_INIT = 0
HIGH_INIT = 0xFFFFFFFF  # 2^32 - 1
TOT = 1 << 15  # Total cumulative frequency (32768)


# ==============================================================================
# Mask and Probability Processing Functions
# ==============================================================================

def apply_mask_and_renormalize(probs: List[float], mask: List[bool]) -> List[float]:
    """
    Apply mask and renormalize probabilities
    
    Args:
        probs: List of 4 floats representing probabilities for A,T,G,C
        mask: List of 4 booleans, True for allowed, False for forbidden
    
    Returns:
        Renormalized probability list (length 4, masked positions are 0)
    
    Raises:
        ValueError: If all mask[i] are False
    """
    if len(probs) != 4 or len(mask) != 4:
        raise ValueError("probs and mask length must be 4")
    
    # Check at least one base is allowed
    if not any(mask):
        raise ValueError("At least one base must be allowed (mask cannot be all False)")
    
    # Apply mask
    masked_probs = [p if m else 0.0 for p, m in zip(probs, mask)]
    
    # Calculate total
    total = sum(masked_probs)
    
    # Renormalize
    if total > 0:
        normalized_probs = [p / total for p in masked_probs]
    else:
        # If total is 0 (all allowed bases have 0 probability), distribute uniformly
        allowed_count = sum(mask)
        normalized_probs = [1.0 / allowed_count if m else 0.0 for m in mask]
    
    return normalized_probs


def probs_to_freqs_with_mask(
    probs: List[float],
    mask: List[bool],
    tot: int = TOT,
    eos_freq: int = 1
) -> List[int]:
    """
    Convert probabilities to frequencies (for arithmetic coding)
    
    Args:
        probs: Renormalized probability list (from apply_mask_and_renormalize)
        mask: Mask list
        tot: Total cumulative frequency
        eos_freq: Frequency for EOS symbol
    
    Returns:
        Frequency list, length is allowed_bases + 1 (last one is EOS)
    """
    # Calculate number of allowed bases
    allowed_count = sum(mask)
    
    if allowed_count == 0:
        raise ValueError("At least one base must be allowed")
    
    # Calculate available frequency for bases
    available_freq = tot - eos_freq
    
    # Only allocate frequency for allowed bases
    freqs = []
    for i, (p, m) in enumerate(zip(probs, mask)):
        if m:
            freq = max(1, int(p * available_freq))
            freqs.append(freq)
    
    # Adjust frequencies to ensure total is correct
    total_freq = sum(freqs)
    target_freq = available_freq
    
    if total_freq != target_freq and len(freqs) > 0:
        diff = target_freq - total_freq
        # Add difference to largest frequency
        max_idx = freqs.index(max(freqs))
        freqs[max_idx] += diff
    
    # Add EOS frequency
    freqs.append(eos_freq)
    
    return freqs


def freqs_to_cumfreq(freqs: List[int]) -> List[int]:
    """Convert frequencies to cumulative frequencies"""
    cumfreq = [0]
    for f in freqs:
        cumfreq.append(cumfreq[-1] + f)
    return cumfreq


def get_allowed_indices(mask: List[bool]) -> List[int]:
    """Get list of allowed base indices"""
    return [i for i, m in enumerate(mask) if m]


def symbol_to_masked_index(symbol: int, mask: List[bool]) -> int:
    """Convert original symbol index to masked index"""
    if not mask[symbol]:
        raise ValueError(f"Symbol {symbol} is forbidden by mask")
    
    allowed = get_allowed_indices(mask)
    return allowed.index(symbol)


def masked_index_to_symbol(masked_idx: int, mask: List[bool]) -> int:
    """Convert masked index to original symbol index"""
    allowed = get_allowed_indices(mask)
    if masked_idx >= len(allowed):
        return -1  # EOS
    return allowed[masked_idx]


# ==============================================================================
# Masked Arithmetic Encoder
# ==============================================================================

class MaskedArithmeticEncoder:
    """Arithmetic encoder with masking support"""
    
    def __init__(self):
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.pending_bits = 0
        self.bitstream = []
    
    def reset(self):
        """Reset encoder state"""
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.pending_bits = 0
        self.bitstream = []
    
    def _output_bit(self, bit: int):
        """Output one bit and handle pending bits"""
        self.bitstream.append(bit)
        for _ in range(self.pending_bits):
            self.bitstream.append(1 - bit)
        self.pending_bits = 0
    
    def _renormalize(self):
        """Renormalize"""
        while True:
            if self.high < HALF:
                self._output_bit(0)
                self.low = self.low * 2
                self.high = self.high * 2 + 1
            elif self.low >= HALF:
                self._output_bit(1)
                self.low = (self.low - HALF) * 2
                self.high = (self.high - HALF) * 2 + 1
            elif self.low >= QUARTER and self.high < THREE_QUARTER:
                self.pending_bits += 1
                self.low = (self.low - QUARTER) * 2
                self.high = (self.high - QUARTER) * 2 + 1
            else:
                break
    
    def encode_symbol(self, symbol_idx: int, cumfreq: List[int], tot: int):
        """Encode one symbol"""
        range_size = self.high - self.low + 1
        
        symbol_low = cumfreq[symbol_idx]
        symbol_high = cumfreq[symbol_idx + 1]
        
        new_low = self.low + (range_size * symbol_low) // tot
        new_high = self.low + (range_size * symbol_high) // tot - 1
        
        if new_low > new_high:
            new_high = new_low
        
        self.low = new_low
        self.high = new_high
        
        self._renormalize()
    
    def finish(self):
        """Finish encoding, output final bits"""
        self.pending_bits += 1
        if self.low < QUARTER:
            self._output_bit(0)
        else:
            self._output_bit(1)
    
    def get_bitstream(self) -> List[int]:
        """Get bitstream"""
        return self.bitstream


class MaskedArithmeticDecoder:
    """Arithmetic decoder with masking support"""
    
    def __init__(self):
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.code = 0
        self.bitstream = []
        self.bit_idx = 0
    
    def reset(self):
        """Reset decoder state"""
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.code = 0
        self.bit_idx = 0
    
    def _read_bit(self) -> int:
        """Read one bit"""
        if self.bit_idx < len(self.bitstream):
            bit = self.bitstream[self.bit_idx]
            self.bit_idx += 1
            return bit
        return 0
    
    def initialize(self, bitstream: List[int]):
        """Initialize decoder"""
        self.bitstream = bitstream
        self.bit_idx = 0
        self.low = LOW_INIT
        self.high = HIGH_INIT
        
        # Read first 32 bits
        self.code = 0
        for _ in range(32):
            self.code = self.code * 2 + self._read_bit()
    
    def _renormalize(self):
        """Renormalize"""
        while True:
            if self.high < HALF:
                self.low = self.low * 2
                self.high = self.high * 2 + 1
                self.code = self.code * 2 + self._read_bit()
            elif self.low >= HALF:
                self.low = (self.low - HALF) * 2
                self.high = (self.high - HALF) * 2 + 1
                self.code = (self.code - HALF) * 2 + self._read_bit()
            elif self.low >= QUARTER and self.high < THREE_QUARTER:
                self.low = (self.low - QUARTER) * 2
                self.high = (self.high - QUARTER) * 2 + 1
                self.code = (self.code - QUARTER) * 2 + self._read_bit()
            else:
                break
    
    def decode_symbol(self, cumfreq: List[int], tot: int) -> int:
        """Decode one symbol"""
        range_size = self.high - self.low + 1
        
        # Find which symbol's interval code falls into
        num_symbols = len(cumfreq) - 1
        symbol = num_symbols - 1  # Default to last one (EOS)
        
        for i in range(num_symbols - 1, -1, -1):
            symbol_low_bound = self.low + (range_size * cumfreq[i]) // tot
            symbol_high_bound = self.low + (range_size * cumfreq[i + 1]) // tot - 1
            
            if symbol_low_bound > symbol_high_bound:
                symbol_high_bound = symbol_low_bound
            
            if symbol_low_bound <= self.code <= symbol_high_bound:
                symbol = i
                break
        
        # Update interval
        symbol_low = cumfreq[symbol]
        symbol_high = cumfreq[symbol + 1]
        
        new_low = self.low + (range_size * symbol_low) // tot
        new_high = self.low + (range_size * symbol_high) // tot - 1
        
        if new_low > new_high:
            new_high = new_low
        
        self.low = new_low
        self.high = new_high
        
        self._renormalize()
        
        return symbol


# ==============================================================================
# Main Encode/Decode Functions
# ==============================================================================

def encode_with_mask(
    symbols: List[int],
    masks: List[List[bool]],
    base_probs: List[float]
) -> List[int]:
    """
    Encode with masking
    
    Args:
        symbols: Symbol list, each symbol ∈ {0,1,2,3} representing A,T,G,C
        masks: Mask list, same length as symbols
        base_probs: Base probability distribution [p_A, p_T, p_G, p_C]
    
    Returns:
        Bitstream (list of integers)
    """
    if len(symbols) != len(masks):
        raise ValueError("symbols and masks must have same length")
    
    if len(symbols) == 0:
        # Empty sequence, only encode EOS
        encoder = MaskedArithmeticEncoder()
        # Use default mask (all allowed)
        default_mask = [True] * 4
        masked_probs = apply_mask_and_renormalize(base_probs, default_mask)
        freqs = probs_to_freqs_with_mask(masked_probs, default_mask)
        cumfreq = freqs_to_cumfreq(freqs)
        # EOS index is the last one
        eos_idx = len(freqs) - 1
        encoder.encode_symbol(eos_idx, cumfreq, TOT)
        encoder.finish()
        return encoder.get_bitstream()
    
    encoder = MaskedArithmeticEncoder()
    
    # Encode each symbol
    for i, (symbol, mask) in enumerate(zip(symbols, masks)):
        if not mask[symbol]:
            raise ValueError(f"Symbol {symbol} at position {i} is forbidden by mask")
        
        allowed_count = sum(mask)
        
        if allowed_count == 1:
            # Only one allowed base, no need to encode
            # Decoder knows the mask, can determine directly
            continue
        
        # Apply mask and renormalize
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # Convert original symbol index to masked index
        masked_idx = symbol_to_masked_index(symbol, mask)
        
        # Encode
        encoder.encode_symbol(masked_idx, cumfreq, TOT)
    
    # Encode EOS
    # Use last mask (or default mask)
    last_mask = masks[-1] if masks else [True] * 4
    allowed_count = sum(last_mask)
    
    if allowed_count > 1:
        masked_probs = apply_mask_and_renormalize(base_probs, last_mask)
        freqs = probs_to_freqs_with_mask(masked_probs, last_mask)
        cumfreq = freqs_to_cumfreq(freqs)
        eos_idx = len(freqs) - 1
        encoder.encode_symbol(eos_idx, cumfreq, TOT)
    
    encoder.finish()
    return encoder.get_bitstream()


def decode_with_mask(
    bits: List[int],
    masks: List[List[bool]],
    base_probs: List[float],
    max_symbols: Optional[int] = None
) -> List[int]:
    """
    Decode with masking
    
    Args:
        bits: Bitstream
        masks: Mask list
        base_probs: Base probability distribution
        max_symbols: Maximum number of symbols to decode (optional, for safety)
    
    Returns:
        Decoded symbol list
    """
    if max_symbols is None:
        max_symbols = len(masks) + 1  # Allow slightly more
    
    if len(masks) == 0:
        # Empty masks, decode EOS
        decoder = MaskedArithmeticDecoder()
        decoder.initialize(bits)
        
        default_mask = [True] * 4
        masked_probs = apply_mask_and_renormalize(base_probs, default_mask)
        freqs = probs_to_freqs_with_mask(masked_probs, default_mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        symbol = decoder.decode_symbol(cumfreq, TOT)
        eos_idx = len(freqs) - 1
        
        if symbol == eos_idx:
            return []
        else:
            raise ValueError("Decode error: expected EOS but got other symbol")
    
    decoder = MaskedArithmeticDecoder()
    decoder.initialize(bits)
    
    decoded_symbols = []
    
    for i in range(max_symbols):
        if i >= len(masks):
            # Exceeded mask range, check for EOS
            break
        
        mask = masks[i]
        allowed_count = sum(mask)
        
        if allowed_count == 1:
            # Only one allowed base, determine directly
            allowed = get_allowed_indices(mask)
            decoded_symbols.append(allowed[0])
            continue
        
        # Apply mask and renormalize
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
        symbol = masked_index_to_symbol(masked_idx, mask)
        if symbol == -1:
            break  # EOS
        
        decoded_symbols.append(symbol)
    
    return decoded_symbols


# ==============================================================================
# Test Functions
# ==============================================================================

def test_apply_mask_and_renormalize():
    """Test mask and renormalize"""
    print("Testing apply_mask_and_renormalize:")
    
    # Test 1: All allowed
    probs = [0.25, 0.25, 0.25, 0.25]
    mask = [True, True, True, True]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10, f"Sum should be 1, got {sum(result)}"
    print(f"  All allowed: {result} ✓")
    
    # Test 2: Only 2 allowed
    mask = [True, False, True, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10
    assert result[1] == 0 and result[3] == 0
    print(f"  Allow A,G: {result} ✓")
    
    # Test 3: Only 1 allowed
    mask = [False, True, False, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert result[1] == 1.0
    print(f"  Only allow T: {result} ✓")
    
    # Test 4: Uneven probabilities
    probs = [0.5, 0.3, 0.15, 0.05]
    mask = [True, True, False, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10
    assert result[2] == 0 and result[3] == 0
    print(f"  Uneven probabilities: {result} ✓")
    
    # Test 5: All forbidden should raise exception
    try:
        mask = [False, False, False, False]
        apply_mask_and_renormalize(probs, mask)
        assert False, "Should raise exception"
    except ValueError:
        print("  All forbidden raises exception ✓")
    
    print("  apply_mask_and_renormalize test passed ✓")


def test_single_allowed():
    """Test single allowed base case"""
    print("\nTesting single allowed base case:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Test each single base case
    for allowed_base in range(4):
        mask = [i == allowed_base for i in range(4)]
        symbols = [allowed_base] * 10  # 10 identical symbols
        masks = [mask] * 10
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
        print(f"  Only allow base {allowed_base}: bitstream length={len(bits)}, match ✓")
    
    print("  Single base test passed ✓")


def test_two_allowed():
    """Test two allowed bases case"""
    print("\nTesting two allowed bases case:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Test various 2-base combinations
    test_cases = [
        ([True, True, False, False], [0, 1, 0, 1, 0]),  # A, T
        ([True, False, True, False], [0, 2, 0, 2, 2]),  # A, G
        ([False, True, False, True], [1, 3, 1, 3, 1]),  # T, C
        ([False, False, True, True], [2, 3, 2, 2, 3]),  # G, C
    ]
    
    for mask, symbols in test_cases:
        masks = [mask] * len(symbols)
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
        print(f"  Mask {mask}: bitstream length={len(bits)}, match ✓")
    
    print("  2-base test passed ✓")


def test_three_allowed():
    """Test three allowed bases case"""
    print("\nTesting three allowed bases case:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Test various 3-base combinations
    test_cases = [
        ([True, True, True, False], [0, 1, 2, 0, 1]),   # A, T, G
        ([True, True, False, True], [0, 1, 3, 0, 3]),   # A, T, C
        ([True, False, True, True], [0, 2, 3, 0, 2]),   # A, G, C
        ([False, True, True, True], [1, 2, 3, 1, 2]),   # T, G, C
    ]
    
    for mask, symbols in test_cases:
        masks = [mask] * len(symbols)
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
        print(f"  Mask {mask}: bitstream length={len(bits)}, match ✓")
    
    print("  3-base test passed ✓")


def test_four_allowed():
    """Test four allowed bases case"""
    print("\nTesting four allowed bases case:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    mask = [True, True, True, True]
    
    symbols = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]
    masks = [mask] * len(symbols)
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
    print(f"  All allowed: bitstream length={len(bits)}, match ✓")
    
    print("  4-base test passed ✓")


def test_mixed_masks():
    """Test mixed masks case"""
    print("\nTesting mixed masks:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # Masks change with steps
    masks = [
        [True, True, True, True],    # 4 allowed
        [True, True, False, False],  # 2 allowed
        [True, False, False, False], # 1 allowed
        [True, True, True, False],   # 3 allowed
        [False, True, True, True],   # 3 allowed
    ]
    symbols = [0, 1, 0, 2, 1]  # Each symbol must be allowed by corresponding mask
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
    print(f"  Mixed masks: bitstream length={len(bits)}, match ✓")
    
    print("  Mixed masks test passed ✓")


def test_empty_sequence():
    """Test empty sequence"""
    print("\nTesting empty sequence:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    symbols = []
    masks = []
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
    print(f"  Empty sequence: bitstream length={len(bits)}, match ✓")
    
    print("  Empty sequence test passed ✓")


def test_random_exhaustive():
    """Random exhaustive test"""
    print("\nRandom exhaustive test (10000 steps):")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    num_trials = 100
    steps_per_trial = 100
    total_steps = num_trials * steps_per_trial
    
    passed = 0
    failed = 0
    
    for trial in range(num_trials):
        # Randomly generate masks and symbols
        masks = []
        symbols = []
        
        for _ in range(steps_per_trial):
            # Randomly generate mask (at least 1 True)
            while True:
                mask = [random.random() > 0.3 for _ in range(4)]
                if any(mask):
                    break
            
            # Randomly choose from allowed bases
            allowed = get_allowed_indices(mask)
            symbol = random.choice(allowed)
            
            masks.append(mask)
            symbols.append(symbol)
        
        try:
            bits = encode_with_mask(symbols, masks, base_probs)
            decoded = decode_with_mask(bits, masks, base_probs)
            
            if decoded == symbols:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(f"  Failed case {trial}:")
                    print(f"    Input length: {len(symbols)}")
                    print(f"    Decoded length: {len(decoded)}")
                    # Find first mismatch
                    for i in range(min(len(symbols), len(decoded))):
                        if symbols[i] != decoded[i]:
                            print(f"    First mismatch position: {i}")
                            print(f"    Mask: {masks[i]}")
                            print(f"    Expected: {symbols[i]}, Decoded: {decoded[i]}")
                            break
        except Exception as e:
            failed += 1
            if failed <= 5:
                print(f"  Exception case {trial}: {e}")
    
    print(f"  Total trials: {num_trials}")
    print(f"  Total steps: {total_steps}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    print(f"  Success rate: {passed/num_trials*100:.2f}%")
    
    assert failed == 0, f"Random test failed: {failed} cases failed"
    print("  Random exhaustive test passed ✓")


def test_uneven_probs():
    """Test uneven probabilities"""
    print("\nTesting uneven probabilities:")
    
    base_probs = [0.5, 0.3, 0.15, 0.05]
    
    # Test various masks
    test_cases = [
        ([True, True, True, True], [0, 0, 0, 1, 1, 2, 3]),
        ([True, True, False, False], [0, 0, 1, 0, 1]),
        ([False, False, True, True], [2, 3, 2, 2, 3]),
    ]
    
    for mask, symbols in test_cases:
        masks = [mask] * len(symbols)
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"Decode failed: {decoded} != {symbols}"
        print(f"  Mask {mask}: bitstream length={len(bits)}, match ✓")
    
    print("  Uneven probabilities test passed ✓")


def run_all_tests():
    """Run all tests"""
    print("=" * 70)
    print("Masked Arithmetic Encoder Tests")
    print("=" * 70)
    
    test_apply_mask_and_renormalize()
    test_single_allowed()
    test_two_allowed()
    test_three_allowed()
    test_four_allowed()
    test_mixed_masks()
    test_empty_sequence()
    test_uneven_probs()
    test_random_exhaustive()
    
    print("\n" + "=" * 70)
    print("🎉 All tests passed!")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
