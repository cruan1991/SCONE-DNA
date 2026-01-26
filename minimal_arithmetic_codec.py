#!/usr/bin/env python
"""
Standard 32-bit Integer Arithmetic Encoder and Decoder

Implements standard arithmetic coding algorithm:
- Uses 32-bit integer arithmetic
- Inclusive high boundary
- Classic E1/E2/E3 renormalization
- Explicit EOS symbol for termination
- Cumulative frequency model
"""


class StandardArithmeticEncoder:
    """Standard 32-bit integer arithmetic encoder"""
    
    # Constants
    HALF = 0x80000000  # 2^31
    QUARTER = 0x40000000  # 2^30
    THREE_QUARTER = 0xC0000000  # 3 * 2^30
    LOW = 0
    HIGH = 0xFFFFFFFF  # 2^32 - 1
    TOT = 1 << 15  # Total cumulative frequency (32768)
    
    def __init__(self):
        """Initialize encoder"""
        pass
    
    def _build_cumfreq(self, probs):
        """
        Build cumulative frequency table
        
        Args:
            probs: Probability list [p0, p1, ..., pM-1]
        
        Returns:
            cumfreq: Cumulative frequency list [0, f0, f0+f1, ..., TOT]
        """
        if abs(sum(probs) - 1.0) > 1e-5:
            raise ValueError(f"Probabilities must sum to 1.0, got {sum(probs)}")
        
        # Normalize probabilities to handle floating point errors
        total = sum(probs)
        if abs(total - 1.0) > 1e-10:
            probs = [p / total for p in probs]
        
        # Convert probabilities to frequencies (integers)
        freqs = []
        for p in probs:
            freq = max(1, int(p * self.TOT))  # At least 1 to avoid zero frequency
            freqs.append(freq)
        
        # Adjust frequencies to ensure sum equals TOT
        total_freq = sum(freqs)
        if total_freq != self.TOT:
            # Scale proportionally
            scale = self.TOT / total_freq
            freqs = [max(1, int(f * scale)) for f in freqs]
            # Fine-tune to ensure correct sum
            total_freq = sum(freqs)
            diff = self.TOT - total_freq
            if diff != 0:
                # Add difference to largest frequency
                max_idx = freqs.index(max(freqs))
                freqs[max_idx] += diff
        
        # Build cumulative frequency
        cumfreq = [0]
        for f in freqs:
            cumfreq.append(cumfreq[-1] + f)
        
        # Ensure last element equals TOT
        cumfreq[-1] = self.TOT
        
        return cumfreq
    
    def encode(self, symbols, probs):
        """
        Encode symbol sequence
        
        Args:
            symbols: Symbol list (integers, 0 to M-1)
            probs: Probability list [p0, p1, ..., pM-1]
        
        Returns:
            bitstream: Bitstream (list of integers, 0 or 1)
        """
        M = len(probs)
        
        # Build cumulative frequency (including EOS symbol)
        # EOS symbol probability set to 1/10 of minimum probability, but ensure enough frequency
        min_prob = min(probs) if probs else 0.01
        eos_prob = min(min_prob / 10, 0.001)
        # Ensure EOS symbol has at least 64 frequency (about 0.002 of TOT)
        min_eos_freq = 64
        min_eos_prob = min_eos_freq / self.TOT
        eos_prob = max(eos_prob, min_eos_prob)
        
        # Adjust probabilities to include EOS
        scale = 1.0 / (1.0 + eos_prob)
        adjusted_probs = [p * scale for p in probs]
        adjusted_probs.append(eos_prob)
        
        # Build cumulative frequency
        cumfreq = self._build_cumfreq(adjusted_probs)
        
        # Initialize interval
        low = self.LOW
        high = self.HIGH
        pending_bits = 0
        
        # Store output bits
        bitstream = []
        
        # Encode each symbol
        for symbol in symbols:
            if symbol < 0 or symbol >= M:
                raise ValueError(f"Symbol {symbol} out of range [0, {M-1}]")
            
            # Calculate interval size (inclusive high, so +1)
            range_size = high - low + 1
            
            # Determine new interval based on symbol
            symbol_low = cumfreq[symbol]
            symbol_high = cumfreq[symbol + 1]
            
            new_low = low + (range_size * symbol_low) // self.TOT
            new_high = low + (range_size * symbol_high) // self.TOT - 1
            
            # Verify invariant: ensure low <= high
            if new_low > new_high:
                # If this happens, adjust new_high
                new_high = new_low
            
            # Update interval
            low = new_low
            high = new_high
            
            # Verify invariants
            assert 0 <= low <= high <= self.HIGH, \
                f"Interval invariant violated: low={low}, high={high}"
            assert high - low + 1 > 0, \
                f"Invalid interval size: range={high - low + 1}"
            
            # Renormalize (E1/E2/E3)
            while True:
                if high < self.HALF:
                    # E1: Output 0, handle pending bits
                    bitstream.append(0)
                    for _ in range(pending_bits):
                        bitstream.append(1)
                    pending_bits = 0
                    low = low * 2
                    high = high * 2 + 1
                elif low >= self.HALF:
                    # E2: Output 1, handle pending bits
                    bitstream.append(1)
                    for _ in range(pending_bits):
                        bitstream.append(0)
                    pending_bits = 0
                    low = (low - self.HALF) * 2
                    high = (high - self.HALF) * 2 + 1
                elif low >= self.QUARTER and high < self.THREE_QUARTER:
                    # E3: Underflow, increase pending bits
                    pending_bits += 1
                    low = (low - self.QUARTER) * 2
                    high = (high - self.QUARTER) * 2 + 1
                else:
                    break
        
        # Encode EOS symbol
        eos_symbol = M
        range_size = high - low + 1
        symbol_low = cumfreq[eos_symbol]
        symbol_high = cumfreq[eos_symbol + 1]
        
        new_low = low + (range_size * symbol_low) // self.TOT
        new_high = low + (range_size * symbol_high) // self.TOT - 1
        
        # Verify invariant
        if new_low > new_high:
            new_high = new_low
        
        low = new_low
        high = new_high
        
        # Verify invariants
        assert 0 <= low <= high <= self.HIGH, \
            f"Interval invariant violated after EOS: low={low}, high={high}"
        assert high - low + 1 > 0, \
            f"Invalid interval size after EOS: range={high - low + 1}"
        
        # Final renormalization (output all possible bits until cannot continue)
        while True:
            if high < self.HALF:
                bitstream.append(0)
                for _ in range(pending_bits):
                    bitstream.append(1)
                pending_bits = 0
                low = low * 2
                high = high * 2 + 1
            elif low >= self.HALF:
                bitstream.append(1)
                for _ in range(pending_bits):
                    bitstream.append(0)
                pending_bits = 0
                low = (low - self.HALF) * 2
                high = (high - self.HALF) * 2 + 1
            elif low >= self.QUARTER and high < self.THREE_QUARTER:
                pending_bits += 1
                low = (low - self.QUARTER) * 2
                high = (high - self.QUARTER) * 2 + 1
            else:
                break
        
        # finish(): Output pending_bits + 2 bits to disambiguate final interval
        # Standard method: increment pending_bit, then output
        # This ensures decoder can uniquely identify final interval
        pending_bits += 1
        if low < self.QUARTER:
            bitstream.append(0)
            for _ in range(pending_bits):
                bitstream.append(1)
        else:
            bitstream.append(1)
            for _ in range(pending_bits):
                bitstream.append(0)
        
        return bitstream


class StandardArithmeticDecoder:
    """Standard 32-bit integer arithmetic decoder"""
    
    # Constants (same as encoder)
    HALF = 0x80000000
    QUARTER = 0x40000000
    THREE_QUARTER = 0xC0000000
    LOW = 0
    HIGH = 0xFFFFFFFF
    TOT = 1 << 15
    
    def __init__(self):
        """Initialize decoder"""
        pass
    
    def _build_cumfreq(self, probs):
        """Build cumulative frequency table (same as encoder)"""
        if abs(sum(probs) - 1.0) > 1e-5:
            raise ValueError(f"Probabilities must sum to 1.0, got {sum(probs)}")
        
        # Normalize probabilities to handle floating point errors
        total = sum(probs)
        if abs(total - 1.0) > 1e-10:
            probs = [p / total for p in probs]
        
        freqs = []
        for p in probs:
            freq = max(1, int(p * self.TOT))
            freqs.append(freq)
        
        total_freq = sum(freqs)
        if total_freq != self.TOT:
            scale = self.TOT / total_freq
            freqs = [max(1, int(f * scale)) for f in freqs]
            total_freq = sum(freqs)
            diff = self.TOT - total_freq
            if diff != 0:
                max_idx = freqs.index(max(freqs))
                freqs[max_idx] += diff
        
        cumfreq = [0]
        for f in freqs:
            cumfreq.append(cumfreq[-1] + f)
        
        cumfreq[-1] = self.TOT
        
        return cumfreq
    
    def _read_bit(self, bitstream, bit_idx):
        """Read one bit, return 0 if bitstream exhausted"""
        if bit_idx < len(bitstream):
            return bitstream[bit_idx]
        return 0
    
    def decode(self, bitstream, probs):
        """
        Decode bitstream
        
        Args:
            bitstream: Bitstream (list of integers, 0 or 1)
            probs: Probability list [p0, p1, ..., pM-1] (must be same as encoding)
        
        Returns:
            decoded_symbols: Decoded symbol list (excluding EOS)
        """
        M = len(probs)
        
        # Build cumulative frequency (including EOS symbol)
        min_prob = min(probs) if probs else 0.01
        eos_prob = min(min_prob / 10, 0.001)
        # Ensure EOS symbol has at least 64 frequency (about 0.002 of TOT)
        min_eos_freq = 64
        min_eos_prob = min_eos_freq / self.TOT
        eos_prob = max(eos_prob, min_eos_prob)
        scale = 1.0 / (1.0 + eos_prob)
        adjusted_probs = [p * scale for p in probs]
        adjusted_probs.append(eos_prob)
        
        cumfreq = self._build_cumfreq(adjusted_probs)
        
        # Initialize code (value): Read first 32 bits from bitstream (pad with 0 if insufficient)
        code = 0
        bit_idx = 0
        for _ in range(32):
            bit = self._read_bit(bitstream, bit_idx)
            code = code * 2 + bit
            bit_idx += 1
        
        # Initialize interval
        low = self.LOW
        high = self.HIGH
        
        decoded_symbols = []
        
        # Decode symbols until EOS
        while True:
            # Calculate interval size
            range_size = high - low + 1
            
            # Calculate position of current code in cumulative frequency
            # Use same interval calculation method as encoder
            # Find which symbol's interval code falls into
            symbol = len(cumfreq) - 2  # Default to EOS
            
            # Check from back to front, prioritize EOS symbol
            for i in range(len(cumfreq) - 2, -1, -1):
                # Calculate symbol i's interval bounds (exactly same as encoder)
                symbol_low_bound = low + (range_size * cumfreq[i]) // self.TOT
                symbol_high_bound = low + (range_size * cumfreq[i + 1]) // self.TOT - 1
                
                # Ensure valid interval
                if symbol_low_bound > symbol_high_bound:
                    symbol_high_bound = symbol_low_bound
                
                # Check if code is in this interval
                if symbol_low_bound <= code <= symbol_high_bound:
                    symbol = i
                    break
            
            # Ensure symbol in valid range
            if symbol >= len(cumfreq) - 1:
                symbol = len(cumfreq) - 2  # Force to EOS
            
            # Check if EOS symbol
            if symbol == M:
                break
            
            decoded_symbols.append(symbol)
            
            # Update interval
            symbol_low = cumfreq[symbol]
            symbol_high = cumfreq[symbol + 1]
            
            new_low = low + (range_size * symbol_low) // self.TOT
            new_high = low + (range_size * symbol_high) // self.TOT - 1
            
            # Verify invariant: ensure low <= high
            if new_low > new_high:
                new_high = new_low
            
            low = new_low
            high = new_high
            
            # Verify invariants
            assert 0 <= low <= high <= self.HIGH, \
                f"Interval invariant violated: low={low}, high={high}"
            assert high - low + 1 > 0, \
                f"Invalid interval size: range={high - low + 1}"
            
            # Renormalize (symmetric with encoding process)
            while True:
                if high < self.HALF:
                    # E1: MSB is all 0
                    low = low * 2
                    high = high * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = code * 2 + bit
                    bit_idx += 1
                elif low >= self.HALF:
                    # E2: MSB is all 1
                    low = (low - self.HALF) * 2
                    high = (high - self.HALF) * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = (code - self.HALF) * 2 + bit
                    bit_idx += 1
                elif low >= self.QUARTER and high < self.THREE_QUARTER:
                    # E3: Underflow
                    low = (low - self.QUARTER) * 2
                    high = (high - self.QUARTER) * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = (code - self.QUARTER) * 2 + bit
                    bit_idx += 1
                else:
                    break
        
        return decoded_symbols


def test_arithmetic_codec():
    """Comprehensive test of arithmetic encoder and decoder"""
    print("="*70)
    print("Standard 32-bit Integer Arithmetic Encoder/Decoder Tests")
    print("="*70)
    
    encoder = StandardArithmeticEncoder()
    decoder = StandardArithmeticDecoder()
    
    # Test cases
    test_cases = [
        {
            'name': 'Test 1: Simple sequence [0,1,2,3]',
            'symbols': [0, 1, 2, 3],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': 'Test 2: Repeated sequence',
            'symbols': [0, 0, 1, 1, 2, 2, 3, 3],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': 'Test 3: Long sequence (40 symbols)',
            'symbols': [0, 1, 2, 3] * 10,
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': 'Test 4: Single symbol',
            'symbols': [2],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': 'Test 5: Uneven probabilities',
            'symbols': [0, 1, 2, 3],
            'probs': [0.5, 0.25, 0.15, 0.1]
        },
        {
            'name': 'Test 6: Highly skewed probabilities',
            'symbols': [0, 0, 0, 0, 1, 1, 2, 3],
            'probs': [0.7, 0.2, 0.07, 0.03]
        },
        {
            'name': 'Test 7: Random sequence (100 symbols)',
            'symbols': None,  # Will be generated in test
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': 'Test 8: Previously failed last symbol test',
            'symbols': [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
            'probs': [0.25, 0.25, 0.25, 0.25]
        }
    ]
    
    all_passed = True
    total_compression = 0.0
    test_count = 0
    
    import random
    random.seed(42)
    
    for test_case in test_cases:
        print(f"\n{test_case['name']}")
        print("-" * 70)
        
        # Generate random sequence (if needed)
        if test_case['symbols'] is None:
            symbols = [random.randint(0, len(test_case['probs'])-1) 
                      for _ in range(100)]
        else:
            symbols = test_case['symbols']
        
        probs = test_case['probs']
        
        print(f"Input symbol count: {len(symbols)}")
        print(f"Input symbols (first 20): {symbols[:20]}")
        print(f"Probability distribution: {probs}")
        
        # Encode
        bitstream = encoder.encode(symbols, probs)
        print(f"Bitstream length: {len(bitstream)} bits")
        print(f"Bitstream (first 30 bits): {bitstream[:30]}")
        
        # Calculate compression statistics
        original_bits = len(symbols) * 2  # 2 bits per symbol (4 symbols)
        compression_ratio = original_bits / len(bitstream) if len(bitstream) > 0 else 0
        bits_per_symbol = len(bitstream) / len(symbols) if len(symbols) > 0 else 0
        print(f"Original bits: {original_bits} bits")
        print(f"Compression ratio: {compression_ratio:.3f}")
        print(f"Bits per symbol: {bits_per_symbol:.3f}")
        
        total_compression += compression_ratio
        test_count += 1
        
        # Decode
        decoded_symbols = decoder.decode(bitstream, probs)
        print(f"Decoded symbol count: {len(decoded_symbols)}")
        print(f"Decoded symbols (first 20): {decoded_symbols[:20]}")
        
        # Verify
        match = decoded_symbols == symbols
        print(f"Match: {'✓ Passed' if match else '✗ Failed'}")
        
        if not match:
            print(f"Mismatch details:")
            mismatches = [i for i in range(min(len(symbols), len(decoded_symbols))) 
                         if symbols[i] != decoded_symbols[i]]
            if len(mismatches) > 0:
                print(f"  Mismatch count: {len(mismatches)}")
                print(f"  First 10 mismatch positions: {mismatches[:10]}")
                for idx in mismatches[:5]:
                    print(f"    Position {idx}: expected={symbols[idx]}, decoded={decoded_symbols[idx]}")
            if len(decoded_symbols) != len(symbols):
                print(f"  Length mismatch: expected={len(symbols)}, decoded={len(decoded_symbols)}")
            all_passed = False
        
        # Assert correctness
        assert match, f"Test failed: {test_case['name']}"
        assert len(decoded_symbols) == len(symbols), \
            f"Length mismatch: expected={len(symbols)}, decoded={len(decoded_symbols)}"
    
    print("\n" + "="*70)
    print("Test Summary")
    print("="*70)
    print(f"Total tests: {test_count}")
    print(f"Average compression ratio: {total_compression / test_count:.3f}")
    
    if all_passed:
        print("✓ All tests passed!")
    else:
        print("✗ Some tests failed")
    print("="*70)
    
    return all_passed


def verify_invariants(low, high, cumfreq, TOT):
    """Verify invariants"""
    errors = []
    
    # Verify interval
    if not (0 <= low <= high <= 0xFFFFFFFF):
        errors.append(f"Invalid interval: low={low}, high={high}")
    
    range_size = high - low + 1
    if range_size <= 0:
        errors.append(f"Invalid interval size: range={range_size}")
    
    # Verify cumulative frequency
    if len(cumfreq) < 2:
        errors.append(f"Invalid cumfreq length: {len(cumfreq)}")
    
    if cumfreq[0] != 0:
        errors.append(f"First cumfreq element must be 0: {cumfreq[0]}")
    
    if cumfreq[-1] != TOT:
        errors.append(f"Last cumfreq element must be TOT: {cumfreq[-1]}, TOT={TOT}")
    
    # Verify strictly increasing
    for i in range(len(cumfreq) - 1):
        if cumfreq[i] >= cumfreq[i + 1]:
            errors.append(f"Cumfreq not strictly increasing: cumfreq[{i}]={cumfreq[i]} >= cumfreq[{i+1}]={cumfreq[i+1]}")
    
    return errors


def test_m2_exhaustive():
    """Exhaustive test for M=2 symbols"""
    print("="*70)
    print("M=2 Symbol Exhaustive Test")
    print("="*70)
    
    encoder = StandardArithmeticEncoder()
    decoder = StandardArithmeticDecoder()
    
    import random
    random.seed(42)
    
    TOT = StandardArithmeticEncoder.TOT
    num_trials = 10000
    passed = 0
    failed_cases = []
    
    print(f"Running {num_trials} tests...")
    
    for trial in range(num_trials):
        # Random sequence length (0 to 500)
        seq_len = random.randint(0, 500)
        
        # Random probability split: p in {1..TOT-1}, [p, TOT-p]
        p = random.randint(1, TOT - 1)
        prob0 = p / TOT
        prob1 = (TOT - p) / TOT
        probs = [prob0, prob1]
        
        # Generate random sequence
        symbols = [random.randint(0, 1) for _ in range(seq_len)]
        
        try:
            # Encode
            bitstream = encoder.encode(symbols, probs)
            
            # Verify encoder invariants (during encoding)
            # Cannot directly access internal state here, so skip
            
            # Decode
            decoded_symbols = decoder.decode(bitstream, probs)
            
            # Verify
            if decoded_symbols != symbols:
                # Find first mismatch position
                min_len = min(len(symbols), len(decoded_symbols))
                first_mismatch = None
                for i in range(min_len):
                    if symbols[i] != decoded_symbols[i]:
                        first_mismatch = i
                        break
                
                if first_mismatch is None and len(decoded_symbols) != len(symbols):
                    first_mismatch = min_len
                
                failed_cases.append({
                    'trial': trial,
                    'seq_len': seq_len,
                    'probs': probs,
                    'p': p,
                    'symbols': symbols,
                    'decoded': decoded_symbols,
                    'bitstream_len': len(bitstream),
                    'first_mismatch': first_mismatch
                })
                
                # Only save first 10 failed cases
                if len(failed_cases) <= 10:
                    print(f"\n❌ Test {trial} failed:")
                    print(f"  Sequence length: {seq_len}")
                    print(f"  Probabilities: {probs} (p={p})")
                    print(f"  Input symbols (first 20): {symbols[:20]}")
                    print(f"  Decoded symbols (first 20): {decoded_symbols[:20]}")
                    print(f"  Bitstream length: {len(bitstream)}")
                    print(f"  First mismatch position: {first_mismatch}")
            else:
                passed += 1
                
        except Exception as e:
            failed_cases.append({
                'trial': trial,
                'seq_len': seq_len,
                'probs': probs,
                'p': p,
                'symbols': symbols,
                'error': str(e)
            })
            if len(failed_cases) <= 10:
                print(f"\n❌ Test {trial} exception: {e}")
    
    print(f"\n" + "="*70)
    print(f"Test Results:")
    print(f"  Total tests: {num_trials}")
    print(f"  Passed: {passed}")
    print(f"  Failed: {len(failed_cases)}")
    print(f"  Success rate: {passed/num_trials*100:.2f}%")
    
    if failed_cases:
        print(f"\nFailed case details (first {min(5, len(failed_cases))}):")
        for i, case in enumerate(failed_cases[:5]):
            print(f"\nCase {i+1}:")
            print(f"  Trial number: {case['trial']}")
            print(f"  Sequence length: {case['seq_len']}")
            print(f"  Probabilities: {case['probs']} (p={case['p']})")
            if 'error' in case:
                print(f"  Error: {case['error']}")
            else:
                print(f"  Input symbols: {case['symbols']}")
                print(f"  Decoded symbols: {case['decoded']}")
                print(f"  Bitstream length: {case['bitstream_len']}")
                print(f"  First mismatch position: {case['first_mismatch']}")
                
                # Print minimal reproducible case
                print(f"\n  Minimal reproducible case:")
                print(f"    probs = {case['probs']}")
                print(f"    symbols = {case['symbols']}")
                print(f"    bitstream_len = {case['bitstream_len']}")
                if case['first_mismatch'] is not None:
                    print(f"    first_mismatch_index = {case['first_mismatch']}")
    
    print("="*70)
    
    assert len(failed_cases) == 0, f"M=2 test failed: {len(failed_cases)} cases failed"
    print("✓ M=2 exhaustive test all passed!")
    
    return len(failed_cases) == 0


if __name__ == "__main__":
    try:
        # Run standard tests first
        print("Running standard tests...")
        success = test_arithmetic_codec()
        if not success:
            exit(1)
        
        # Run M=2 exhaustive test
        print("\n" + "="*70)
        print("Running M=2 exhaustive test...")
        m2_success = test_m2_exhaustive()
        if not m2_success:
            exit(1)
        
        print("\n" + "="*70)
        print("🎉 All tests passed!")
        print("="*70)
        
    except AssertionError as e:
        print(f"\n❌ Assertion failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    except Exception as e:
        print(f"\n❌ Error occurred: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
