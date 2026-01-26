#!/usr/bin/env python
"""
带静态掩码层的DNA碱基算术编码器

在标准算术编码器基础上添加掩码支持，用于FSM约束下的DNA编码。
支持动态掩码，允许集大小可以是4、3、2或1。
"""

from typing import List, Tuple, Optional


# ==============================================================================
# 常量定义
# ==============================================================================

# DNA碱基映射
BASE_TO_INDEX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INDEX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
NUM_BASES = 4

# 算术编码常量
HALF = 0x80000000  # 2^31
QUARTER = 0x40000000  # 2^30
THREE_QUARTER = 0xC0000000  # 3 * 2^30
LOW_INIT = 0
HIGH_INIT = 0xFFFFFFFF  # 2^32 - 1
TOT = 1 << 15  # 累积频率总数 (32768)


# ==============================================================================
# 掩码和概率处理函数
# ==============================================================================

def apply_mask_and_renormalize(probs: List[float], mask: List[bool]) -> List[float]:
    """
    应用掩码并重新归一化概率
    
    参数:
        probs: 4个浮点数的列表，表示A,T,G,C的概率
        mask: 4个布尔值的列表，True表示允许，False表示禁止
    
    返回:
        重新归一化后的概率列表（长度4，被掩码的位置为0）
    
    异常:
        ValueError: 如果所有mask[i]都为False
    """
    if len(probs) != 4 or len(mask) != 4:
        raise ValueError("probs和mask长度必须为4")
    
    # 检查是否至少有一个允许的碱基
    if not any(mask):
        raise ValueError("至少需要一个允许的碱基（mask不能全为False）")
    
    # 应用掩码
    masked_probs = [p if m else 0.0 for p, m in zip(probs, mask)]
    
    # 计算总和
    total = sum(masked_probs)
    
    # 重新归一化
    if total > 0:
        normalized_probs = [p / total for p in masked_probs]
    else:
        # 如果总和为0（所有允许的碱基概率都是0），均匀分配
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
    将概率转换为频率（用于算术编码）
    
    参数:
        probs: 重新归一化后的概率列表（来自apply_mask_and_renormalize）
        mask: 掩码列表
        tot: 累积频率总数
        eos_freq: EOS符号的频率
    
    返回:
        频率列表，长度为允许碱基数+1（最后一个是EOS）
    """
    # 计算允许的碱基数量
    allowed_count = sum(mask)
    
    if allowed_count == 0:
        raise ValueError("至少需要一个允许的碱基")
    
    # 计算可用于碱基的总频率
    available_freq = tot - eos_freq
    
    # 只为允许的碱基分配频率
    freqs = []
    for i, (p, m) in enumerate(zip(probs, mask)):
        if m:
            freq = max(1, int(p * available_freq))
            freqs.append(freq)
    
    # 调整频率以确保总和正确
    total_freq = sum(freqs)
    target_freq = available_freq
    
    if total_freq != target_freq and len(freqs) > 0:
        diff = target_freq - total_freq
        # 将差值加到最大频率上
        max_idx = freqs.index(max(freqs))
        freqs[max_idx] += diff
    
    # 添加EOS频率
    freqs.append(eos_freq)
    
    return freqs


def freqs_to_cumfreq(freqs: List[int]) -> List[int]:
    """将频率转换为累积频率"""
    cumfreq = [0]
    for f in freqs:
        cumfreq.append(cumfreq[-1] + f)
    return cumfreq


def get_allowed_indices(mask: List[bool]) -> List[int]:
    """获取允许的碱基索引列表"""
    return [i for i, m in enumerate(mask) if m]


def symbol_to_masked_index(symbol: int, mask: List[bool]) -> int:
    """将原始符号索引转换为掩码后的索引"""
    if not mask[symbol]:
        raise ValueError(f"符号{symbol}被掩码禁止")
    
    allowed = get_allowed_indices(mask)
    return allowed.index(symbol)


def masked_index_to_symbol(masked_idx: int, mask: List[bool]) -> int:
    """将掩码后的索引转换为原始符号索引"""
    allowed = get_allowed_indices(mask)
    if masked_idx >= len(allowed):
        return -1  # EOS
    return allowed[masked_idx]


# ==============================================================================
# 带掩码的算术编码器
# ==============================================================================

class MaskedArithmeticEncoder:
    """带掩码的算术编码器"""
    
    def __init__(self):
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.pending_bits = 0
        self.bitstream = []
    
    def reset(self):
        """重置编码器状态"""
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.pending_bits = 0
        self.bitstream = []
    
    def _output_bit(self, bit: int):
        """输出一位并处理待处理位"""
        self.bitstream.append(bit)
        for _ in range(self.pending_bits):
            self.bitstream.append(1 - bit)
        self.pending_bits = 0
    
    def _renormalize(self):
        """重新归一化"""
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
        """编码一个符号"""
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
        """完成编码，输出最终位"""
        self.pending_bits += 1
        if self.low < QUARTER:
            self._output_bit(0)
        else:
            self._output_bit(1)
    
    def get_bitstream(self) -> List[int]:
        """获取位流"""
        return self.bitstream


class MaskedArithmeticDecoder:
    """带掩码的算术解码器"""
    
    def __init__(self):
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.code = 0
        self.bitstream = []
        self.bit_idx = 0
    
    def reset(self):
        """重置解码器状态"""
        self.low = LOW_INIT
        self.high = HIGH_INIT
        self.code = 0
        self.bit_idx = 0
    
    def _read_bit(self) -> int:
        """读取一位"""
        if self.bit_idx < len(self.bitstream):
            bit = self.bitstream[self.bit_idx]
            self.bit_idx += 1
            return bit
        return 0
    
    def initialize(self, bitstream: List[int]):
        """初始化解码器"""
        self.bitstream = bitstream
        self.bit_idx = 0
        self.low = LOW_INIT
        self.high = HIGH_INIT
        
        # 读取前32位
        self.code = 0
        for _ in range(32):
            self.code = self.code * 2 + self._read_bit()
    
    def _renormalize(self):
        """重新归一化"""
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
        """解码一个符号"""
        range_size = self.high - self.low + 1
        
        # 找到code落在哪个符号的区间内
        num_symbols = len(cumfreq) - 1
        symbol = num_symbols - 1  # 默认为最后一个（EOS）
        
        for i in range(num_symbols - 1, -1, -1):
            symbol_low_bound = self.low + (range_size * cumfreq[i]) // tot
            symbol_high_bound = self.low + (range_size * cumfreq[i + 1]) // tot - 1
            
            if symbol_low_bound > symbol_high_bound:
                symbol_high_bound = symbol_low_bound
            
            if symbol_low_bound <= self.code <= symbol_high_bound:
                symbol = i
                break
        
        # 更新区间
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
# 主要编码/解码函数
# ==============================================================================

def encode_with_mask(
    symbols: List[int],
    masks: List[List[bool]],
    base_probs: List[float]
) -> List[int]:
    """
    带掩码的编码
    
    参数:
        symbols: 符号列表，每个符号∈{0,1,2,3}表示A,T,G,C
        masks: 掩码列表，与symbols等长
        base_probs: 基础概率分布[p_A, p_T, p_G, p_C]
    
    返回:
        位流（整数列表）
    """
    if len(symbols) != len(masks):
        raise ValueError("symbols和masks长度必须相同")
    
    if len(symbols) == 0:
        # 空序列，只编码EOS
        encoder = MaskedArithmeticEncoder()
        # 使用默认掩码（全部允许）
        default_mask = [True] * 4
        masked_probs = apply_mask_and_renormalize(base_probs, default_mask)
        freqs = probs_to_freqs_with_mask(masked_probs, default_mask)
        cumfreq = freqs_to_cumfreq(freqs)
        # EOS索引是最后一个
        eos_idx = len(freqs) - 1
        encoder.encode_symbol(eos_idx, cumfreq, TOT)
        encoder.finish()
        return encoder.get_bitstream()
    
    encoder = MaskedArithmeticEncoder()
    
    # 编码每个符号
    for i, (symbol, mask) in enumerate(zip(symbols, masks)):
        if not mask[symbol]:
            raise ValueError(f"位置{i}的符号{symbol}被掩码禁止")
        
        allowed_count = sum(mask)
        
        if allowed_count == 1:
            # 只有一个允许的碱基，不需要编码
            # 解码器知道掩码，可以直接确定
            continue
        
        # 应用掩码并重新归一化
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # 将原始符号索引转换为掩码后的索引
        masked_idx = symbol_to_masked_index(symbol, mask)
        
        # 编码
        encoder.encode_symbol(masked_idx, cumfreq, TOT)
    
    # 编码EOS
    # 使用最后一个掩码（或默认掩码）
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
    带掩码的解码
    
    参数:
        bits: 位流
        masks: 掩码列表
        base_probs: 基础概率分布
        max_symbols: 最大解码符号数（可选，用于安全）
    
    返回:
        解码后的符号列表
    """
    if max_symbols is None:
        max_symbols = len(masks) + 1  # 允许稍微多一点
    
    if len(masks) == 0:
        # 空掩码，解码EOS
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
            raise ValueError("解码错误：期望EOS但得到其他符号")
    
    decoder = MaskedArithmeticDecoder()
    decoder.initialize(bits)
    
    decoded_symbols = []
    
    for i in range(max_symbols):
        if i >= len(masks):
            # 超出掩码范围，检查EOS
            break
        
        mask = masks[i]
        allowed_count = sum(mask)
        
        if allowed_count == 1:
            # 只有一个允许的碱基，直接确定
            allowed = get_allowed_indices(mask)
            decoded_symbols.append(allowed[0])
            continue
        
        # 应用掩码并重新归一化
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # 解码
        masked_idx = decoder.decode_symbol(cumfreq, TOT)
        
        # EOS检查
        eos_idx = len(freqs) - 1
        if masked_idx == eos_idx:
            break
        
        # 将掩码后的索引转换为原始符号索引
        symbol = masked_index_to_symbol(masked_idx, mask)
        if symbol == -1:
            break  # EOS
        
        decoded_symbols.append(symbol)
    
    return decoded_symbols


# ==============================================================================
# 测试函数
# ==============================================================================

def test_apply_mask_and_renormalize():
    """测试掩码和重新归一化"""
    print("测试 apply_mask_and_renormalize:")
    
    # 测试1：全部允许
    probs = [0.25, 0.25, 0.25, 0.25]
    mask = [True, True, True, True]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10, f"总和应为1，实际为{sum(result)}"
    print(f"  全部允许: {result} ✓")
    
    # 测试2：只允许2个
    mask = [True, False, True, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10
    assert result[1] == 0 and result[3] == 0
    print(f"  允许A,G: {result} ✓")
    
    # 测试3：只允许1个
    mask = [False, True, False, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert result[1] == 1.0
    print(f"  只允许T: {result} ✓")
    
    # 测试4：不均匀概率
    probs = [0.5, 0.3, 0.15, 0.05]
    mask = [True, True, False, False]
    result = apply_mask_and_renormalize(probs, mask)
    assert abs(sum(result) - 1.0) < 1e-10
    assert result[2] == 0 and result[3] == 0
    print(f"  不均匀概率: {result} ✓")
    
    # 测试5：全部禁止应抛出异常
    try:
        mask = [False, False, False, False]
        apply_mask_and_renormalize(probs, mask)
        assert False, "应该抛出异常"
    except ValueError:
        print("  全部禁止抛出异常 ✓")
    
    print("  apply_mask_and_renormalize 测试通过 ✓")


def test_single_allowed():
    """测试只允许1个碱基的情况"""
    print("\n测试只允许1个碱基的情况:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 测试各种单碱基情况
    for allowed_base in range(4):
        mask = [i == allowed_base for i in range(4)]
        symbols = [allowed_base] * 10  # 10个相同的符号
        masks = [mask] * 10
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
        print(f"  只允许碱基{allowed_base}: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  单碱基测试通过 ✓")


def test_two_allowed():
    """测试允许2个碱基的情况"""
    print("\n测试允许2个碱基的情况:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 测试各种2碱基组合
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
        
        assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
        print(f"  掩码{mask}: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  2碱基测试通过 ✓")


def test_three_allowed():
    """测试允许3个碱基的情况"""
    print("\n测试允许3个碱基的情况:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 测试各种3碱基组合
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
        
        assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
        print(f"  掩码{mask}: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  3碱基测试通过 ✓")


def test_four_allowed():
    """测试允许4个碱基的情况"""
    print("\n测试允许4个碱基的情况:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    mask = [True, True, True, True]
    
    symbols = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]
    masks = [mask] * len(symbols)
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
    print(f"  全部允许: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  4碱基测试通过 ✓")


def test_mixed_masks():
    """测试混合掩码的情况"""
    print("\n测试混合掩码:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 掩码随步骤变化
    masks = [
        [True, True, True, True],    # 4个允许
        [True, True, False, False],  # 2个允许
        [True, False, False, False], # 1个允许
        [True, True, True, False],   # 3个允许
        [False, True, True, True],   # 3个允许
    ]
    symbols = [0, 1, 0, 2, 1]  # 每个符号必须被对应掩码允许
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
    print(f"  混合掩码: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  混合掩码测试通过 ✓")


def test_empty_sequence():
    """测试空序列"""
    print("\n测试空序列:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    symbols = []
    masks = []
    
    bits = encode_with_mask(symbols, masks, base_probs)
    decoded = decode_with_mask(bits, masks, base_probs)
    
    assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
    print(f"  空序列: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  空序列测试通过 ✓")


def test_random_exhaustive():
    """随机详尽测试"""
    print("\n随机详尽测试 (10000步):")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    num_trials = 100
    steps_per_trial = 100
    total_steps = num_trials * steps_per_trial
    
    passed = 0
    failed = 0
    
    for trial in range(num_trials):
        # 随机生成掩码和符号
        masks = []
        symbols = []
        
        for _ in range(steps_per_trial):
            # 随机生成掩码（至少1个True）
            while True:
                mask = [random.random() > 0.3 for _ in range(4)]
                if any(mask):
                    break
            
            # 从允许的碱基中随机选择
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
                    print(f"  失败案例 {trial}:")
                    print(f"    输入长度: {len(symbols)}")
                    print(f"    解码长度: {len(decoded)}")
                    # 找到第一个不匹配
                    for i in range(min(len(symbols), len(decoded))):
                        if symbols[i] != decoded[i]:
                            print(f"    第一个不匹配位置: {i}")
                            print(f"    掩码: {masks[i]}")
                            print(f"    期望: {symbols[i]}, 解码: {decoded[i]}")
                            break
        except Exception as e:
            failed += 1
            if failed <= 5:
                print(f"  异常案例 {trial}: {e}")
    
    print(f"  总试验数: {num_trials}")
    print(f"  总步骤数: {total_steps}")
    print(f"  通过: {passed}")
    print(f"  失败: {failed}")
    print(f"  成功率: {passed/num_trials*100:.2f}%")
    
    assert failed == 0, f"随机测试失败: {failed}个案例失败"
    print("  随机详尽测试通过 ✓")


def test_uneven_probs():
    """测试不均匀概率"""
    print("\n测试不均匀概率:")
    
    base_probs = [0.5, 0.3, 0.15, 0.05]
    
    # 测试各种掩码
    test_cases = [
        ([True, True, True, True], [0, 0, 0, 1, 1, 2, 3]),
        ([True, True, False, False], [0, 0, 1, 0, 1]),
        ([False, False, True, True], [2, 3, 2, 2, 3]),
    ]
    
    for mask, symbols in test_cases:
        masks = [mask] * len(symbols)
        
        bits = encode_with_mask(symbols, masks, base_probs)
        decoded = decode_with_mask(bits, masks, base_probs)
        
        assert decoded == symbols, f"解码失败：{decoded} != {symbols}"
        print(f"  掩码{mask}: 位流长度={len(bits)}, 匹配 ✓")
    
    print("  不均匀概率测试通过 ✓")


def run_all_tests():
    """运行所有测试"""
    print("=" * 70)
    print("带掩码的算术编码器测试")
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
    print("🎉 所有测试通过！")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
