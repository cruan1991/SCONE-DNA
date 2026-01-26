#!/usr/bin/env python
"""
SCONE FSM算术编码器

结合FSM约束控制器和带掩码的算术编码器，实现DNA编码。
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
# 主要编码/解码函数
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
    使用FSM约束编码潜在符号序列
    
    参数:
        latent_symbols: 潜在符号列表（每个符号∈{0,1,2,3}表示要编码的碱基索引）
        base_probs: 基础概率分布[p_A, p_T, p_G, p_C]
        gc_window: GC窗口大小
        gc_low: GC含量下限
        gc_high: GC含量上限
        max_homopolymer: 最大同聚物长度
    
    返回:
        (bitstream, dna_string): 位流和DNA序列字符串
    """
    # 创建FSM约束控制器
    fsm = FSMConstraint(
        gc_window=gc_window,
        gc_low=gc_low,
        gc_high=gc_high,
        max_homopolymer=max_homopolymer
    )
    
    # 创建算术编码器
    encoder = MaskedArithmeticEncoder()
    
    # DNA序列
    dna_bases = []
    
    # 编码每个符号
    for symbol in latent_symbols:
        if symbol < 0 or symbol >= NUM_BASES:
            raise ValueError(f"无效的符号: {symbol}，必须在[0, {NUM_BASES-1}]范围内")
        
        # 获取当前掩码
        mask = fsm.get_mask()
        
        # 检查符号是否被允许
        if not mask[symbol]:
            raise ValueError(f"符号{symbol}在当前FSM状态下被禁止")
        
        # 总是使用算术编码（包括allowed_count==1的情况，为了与解码器保持对称）
        masked_probs = apply_mask_and_renormalize(base_probs, mask)
        freqs = probs_to_freqs_with_mask(masked_probs, mask)
        cumfreq = freqs_to_cumfreq(freqs)
        
        # 将原始符号索引转换为掩码后的索引
        masked_idx = symbol_to_masked_index(symbol, mask)
        
        # 编码
        encoder.encode_symbol(masked_idx, cumfreq, TOT)
        
        # 更新FSM状态
        fsm.update(symbol)
        
        # 添加到DNA序列
        dna_bases.append(INDEX_TO_BASE[symbol])
    
    # 编码EOS
    # 使用当前FSM状态的掩码
    # 注意：总是编码EOS，即使allowed_count==1，因为解码器需要知道何时停止
    mask = fsm.get_mask()
    masked_probs = apply_mask_and_renormalize(base_probs, mask)
    freqs = probs_to_freqs_with_mask(masked_probs, mask)
    cumfreq = freqs_to_cumfreq(freqs)
    eos_idx = len(freqs) - 1
    encoder.encode_symbol(eos_idx, cumfreq, TOT)
    
    # 完成编码
    encoder.finish()
    
    # 返回位流和DNA序列
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
    使用FSM约束解码位流
    
    参数:
        bits: 位流
        base_probs: 基础概率分布
        gc_window: GC窗口大小
        gc_low: GC含量下限
        gc_high: GC含量上限
        max_homopolymer: 最大同聚物长度
        max_symbols: 最大解码符号数（安全限制）
    
    返回:
        (symbols, dna_string): 解码后的符号列表和DNA序列字符串
    """
    if max_symbols is None:
        max_symbols = 100000  # 默认最大限制
    
    # 创建FSM约束控制器
    fsm = FSMConstraint(
        gc_window=gc_window,
        gc_low=gc_low,
        gc_high=gc_high,
        max_homopolymer=max_homopolymer
    )
    
    # 创建算术解码器
    decoder = MaskedArithmeticDecoder()
    decoder.initialize(bits)
    
    # 解码结果
    decoded_symbols = []
    dna_bases = []
    
    # 解码循环
    for _ in range(max_symbols):
        # 获取当前掩码
        mask = fsm.get_mask()
        
        # 总是使用算术解码（即使allowed_count==1，因为需要检查EOS）
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
        base = masked_index_to_symbol(masked_idx, mask)
        if base == -1:
            break  # EOS
        
        # 更新FSM状态
        fsm.update(base)
        
        # 添加到结果
        decoded_symbols.append(base)
        dna_bases.append(INDEX_TO_BASE[base])
    
    # 返回解码结果
    dna_string = ''.join(dna_bases)
    
    return decoded_symbols, dna_string


# ==============================================================================
# 辅助函数
# ==============================================================================

def calculate_gc_content(dna_string: str) -> float:
    """计算DNA序列的GC含量"""
    if not dna_string:
        return 0.0
    gc_count = sum(1 for base in dna_string if base in 'GC')
    return gc_count / len(dna_string)


def calculate_max_homopolymer(dna_string: str) -> int:
    """计算DNA序列中最长的同聚物长度"""
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
    计算DNA序列的窗口GC含量范围
    
    返回:
        (min_gc, max_gc): 最小和最大窗口GC含量
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
# 测试函数
# ==============================================================================

def test_basic_encode_decode():
    """测试基本编码解码"""
    print("测试基本编码解码:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 测试简单序列
    symbols = [0, 1, 2, 3, 0, 1, 2, 3]  # ATGCATGC
    
    bits, dna = encode_fsm(symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=len(symbols)+10)
    
    print(f"  输入符号: {symbols}")
    print(f"  编码DNA: {dna}")
    print(f"  位流长度: {len(bits)}")
    print(f"  解码符号: {decoded}")
    print(f"  解码DNA: {decoded_dna}")
    
    assert decoded == symbols, f"解码失败: {decoded} != {symbols}"
    assert decoded_dna == dna, f"DNA不匹配: {decoded_dna} != {dna}"
    print("  基本编码解码测试通过 ✓")


def test_gc_control():
    """测试GC控制"""
    print("\n测试GC控制:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    gc_window = 20
    gc_low = 0.40
    gc_high = 0.60
    
    # 生成随机序列（从允许的碱基中选择）
    fsm = FSMConstraint(gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10)
    
    symbols = []
    for _ in range(100):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # 编码
    bits, dna = encode_fsm(symbols, base_probs, gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10)
    
    # 检查GC含量
    min_gc, max_gc = calculate_windowed_gc(dna, gc_window)
    overall_gc = calculate_gc_content(dna)
    
    print(f"  序列长度: {len(dna)}")
    print(f"  总体GC含量: {overall_gc:.2%}")
    print(f"  窗口GC范围: [{min_gc:.2%}, {max_gc:.2%}]")
    print(f"  目标范围: [{gc_low:.2%}, {gc_high:.2%}]")
    
    # 对于足够长的序列，窗口GC应该在范围内
    if len(dna) >= gc_window:
        # 允许一些边界容差
        assert min_gc >= gc_low - 0.1, f"最小GC太低: {min_gc:.2%}"
        assert max_gc <= gc_high + 0.1, f"最大GC太高: {max_gc:.2%}"
    
    # 解码
    decoded, decoded_dna = decode_fsm(bits, base_probs, gc_window=gc_window, gc_low=gc_low, gc_high=gc_high, max_homopolymer=10, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"解码失败"
    print("  GC控制测试通过 ✓")


def test_homopolymer_control():
    """测试同聚物控制"""
    print("\n测试同聚物控制:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    max_homopolymer = 3
    
    # 生成随机序列（从允许的碱基中选择）
    fsm = FSMConstraint(gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer)
    
    symbols = []
    for _ in range(200):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # 编码
    bits, dna = encode_fsm(symbols, base_probs, gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer)
    
    # 检查同聚物长度
    max_run = calculate_max_homopolymer(dna)
    
    print(f"  序列长度: {len(dna)}")
    print(f"  最大同聚物长度: {max_run}")
    print(f"  目标限制: {max_homopolymer}")
    
    assert max_run <= max_homopolymer, f"同聚物长度超过限制: {max_run} > {max_homopolymer}"
    
    # 解码
    decoded, decoded_dna = decode_fsm(bits, base_probs, gc_window=100, gc_low=0.0, gc_high=1.0, max_homopolymer=max_homopolymer, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"解码失败"
    print("  同聚物控制测试通过 ✓")


def test_reversibility():
    """测试可逆性"""
    print("\n测试可逆性:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    passed = 0
    failed = 0
    
    # 测试不同参数组合
    test_configs = [
        {'gc_window': 20, 'gc_low': 0.45, 'gc_high': 0.55, 'max_homopolymer': 3},
        {'gc_window': 10, 'gc_low': 0.40, 'gc_high': 0.60, 'max_homopolymer': 2},
        {'gc_window': 30, 'gc_low': 0.35, 'gc_high': 0.65, 'max_homopolymer': 4},
    ]
    
    for config in test_configs:
        for trial in range(10):
            # 生成随机序列
            fsm = FSMConstraint(**config)
            
            seq_len = random.randint(10, 100)
            symbols = []
            
            for _ in range(seq_len):
                mask = fsm.get_mask()
                allowed = [i for i, m in enumerate(mask) if m]
                symbol = random.choice(allowed)
                symbols.append(symbol)
                fsm.update(symbol)
            
            # 编码
            bits, dna = encode_fsm(symbols, base_probs, **config)
            
            # 解码
            decoded, decoded_dna = decode_fsm(bits, base_probs, **config, max_symbols=len(symbols)+10)
            
            if decoded == symbols:
                passed += 1
            else:
                failed += 1
                if failed <= 3:
                    print(f"  失败案例:")
                    print(f"    配置: {config}")
                    print(f"    输入长度: {len(symbols)}")
                    print(f"    解码长度: {len(decoded)}")
    
    print(f"  通过: {passed}, 失败: {failed}")
    assert failed == 0, f"可逆性测试失败: {failed}个案例"
    print("  可逆性测试通过 ✓")


def test_random_1000_steps():
    """测试随机1000步序列"""
    print("\n测试随机1000步序列:")
    
    import random
    random.seed(42)
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    config = {'gc_window': 20, 'gc_low': 0.45, 'gc_high': 0.55, 'max_homopolymer': 3}
    
    # 生成1000步序列
    fsm = FSMConstraint(**config)
    
    symbols = []
    for _ in range(1000):
        mask = fsm.get_mask()
        allowed = [i for i, m in enumerate(mask) if m]
        symbol = random.choice(allowed)
        symbols.append(symbol)
        fsm.update(symbol)
    
    # 编码
    bits, dna = encode_fsm(symbols, base_probs, **config)
    
    print(f"  序列长度: {len(dna)}")
    print(f"  位流长度: {len(bits)}")
    print(f"  每碱基位数: {len(bits)/len(dna):.3f}")
    
    # 检查约束
    gc = calculate_gc_content(dna)
    max_run = calculate_max_homopolymer(dna)
    min_gc, max_gc = calculate_windowed_gc(dna, config['gc_window'])
    
    print(f"  总体GC含量: {gc:.2%}")
    print(f"  窗口GC范围: [{min_gc:.2%}, {max_gc:.2%}]")
    print(f"  最大同聚物: {max_run}")
    
    assert max_run <= config['max_homopolymer'], f"同聚物超限: {max_run}"
    
    # 解码
    decoded, decoded_dna = decode_fsm(bits, base_probs, **config, max_symbols=len(symbols)+10)
    
    assert decoded == symbols, f"解码失败: 长度{len(decoded)} vs {len(symbols)}"
    assert decoded_dna == dna, f"DNA不匹配"
    
    print("  随机1000步测试通过 ✓")


def test_edge_cases():
    """测试边界情况"""
    print("\n测试边界情况:")
    
    base_probs = [0.25, 0.25, 0.25, 0.25]
    
    # 测试空序列
    symbols = []
    bits, dna = encode_fsm(symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=10)
    assert decoded == symbols, f"空序列解码失败"
    print("  空序列 ✓")
    
    # 测试单个符号
    for s in range(4):
        symbols = [s]
        bits, dna = encode_fsm(symbols, base_probs)
        decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=10)
        assert decoded == symbols, f"单符号{s}解码失败"
    print("  单符号 ✓")
    
    # 测试不均匀概率
    base_probs = [0.5, 0.3, 0.15, 0.05]
    symbols = [0, 0, 1, 2, 3, 0, 1, 0]
    
    # 生成有效序列
    fsm = FSMConstraint()
    valid_symbols = []
    for s in symbols:
        mask = fsm.get_mask()
        if mask[s]:
            valid_symbols.append(s)
            fsm.update(s)
    
    bits, dna = encode_fsm(valid_symbols, base_probs)
    decoded, decoded_dna = decode_fsm(bits, base_probs, max_symbols=len(valid_symbols)+10)
    assert decoded == valid_symbols, f"不均匀概率解码失败"
    print("  不均匀概率 ✓")
    
    print("  边界情况测试通过 ✓")


def run_all_tests():
    """运行所有测试"""
    print("=" * 70)
    print("SCONE FSM算术编码器测试")
    print("=" * 70)
    
    test_basic_encode_decode()
    test_gc_control()
    test_homopolymer_control()
    test_reversibility()
    test_random_1000_steps()
    test_edge_cases()
    
    print("\n" + "=" * 70)
    print("🎉 所有测试通过！")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
