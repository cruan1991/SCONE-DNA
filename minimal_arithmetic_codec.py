#!/usr/bin/env python
"""
标准32位整数算术编码器和解码器

实现标准算术编码算法：
- 使用32位整数算术
- 包含high（inclusive）
- 经典E1/E2/E3重新归一化
- 显式EOS符号用于终止
- 累积频率模型
"""


class StandardArithmeticEncoder:
    """标准32位整数算术编码器"""
    
    # 常量定义
    HALF = 0x80000000  # 2^31
    QUARTER = 0x40000000  # 2^30
    THREE_QUARTER = 0xC0000000  # 3 * 2^30
    LOW = 0
    HIGH = 0xFFFFFFFF  # 2^32 - 1
    TOT = 1 << 15  # 累积频率总数（32768）
    
    def __init__(self):
        """初始化编码器"""
        pass
    
    def _build_cumfreq(self, probs):
        """
        构建累积频率表
        
        参数:
            probs: 概率列表 [p0, p1, ..., pM-1]
        
        返回:
            cumfreq: 累积频率列表 [0, f0, f0+f1, ..., TOT]
        """
        if abs(sum(probs) - 1.0) > 1e-5:
            raise ValueError(f"概率总和必须为1.0，当前为{sum(probs)}")
        
        # 归一化概率以处理浮点数误差
        total = sum(probs)
        if abs(total - 1.0) > 1e-10:
            probs = [p / total for p in probs]
        
        # 将概率转换为频率（整数）
        freqs = []
        for p in probs:
            freq = max(1, int(p * self.TOT))  # 至少为1，避免零频率
            freqs.append(freq)
        
        # 调整频率以确保总和等于TOT
        total_freq = sum(freqs)
        if total_freq != self.TOT:
            # 按比例调整
            scale = self.TOT / total_freq
            freqs = [max(1, int(f * scale)) for f in freqs]
            # 微调以确保总和正确
            total_freq = sum(freqs)
            diff = self.TOT - total_freq
            if diff != 0:
                # 将差值加到最大频率上
                max_idx = freqs.index(max(freqs))
                freqs[max_idx] += diff
        
        # 构建累积频率
        cumfreq = [0]
        for f in freqs:
            cumfreq.append(cumfreq[-1] + f)
        
        # 确保最后一个等于TOT
        cumfreq[-1] = self.TOT
        
        return cumfreq
    
    def encode(self, symbols, probs):
        """
        编码符号序列
        
        参数:
            symbols: 符号列表（整数，0到M-1）
            probs: 概率列表 [p0, p1, ..., pM-1]
        
        返回:
            bitstream: 位流（整数列表，0或1）
        """
        M = len(probs)
        
        # 构建累积频率（包含EOS符号）
        # EOS符号的概率设为最小概率的1/10，但至少保证有足够的频率
        min_prob = min(probs) if probs else 0.01
        eos_prob = min(min_prob / 10, 0.001)
        # 确保EOS符号至少有64的频率（约0.002的TOT）
        min_eos_freq = 64
        min_eos_prob = min_eos_freq / self.TOT
        eos_prob = max(eos_prob, min_eos_prob)
        
        # 调整概率以包含EOS
        scale = 1.0 / (1.0 + eos_prob)
        adjusted_probs = [p * scale for p in probs]
        adjusted_probs.append(eos_prob)
        
        # 构建累积频率
        cumfreq = self._build_cumfreq(adjusted_probs)
        
        # 初始化区间
        low = self.LOW
        high = self.HIGH
        pending_bits = 0
        
        # 存储输出的位
        bitstream = []
        
        # 编码每个符号
        for symbol in symbols:
            if symbol < 0 or symbol >= M:
                raise ValueError(f"符号 {symbol} 超出范围 [0, {M-1}]")
            
            # 计算区间大小（包含high，所以+1）
            range_size = high - low + 1
            
            # 根据符号确定新区间
            symbol_low = cumfreq[symbol]
            symbol_high = cumfreq[symbol + 1]
            
            new_low = low + (range_size * symbol_low) // self.TOT
            new_high = low + (range_size * symbol_high) // self.TOT - 1
            
            # 验证不变量：确保 low <= high
            if new_low > new_high:
                # 如果出现这种情况，调整new_high
                new_high = new_low
            
            # 更新区间
            low = new_low
            high = new_high
            
            # 验证不变量
            assert 0 <= low <= high <= self.HIGH, \
                f"区间不变量违反: low={low}, high={high}"
            assert high - low + 1 > 0, \
                f"区间大小无效: range={high - low + 1}"
            
            # 重新归一化（E1/E2/E3）
            while True:
                if high < self.HALF:
                    # E1: 输出0，处理待处理位
                    bitstream.append(0)
                    for _ in range(pending_bits):
                        bitstream.append(1)
                    pending_bits = 0
                    low = low * 2
                    high = high * 2 + 1
                elif low >= self.HALF:
                    # E2: 输出1，处理待处理位
                    bitstream.append(1)
                    for _ in range(pending_bits):
                        bitstream.append(0)
                    pending_bits = 0
                    low = (low - self.HALF) * 2
                    high = (high - self.HALF) * 2 + 1
                elif low >= self.QUARTER and high < self.THREE_QUARTER:
                    # E3: 下溢，增加待处理位
                    pending_bits += 1
                    low = (low - self.QUARTER) * 2
                    high = (high - self.QUARTER) * 2 + 1
                else:
                    break
        
        # 编码EOS符号
        eos_symbol = M
        range_size = high - low + 1
        symbol_low = cumfreq[eos_symbol]
        symbol_high = cumfreq[eos_symbol + 1]
        
        new_low = low + (range_size * symbol_low) // self.TOT
        new_high = low + (range_size * symbol_high) // self.TOT - 1
        
        # 验证不变量
        if new_low > new_high:
            new_high = new_low
        
        low = new_low
        high = new_high
        
        # 验证不变量
        assert 0 <= low <= high <= self.HIGH, \
            f"EOS编码后区间不变量违反: low={low}, high={high}"
        assert high - low + 1 > 0, \
            f"EOS编码后区间大小无效: range={high - low + 1}"
        
        # 最终重新归一化（输出所有可能的位，直到无法继续）
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
        
        # finish(): 输出pending_bits + 2位以消除最终区间的歧义
        # 标准方法：增加一个pending_bit，然后输出
        # 这确保解码器能够唯一识别最终区间
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
    """标准32位整数算术解码器"""
    
    # 常量定义（与编码器相同）
    HALF = 0x80000000
    QUARTER = 0x40000000
    THREE_QUARTER = 0xC0000000
    LOW = 0
    HIGH = 0xFFFFFFFF
    TOT = 1 << 15
    
    def __init__(self):
        """初始化解码器"""
        pass
    
    def _build_cumfreq(self, probs):
        """构建累积频率表（与编码器相同）"""
        if abs(sum(probs) - 1.0) > 1e-5:
            raise ValueError(f"概率总和必须为1.0，当前为{sum(probs)}")
        
        # 归一化概率以处理浮点数误差
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
        """读取一位，如果位流耗尽则返回0"""
        if bit_idx < len(bitstream):
            return bitstream[bit_idx]
        return 0
    
    def decode(self, bitstream, probs):
        """
        解码位流
        
        参数:
            bitstream: 位流（整数列表，0或1）
            probs: 概率列表 [p0, p1, ..., pM-1]（必须与编码时相同）
        
        返回:
            decoded_symbols: 解码后的符号列表（不包含EOS）
        """
        M = len(probs)
        
        # 构建累积频率（包含EOS符号）
        min_prob = min(probs) if probs else 0.01
        eos_prob = min(min_prob / 10, 0.001)
        # 确保EOS符号至少有64的频率（约0.002的TOT）
        min_eos_freq = 64
        min_eos_prob = min_eos_freq / self.TOT
        eos_prob = max(eos_prob, min_eos_prob)
        scale = 1.0 / (1.0 + eos_prob)
        adjusted_probs = [p * scale for p in probs]
        adjusted_probs.append(eos_prob)
        
        cumfreq = self._build_cumfreq(adjusted_probs)
        
        # 初始化code（value）：从位流开头读取前32位（如果位流不足32位，用0填充）
        code = 0
        bit_idx = 0
        for _ in range(32):
            bit = self._read_bit(bitstream, bit_idx)
            code = code * 2 + bit
            bit_idx += 1
        
        # 初始化区间
        low = self.LOW
        high = self.HIGH
        
        decoded_symbols = []
        
        # 解码符号直到遇到EOS
        while True:
            # 计算区间大小
            range_size = high - low + 1
            
            # 计算当前code在累积频率中的位置
            # 使用与编码器完全相同的区间计算方法
            # 找到code落在哪个符号的区间内
            symbol = len(cumfreq) - 2  # 默认为EOS
            
            # 从后往前检查，优先检查EOS符号
            for i in range(len(cumfreq) - 2, -1, -1):
                # 计算符号i的区间边界（与编码器完全一致）
                symbol_low_bound = low + (range_size * cumfreq[i]) // self.TOT
                symbol_high_bound = low + (range_size * cumfreq[i + 1]) // self.TOT - 1
                
                # 确保区间有效
                if symbol_low_bound > symbol_high_bound:
                    symbol_high_bound = symbol_low_bound
                
                # 检查code是否在这个区间内
                if symbol_low_bound <= code <= symbol_high_bound:
                    symbol = i
                    break
            
            # 确保symbol在有效范围内
            if symbol >= len(cumfreq) - 1:
                symbol = len(cumfreq) - 2  # 强制为EOS
            
            # 检查是否是EOS符号
            if symbol == M:
                break
            
            decoded_symbols.append(symbol)
            
            # 更新区间
            symbol_low = cumfreq[symbol]
            symbol_high = cumfreq[symbol + 1]
            
            new_low = low + (range_size * symbol_low) // self.TOT
            new_high = low + (range_size * symbol_high) // self.TOT - 1
            
            # 验证不变量：确保 low <= high
            if new_low > new_high:
                new_high = new_low
            
            low = new_low
            high = new_high
            
            # 验证不变量
            assert 0 <= low <= high <= self.HIGH, \
                f"区间不变量违反: low={low}, high={high}"
            assert high - low + 1 > 0, \
                f"区间大小无效: range={high - low + 1}"
            
            # 重新归一化（与编码过程对称）
            while True:
                if high < self.HALF:
                    # E1: MSB都是0
                    low = low * 2
                    high = high * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = code * 2 + bit
                    bit_idx += 1
                elif low >= self.HALF:
                    # E2: MSB都是1
                    low = (low - self.HALF) * 2
                    high = (high - self.HALF) * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = (code - self.HALF) * 2 + bit
                    bit_idx += 1
                elif low >= self.QUARTER and high < self.THREE_QUARTER:
                    # E3: 下溢
                    low = (low - self.QUARTER) * 2
                    high = (high - self.QUARTER) * 2 + 1
                    bit = self._read_bit(bitstream, bit_idx)
                    code = (code - self.QUARTER) * 2 + bit
                    bit_idx += 1
                else:
                    break
        
        return decoded_symbols


def test_arithmetic_codec():
    """全面测试算术编码器和解码器"""
    print("="*70)
    print("标准32位整数算术编码器/解码器测试")
    print("="*70)
    
    encoder = StandardArithmeticEncoder()
    decoder = StandardArithmeticDecoder()
    
    # 测试用例
    test_cases = [
        {
            'name': '测试1: 简单序列 [0,1,2,3]',
            'symbols': [0, 1, 2, 3],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': '测试2: 重复序列',
            'symbols': [0, 0, 1, 1, 2, 2, 3, 3],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': '测试3: 长序列（40个符号）',
            'symbols': [0, 1, 2, 3] * 10,
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': '测试4: 单符号',
            'symbols': [2],
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': '测试5: 不均匀概率',
            'symbols': [0, 1, 2, 3],
            'probs': [0.5, 0.25, 0.15, 0.1]
        },
        {
            'name': '测试6: 极端不均匀概率',
            'symbols': [0, 0, 0, 0, 1, 1, 2, 3],
            'probs': [0.7, 0.2, 0.07, 0.03]
        },
        {
            'name': '测试7: 随机序列（100个符号）',
            'symbols': None,  # 将在测试中生成
            'probs': [0.25, 0.25, 0.25, 0.25]
        },
        {
            'name': '测试8: 之前失败的最后一个符号测试',
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
        
        # 生成随机序列（如果需要）
        if test_case['symbols'] is None:
            symbols = [random.randint(0, len(test_case['probs'])-1) 
                      for _ in range(100)]
        else:
            symbols = test_case['symbols']
        
        probs = test_case['probs']
        
        print(f"输入符号数量: {len(symbols)}")
        print(f"输入符号 (前20个): {symbols[:20]}")
        print(f"概率分布: {probs}")
        
        # 编码
        bitstream = encoder.encode(symbols, probs)
        print(f"位流长度: {len(bitstream)} 位")
        print(f"位流 (前30位): {bitstream[:30]}")
        
        # 计算压缩统计
        original_bits = len(symbols) * 2  # 每个符号2位（4个符号）
        compression_ratio = original_bits / len(bitstream) if len(bitstream) > 0 else 0
        bits_per_symbol = len(bitstream) / len(symbols) if len(symbols) > 0 else 0
        print(f"原始位数: {original_bits} 位")
        print(f"压缩比: {compression_ratio:.3f}")
        print(f"每符号位数: {bits_per_symbol:.3f}")
        
        total_compression += compression_ratio
        test_count += 1
        
        # 解码
        decoded_symbols = decoder.decode(bitstream, probs)
        print(f"解码符号数量: {len(decoded_symbols)}")
        print(f"解码符号 (前20个): {decoded_symbols[:20]}")
        
        # 验证
        match = decoded_symbols == symbols
        print(f"匹配: {'✓ 通过' if match else '✗ 失败'}")
        
        if not match:
            print(f"不匹配详情:")
            mismatches = [i for i in range(min(len(symbols), len(decoded_symbols))) 
                         if symbols[i] != decoded_symbols[i]]
            if len(mismatches) > 0:
                print(f"  不匹配位置数量: {len(mismatches)}")
                print(f"  前10个不匹配位置: {mismatches[:10]}")
                for idx in mismatches[:5]:
                    print(f"    位置{idx}: 期望={symbols[idx]}, 解码={decoded_symbols[idx]}")
            if len(decoded_symbols) != len(symbols):
                print(f"  长度不匹配: 期望={len(symbols)}, 解码={len(decoded_symbols)}")
            all_passed = False
        
        # 断言正确性
        assert match, f"测试失败: {test_case['name']}"
        assert len(decoded_symbols) == len(symbols), \
            f"长度不匹配: 期望={len(symbols)}, 解码={len(decoded_symbols)}"
    
    print("\n" + "="*70)
    print("测试总结")
    print("="*70)
    print(f"总测试数: {test_count}")
    print(f"平均压缩比: {total_compression / test_count:.3f}")
    
    if all_passed:
        print("✓ 所有测试通过！")
    else:
        print("✗ 部分测试失败")
    print("="*70)
    
    return all_passed


def verify_invariants(low, high, cumfreq, TOT):
    """验证不变量"""
    errors = []
    
    # 验证区间
    if not (0 <= low <= high <= 0xFFFFFFFF):
        errors.append(f"区间无效: low={low}, high={high}")
    
    range_size = high - low + 1
    if range_size <= 0:
        errors.append(f"区间大小无效: range={range_size}")
    
    # 验证累积频率
    if len(cumfreq) < 2:
        errors.append(f"累积频率长度无效: {len(cumfreq)}")
    
    if cumfreq[0] != 0:
        errors.append(f"累积频率第一个元素必须为0: {cumfreq[0]}")
    
    if cumfreq[-1] != TOT:
        errors.append(f"累积频率最后一个元素必须为TOT: {cumfreq[-1]}, TOT={TOT}")
    
    # 验证严格递增
    for i in range(len(cumfreq) - 1):
        if cumfreq[i] >= cumfreq[i + 1]:
            errors.append(f"累积频率不是严格递增: cumfreq[{i}]={cumfreq[i]} >= cumfreq[{i+1}]={cumfreq[i+1]}")
    
    return errors


def test_m2_exhaustive():
    """M=2符号的详尽测试"""
    print("="*70)
    print("M=2符号详尽测试")
    print("="*70)
    
    encoder = StandardArithmeticEncoder()
    decoder = StandardArithmeticDecoder()
    
    import random
    random.seed(42)
    
    TOT = StandardArithmeticEncoder.TOT
    num_trials = 10000
    passed = 0
    failed_cases = []
    
    print(f"运行 {num_trials} 次测试...")
    
    for trial in range(num_trials):
        # 随机序列长度 (0到500)
        seq_len = random.randint(0, 500)
        
        # 随机概率分割: p in {1..TOT-1}, [p, TOT-p]
        p = random.randint(1, TOT - 1)
        prob0 = p / TOT
        prob1 = (TOT - p) / TOT
        probs = [prob0, prob1]
        
        # 生成随机序列
        symbols = [random.randint(0, 1) for _ in range(seq_len)]
        
        try:
            # 编码
            bitstream = encoder.encode(symbols, probs)
            
            # 验证编码器的不变量（在编码过程中）
            # 这里我们无法直接访问内部状态，所以跳过
            
            # 解码
            decoded_symbols = decoder.decode(bitstream, probs)
            
            # 验证
            if decoded_symbols != symbols:
                # 找到第一个不匹配的位置
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
                
                # 只保存前10个失败案例
                if len(failed_cases) <= 10:
                    print(f"\n❌ 测试 {trial} 失败:")
                    print(f"  序列长度: {seq_len}")
                    print(f"  概率: {probs} (p={p})")
                    print(f"  输入符号 (前20个): {symbols[:20]}")
                    print(f"  解码符号 (前20个): {decoded_symbols[:20]}")
                    print(f"  位流长度: {len(bitstream)}")
                    print(f"  第一个不匹配位置: {first_mismatch}")
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
                print(f"\n❌ 测试 {trial} 异常: {e}")
    
    print(f"\n" + "="*70)
    print(f"测试结果:")
    print(f"  总测试数: {num_trials}")
    print(f"  通过: {passed}")
    print(f"  失败: {len(failed_cases)}")
    print(f"  成功率: {passed/num_trials*100:.2f}%")
    
    if failed_cases:
        print(f"\n失败案例详情 (前{min(5, len(failed_cases))}个):")
        for i, case in enumerate(failed_cases[:5]):
            print(f"\n案例 {i+1}:")
            print(f"  测试编号: {case['trial']}")
            print(f"  序列长度: {case['seq_len']}")
            print(f"  概率: {case['probs']} (p={case['p']})")
            if 'error' in case:
                print(f"  错误: {case['error']}")
            else:
                print(f"  输入符号: {case['symbols']}")
                print(f"  解码符号: {case['decoded']}")
                print(f"  位流长度: {case['bitstream_len']}")
                print(f"  第一个不匹配位置: {case['first_mismatch']}")
                
                # 打印最小可复现案例
                print(f"\n  最小可复现案例:")
                print(f"    probs = {case['probs']}")
                print(f"    symbols = {case['symbols']}")
                print(f"    bitstream_len = {case['bitstream_len']}")
                if case['first_mismatch'] is not None:
                    print(f"    first_mismatch_index = {case['first_mismatch']}")
    
    print("="*70)
    
    assert len(failed_cases) == 0, f"M=2测试失败: {len(failed_cases)}个案例失败"
    print("✓ M=2详尽测试全部通过！")
    
    return len(failed_cases) == 0


if __name__ == "__main__":
    try:
        # 先运行标准测试
        print("运行标准测试...")
        success = test_arithmetic_codec()
        if not success:
            exit(1)
        
        # 运行M=2详尽测试
        print("\n" + "="*70)
        print("运行M=2详尽测试...")
        m2_success = test_m2_exhaustive()
        if not m2_success:
            exit(1)
        
        print("\n" + "="*70)
        print("🎉 所有测试通过！")
        print("="*70)
        
    except AssertionError as e:
        print(f"\n❌ 断言失败: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
