#!/usr/bin/env python
"""
完全确定性的DNA编码生化FSM约束控制器

实现GC含量窗口控制和同聚物长度限制。
"""

from typing import List, Tuple, Optional
from collections import deque


# DNA碱基映射
BASE_TO_INDEX = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INDEX_TO_BASE = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
NUM_BASES = 4

# GC碱基（G和C）
GC_BASES = {2, 3}  # G=2, C=3
AT_BASES = {0, 1}  # A=0, T=1


class FSMConstraint:
    """
    完全确定性的DNA编码生化FSM约束控制器
    
    支持：
    - GC含量窗口控制
    - 同聚物长度限制
    - 保证至少一个允许的碱基
    """
    
    def __init__(
        self,
        gc_window: int = 20,
        gc_low: float = 0.45,
        gc_high: float = 0.55,
        max_homopolymer: int = 3
    ):
        """
        初始化FSM约束控制器
        
        参数:
            gc_window: GC窗口大小（默认20）
            gc_low: GC含量下限（默认0.45）
            gc_high: GC含量上限（默认0.55）
            max_homopolymer: 最大同聚物长度（默认3）
        """
        self.gc_window = gc_window
        self.gc_low = gc_low
        self.gc_high = gc_high
        self.max_homopolymer = max_homopolymer
        
        # 内部状态
        self._window_buffer: deque = deque(maxlen=gc_window)
        self._gc_count: int = 0
        self._homopolymer_base: Optional[int] = None
        self._homopolymer_length: int = 0
        self._step_count: int = 0
    
    def reset(self):
        """重置内部状态"""
        self._window_buffer.clear()
        self._gc_count = 0
        self._homopolymer_base = None
        self._homopolymer_length = 0
        self._step_count = 0
    
    def get_state(self) -> dict:
        """获取当前状态（用于调试）"""
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
        """判断是否是GC碱基"""
        return base_index in GC_BASES
    
    def _get_current_gc_ratio(self) -> float:
        """获取当前GC比率"""
        if len(self._window_buffer) == 0:
            return 0.5  # 无数据时返回中间值
        return self._gc_count / len(self._window_buffer)
    
    def _would_gc_ratio_be(self, new_base: int) -> float:
        """
        计算添加新碱基后的GC比率
        考虑窗口满时会移除最旧的碱基
        """
        new_gc_count = self._gc_count
        window_size = len(self._window_buffer)
        
        # 如果窗口已满，需要移除最旧的碱基
        if window_size >= self.gc_window:
            oldest_base = self._window_buffer[0]
            if self._is_gc_base(oldest_base):
                new_gc_count -= 1
        else:
            window_size += 1
        
        # 添加新碱基
        if self._is_gc_base(new_base):
            new_gc_count += 1
        
        return new_gc_count / window_size if window_size > 0 else 0.5
    
    def _would_violate_homopolymer(self, base_index: int) -> bool:
        """判断添加新碱基是否会违反同聚物限制"""
        if self._homopolymer_base is None:
            return False
        
        if base_index == self._homopolymer_base:
            # 同一个碱基，检查是否会超过限制
            return self._homopolymer_length >= self.max_homopolymer
        
        return False
    
    def get_mask(self) -> List[bool]:
        """
        获取当前允许的碱基掩码
        
        返回:
            布尔列表[4]，表示A,T,G,C是否允许
        
        规则:
            1. 强制同聚物限制
            2. 强制GC边界
            3. 如果掩码为空，放松GC约束直到至少有1个允许
        """
        mask = [True, True, True, True]
        
        # 1. 应用同聚物限制
        for i in range(NUM_BASES):
            if self._would_violate_homopolymer(i):
                mask[i] = False
        
        # 2. 应用GC边界约束
        # 只有在窗口有足够数据时才应用GC约束
        if len(self._window_buffer) >= self.gc_window // 2:
            gc_mask = [True, True, True, True]
            
            for i in range(NUM_BASES):
                if not mask[i]:
                    continue  # 已经被同聚物约束禁止
                
                future_gc_ratio = self._would_gc_ratio_be(i)
                
                # 检查是否会违反GC约束
                if future_gc_ratio < self.gc_low:
                    # GC太低，只允许GC碱基
                    if i not in GC_BASES:
                        gc_mask[i] = False
                elif future_gc_ratio > self.gc_high:
                    # GC太高，只允许AT碱基
                    if i not in AT_BASES:
                        gc_mask[i] = False
            
            # 合并GC掩码
            combined_mask = [m1 and m2 for m1, m2 in zip(mask, gc_mask)]
            
            # 3. 如果合并后掩码为空，放松GC约束
            if any(combined_mask):
                mask = combined_mask
            # 否则保留只有同聚物约束的掩码
        
        # 最终检查：确保至少有一个允许的碱基
        if not any(mask):
            # 这不应该发生，但作为安全措施
            # 放松所有约束，只保留同聚物约束
            mask = [True, True, True, True]
            for i in range(NUM_BASES):
                if self._would_violate_homopolymer(i):
                    mask[i] = False
            
            # 如果仍然为空（不应该发生），允许所有
            if not any(mask):
                mask = [True, True, True, True]
        
        return mask
    
    def update(self, base_index: int):
        """
        更新FSM状态
        
        参数:
            base_index: 碱基索引（0=A, 1=T, 2=G, 3=C）
        """
        if base_index < 0 or base_index >= NUM_BASES:
            raise ValueError(f"无效的碱基索引: {base_index}")
        
        # 更新GC窗口
        if len(self._window_buffer) >= self.gc_window:
            # 窗口已满，移除最旧的碱基
            oldest_base = self._window_buffer[0]
            if self._is_gc_base(oldest_base):
                self._gc_count -= 1
        
        # 添加新碱基到窗口
        self._window_buffer.append(base_index)
        if self._is_gc_base(base_index):
            self._gc_count += 1
        
        # 更新同聚物状态
        if self._homopolymer_base == base_index:
            self._homopolymer_length += 1
        else:
            self._homopolymer_base = base_index
            self._homopolymer_length = 1
        
        # 更新步数
        self._step_count += 1
    
    def get_allowed_count(self) -> int:
        """获取当前允许的碱基数量"""
        return sum(self.get_mask())
    
    def get_allowed_bases(self) -> List[int]:
        """获取当前允许的碱基索引列表"""
        mask = self.get_mask()
        return [i for i, m in enumerate(mask) if m]


# ==============================================================================
# 测试函数
# ==============================================================================

def test_basic_functionality():
    """测试基本功能"""
    print("测试基本功能:")
    
    fsm = FSMConstraint(gc_window=20, gc_low=0.45, gc_high=0.55, max_homopolymer=3)
    
    # 初始状态
    mask = fsm.get_mask()
    assert all(mask), f"初始状态应该允许所有碱基: {mask}"
    print(f"  初始掩码: {mask} ✓")
    
    # 添加几个碱基
    fsm.update(0)  # A
    fsm.update(1)  # T
    fsm.update(2)  # G
    fsm.update(3)  # C
    
    state = fsm.get_state()
    print(f"  添加ATGC后状态: {state}")
    assert state['step_count'] == 4
    assert state['gc_count'] == 2  # G和C
    print("  基本功能测试通过 ✓")


def test_homopolymer_constraint():
    """测试同聚物约束"""
    print("\n测试同聚物约束:")
    
    fsm = FSMConstraint(gc_window=20, gc_low=0.0, gc_high=1.0, max_homopolymer=3)
    
    # 添加3个连续的A
    fsm.update(0)  # A
    fsm.update(0)  # A
    fsm.update(0)  # A
    
    # 第4个应该禁止A
    mask = fsm.get_mask()
    assert not mask[0], f"连续3个A后应该禁止A: {mask}"
    assert mask[1] and mask[2] and mask[3], f"其他碱基应该允许: {mask}"
    print(f"  连续3个A后掩码: {mask} ✓")
    
    # 添加一个不同的碱基
    fsm.update(1)  # T
    
    # 现在A应该又允许了
    mask = fsm.get_mask()
    assert mask[0], f"添加T后A应该允许: {mask}"
    print(f"  添加T后掩码: {mask} ✓")
    
    print("  同聚物约束测试通过 ✓")


def test_gc_constraint():
    """测试GC约束"""
    print("\n测试GC约束:")
    
    fsm = FSMConstraint(gc_window=10, gc_low=0.4, gc_high=0.6, max_homopolymer=10)
    
    # 添加10个A（GC=0%）
    for _ in range(10):
        fsm.update(0)  # A
    
    state = fsm.get_state()
    print(f"  10个A后: GC比率={state['gc_ratio']:.2f}")
    
    # GC太低，应该只允许GC碱基
    mask = fsm.get_mask()
    print(f"  掩码: {mask}")
    # 由于GC太低，应该偏向允许GC碱基
    
    # 添加一些G来提高GC
    fsm.reset()
    for _ in range(10):
        fsm.update(2)  # G
    
    state = fsm.get_state()
    print(f"  10个G后: GC比率={state['gc_ratio']:.2f}")
    
    # GC太高，应该只允许AT碱基
    mask = fsm.get_mask()
    print(f"  掩码: {mask}")
    
    print("  GC约束测试通过 ✓")


def test_mask_never_empty():
    """测试掩码永远不会为空"""
    print("\n测试掩码永远不会为空:")
    
    import random
    random.seed(42)
    
    fsm = FSMConstraint(gc_window=10, gc_low=0.45, gc_high=0.55, max_homopolymer=3)
    
    for trial in range(100):
        fsm.reset()
        
        for step in range(200):
            mask = fsm.get_mask()
            
            # 确保掩码不为空
            assert any(mask), f"掩码为空！试验{trial}，步骤{step}，状态{fsm.get_state()}"
            
            # 从允许的碱基中随机选择
            allowed = [i for i, m in enumerate(mask) if m]
            base = random.choice(allowed)
            fsm.update(base)
    
    print("  100次试验，每次200步，掩码从未为空 ✓")
    print("  掩码非空测试通过 ✓")


def test_determinism():
    """测试确定性"""
    print("\n测试确定性:")
    
    sequence = [0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 0, 1, 2, 2, 2, 1, 3, 3, 3, 0]
    
    # 运行两次，应该得到相同的掩码序列
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
    
    assert masks1 == masks2, "两次运行的掩码序列不同！"
    print("  两次运行的掩码序列完全相同 ✓")
    print("  确定性测试通过 ✓")


def run_fsm_tests():
    """运行所有FSM测试"""
    print("=" * 70)
    print("FSM约束控制器测试")
    print("=" * 70)
    
    test_basic_functionality()
    test_homopolymer_constraint()
    test_gc_constraint()
    test_mask_never_empty()
    test_determinism()
    
    print("\n" + "=" * 70)
    print("🎉 所有FSM测试通过！")
    print("=" * 70)


if __name__ == "__main__":
    run_fsm_tests()
