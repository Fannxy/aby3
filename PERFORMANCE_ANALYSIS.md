# ABY3 Join 和 Group 操作性能分析报告

## 一、监测数据分析

### 1. 小数据场景（<10行）
- **Join Test**: 执行时长约80秒，CPU平均3.19%，最大12.50%
- **Group Test**: 执行时长相对较短

### 2. 大数据场景（1024行，q13）
- **执行时长**: 显著增加（预计数倍到数十倍）
- **性能下降原因**: 多个操作的复杂度随数据量非线性增长

## 二、Join操作性能瓶颈分析

### 2.1 主要步骤
1. **augment_table** (1062-1150行)
2. **oblivious_expand** (1830-2000+行) - 对T1和T2各执行一次
3. **align_table** (2427-2500+行)
4. **数据重组** (2747-2763行)
5. **bool2arith** (2768行)
6. **memcpy操作** (2773-2784行)

### 2.2 关键性能瓶颈

#### 瓶颈1: augment_table中的逐元素复制
```cpp
// 1086-1101行: 嵌套循环，逐元素复制
for(size_t i=0; i<len1; i++){
    for(size_t j=0; j<key_num; j++){
        T_c[0].mShares[0](i, j) = T_1_key[j].mShares[0](i, 0);
        T_c[0].mShares[1](i, j) = T_1_key[j].mShares[1](i, 0);
    }
    // ... 更多逐元素复制
}
```
**问题**: 复杂度O(n*m)，当n=1024时，需要大量逐元素访问
**优化建议**: 
- 使用整列赋值：`T_c[0].col(j) = T_1_key[j].col(0);`（如果T_1_key[j]是单列）
- 使用memcpy批量复制整列数据（如果数据布局允许）
- 重新组织数据结构，使列优先存储

#### 瓶颈2: oblivious_expand中的循环加密操作
```cpp
// 1863-1878行: 对每一行都执行加密操作
for(size_t i=0; i<n; i++){
    // ...
    bool_cipher_selector(pIdx, flag_i, zero_i, s, f_i, enc, eval, runtime);
    bool_cipher_add(pIdx, s, g_i, s, enc, eval, runtime);
}
```
**问题**: 
- 每次循环都触发网络通信（三方安全计算）
- 当n=1024时，需要1024次网络往返
- 加密操作本身开销大

**优化建议**:
- 批量处理多行数据，减少网络往返次数
- 使用SIMD指令优化本地计算部分
- 考虑使用更高效的加密协议

**遍历累加的具体优化思路**（不修改主代码，仅作方案参考）：

循环语义可归纳为：
- `s` 为累积和（64-bit），初值 0；
- 每轮：`f_i = flag_i ? 0 : s`（selector），`fx[i]=f_i`，然后 `s = s + g_i`。
即 `s_i = g[0]+…+g[i-1]`（前缀和），`fx[i] = (1-flag[i]) * s_i`。瓶颈在于每轮 2 次网络调用、共 2n 轮。

**思路 1：批量前缀和 + 批量 selector（推荐）**

1. **一次求出所有前缀和**  
   用「批量前缀和」电路/协议：输入 `g(n×64)`，输出 `s_prefix(0..n-1)`，其中 `s_prefix[i] = g[0]+…+g[i-1]`。  
   - 若底层有现成的前缀和（或 scan）原语，直接调用一次即可。  
   - 若无，可用一个电路表达递推：`s[0]=0, s[i+1]=s[i]+g[i]`，在 Sh3BinaryEvaluator 里用一条「大电路」一次跑完，把 2n 轮压成 1 轮（或很少几轮）。

2. **批量算 fx**  
   `fx[i] = (1-flag[i]) * s_prefix[i]`。  
   - 把 `flag(n×1)`、`s_prefix(n×64)`、`zero(n×64)` 拼成矩阵，用**支持多行的** `bool_cipher_selector`（Basics.h 里已有 `sbMatrix &flag` 多行版本），一次调用得到整列 `fx`，不再按行循环。

示例伪代码（仅说明接口形态，不要求你改现有代码）：

```cpp
// 伪代码：批量版
sbMatrix g(n, 64);           // 已有
sbMatrix flag(n, 1);        // 已有
sbMatrix s_prefix(n, 64);    // 待求：s_prefix[i] = g[0]+...+g[i-1]

// 1) 一次批量前缀和（需在电路/协议层实现或调用现有 scan）
bool_prefix_sum(pIdx, g, s_prefix, enc, eval, runtime);  // 1 轮

// 2) fx[i] = (1-flag[i]) * s_prefix[i]，批量 selector
sbMatrix zero_n(n, 64);
bool_init_false(pIdx, zero_n);
sbMatrix fx(n, 64);
bool_cipher_selector(pIdx, flag, zero_n, s_prefix, fx, enc, eval, runtime);  // 1 轮

// 3) 最终 s 就是 s_prefix 的最后一行（若需要）
s.mShares[0](0,0) = s_prefix.mShares[0](n-1, 0);
s.mShares[1](0,0) = s_prefix.mShares[1](n-1, 0);
// 若协议里 s = s_prefix[n-1] + g[n-1]，再补一次 add 即可
```

这样从 **2n 次网络往返** 降为 **2 次**（前缀和 1 次 + 批量 selector 1 次）。

**思路 2：分块折中（不改协议时的最小改动）**

若暂时不加「批量前缀和」电路，可把 n 行分成块长 B（如 32）：

- 块内仍用原循环算 `s` 和 `fx`，但每块只做 B 轮；
- 块与块之间只传「块末的 s」到下一块。

轮数从 2n 变为 2*(n/B)（块间同步）+ 块内 2B，若 B 固定则约为 O(n/B) 次同步。收益小于思路 1，但实现简单。

**思路 3：树状递推（电路层）**

在电路里显式实现递推链：  
`S[0]=0, S[i+1]=S[i]+G[i], FX[i]=Mux(1-flag[i], 0, S[i])`。  
把整条链做成一个电路，一次 `asyncEvaluate` 跑完，则通信轮次与电路深度相关（通常仍可视为常数轮），而不是 n。需要在前端用 CircuitLibrary 或自定义电路描述「多级加法 + Mux」，并接到现有 bool_cipher_add / selector 的底层原语。

**小结**

| 方案           | 通信轮次     | 实现难度 | 备注                     |
|----------------|--------------|----------|--------------------------|
| 现状（逐行）   | O(n)         | -        | 简单但慢                 |
| 思路 1 批量    | O(1)         | 中       | 需前缀和 + 批量 selector |
| 思路 2 分块    | O(n/B)       | 低       | 不改协议即可用           |
| 思路 3 单电路  | O(1)/O(depth)| 高       | 电路设计与维护成本大     |

实际落地时优先考虑**思路 1**：在 GORAM-Core/Basics 或单独模块里实现 `bool_prefix_sum`（或复用已有 scan），再配合多行 `bool_cipher_selector` 即可在不改上层业务逻辑的前提下，把这段遍历累加从 O(n) 轮压到常数轮。

#### 瓶颈3: 数据重组的三重嵌套循环
```cpp
// 2747-2763行: 三重嵌套循环
for(size_t i=0; i<m; i++){
    for(size_t l=0; l<key_num; l++){
        T_joined_sb_concat.mShares[0](i+l*m, 0) = T_1_expanded[0].mShares[0](i, l);
        // ...
    }
    // ... 更多嵌套循环
}
```
**问题**: 复杂度O(m*cols)，逐元素访问效率低
**优化建议**:
- 使用memcpy批量复制整行或整列
- 重新组织数据布局，减少转置操作

#### 瓶颈4: bool2arith在大矩阵上的开销
```cpp
// 2768行: 对整个连接矩阵执行bool2arith
bool2arith(pIdx, T_joined_sb_concat, T_joined_si_concat, enc, eval, runtime);
```
**问题**: bool2arith在大型矩阵上非常耗时，涉及大量加密操作
**优化建议**:
- 如果可能，减少需要转换的数据量
- 优化bool2arith的实现，使用批量处理

## 三、Group操作性能瓶颈分析

### 3.1 主要步骤
1. **数据拼接** (570-577行)
2. **genPerm** (578行) - 生成排列
3. **applyPerm** (584-585行) - 应用排列（两次）
4. **compare_consecutive_rows_arith** (591行)
5. **逐列AND操作** (613-624行)
6. **bool2arith** (627行)
7. **更多排列操作** (697, 701行)

### 3.2 关键性能瓶颈

#### 瓶颈1: 数据拼接的逐元素复制
```cpp
// 570-577行: 嵌套循环，逐元素复制
for(size_t i=0; i<rows; i++){
    k_v_concat.mShares[0](i, 0) = val.mShares[0](i, 0);
    k_v_concat.mShares[1](i, 0) = val.mShares[1](i, 0);
    for(size_t j=1; j<=key_cols; j++){
        k_v_concat.mShares[0](i, j) = key.mShares[0](i, key_cols-j);
        k_v_concat.mShares[1](i, j) = key.mShares[1](i, key_cols-j);
    }
}
```
**问题**: 复杂度O(rows*cols)，逐元素访问
**优化建议**:
- 使用memcpy批量复制
- 如果key列需要倒序，考虑一次性memcpy然后reverse

#### 瓶颈2: 多次排列操作
```cpp
// 578行: genPerm
genPerm(pIdx, k_v_concat, perm, enc, eval, runtime);
// 584-585行: applyPerm（两次）
applyPerm(pIdx, perm, key, key_g, enc, eval, runtime);
applyPerm(pIdx, perm, val, val_g, enc, eval, runtime);
// 697行: 再次genPerm
genPerm(pIdx, tmp, perm_GN, enc, eval, runtime);
// 701行: 再次applyPerm
applyPerm(pIdx, perm_GN, key_GN, key_out, enc, eval, runtime);
```
**问题**: 
- genPerm和applyPerm都需要网络通信
- 当rows=1024时，排列操作的开销很大
- 总共执行了3次genPerm和3次applyPerm

**优化建议**:
- 合并排列操作，减少网络往返
- 优化排列算法的实现
- 考虑使用更高效的排列网络

#### 瓶颈3: 逐列AND操作
```cpp
// 613-624行: 对每一列都执行AND操作
for(int j=0;j<key_cols;j++){
    // 提取列数据（逐元素）
    for(int i=0;i<rows-1;i++){
        f_vector_col.mShares[0](i, 0) = f_vector.mShares[0](i*key_cols+j, 0);
        f_vector_col.mShares[1](i, 0) = f_vector.mShares[1](i*key_cols+j, 0);
    }
    // 执行AND操作（需要网络通信）
    bool_cipher_and(pIdx, f, f_vector_col, f_new, enc, eval, runtime);
    f = f_new;
}
```
**问题**: 
- 如果有k列，就要执行k次bool_cipher_and
- 每次AND操作都需要网络通信
- 列提取也是逐元素操作

**优化建议**:
- 批量处理多列，减少网络往返
- 优化列提取，使用memcpy（如果数据布局允许）
- 考虑使用SIMD指令并行处理多列

#### 瓶颈4: compare_consecutive_rows_arith
```cpp
// 591行: 比较连续行
compare_consecutive_rows_arith(pIdx, key_g, f_vector, enc, eval, runtime);
```
**问题**: 当rows=1024时，需要比较1023对行，涉及大量加密操作
**优化建议**:
- 优化比较算法，使用批量比较
- 减少不必要的加密操作

## 四、为什么q13（1024行）比小数据慢很多？

### 4.1 复杂度分析
- **小数据（<10行）**: 
  - 大部分操作的复杂度可以视为O(1)或O(n)，其中n很小
  - 网络通信次数少
  - 加密操作次数少

- **大数据（1024行）**:
  - 许多操作是O(n²)或更高复杂度
  - 网络通信次数随n线性或平方增长
  - 加密操作次数大幅增加

### 4.2 具体影响
1. **网络通信开销**: 
   - 小数据: 每次操作通信量小，延迟占主导
   - 大数据: 通信量大幅增加，带宽和延迟都成为瓶颈

2. **加密计算开销**:
   - 小数据: 加密操作次数少
   - 大数据: 加密操作次数成倍增加（如oblivious_expand中的1024次循环）

3. **内存访问模式**:
   - 小数据: 缓存友好
   - 大数据: 缓存不友好，大量逐元素访问导致缓存未命中

## 五、优化建议总结

### 5.1 立即可以实施的优化（高优先级）

#### 1. 使用整列赋值和memcpy替代逐元素复制
**位置**: 
- `augment_table` (1086-1111行) - 使用 `T_c[0].col(j) = T_1_key[j].col(0);`
- `group_by_common` 数据拼接 (570-577行)
- `join` 数据重组 (2747-2763行)
- `group_by_common` 列提取 (616-619行)

**方法**: 
- 识别可以批量复制的数据块
- 使用整列赋值（如果数据结构支持）
- 使用memcpy替代嵌套循环
- 注意数据对齐和大小

#### 2. 批量处理加密操作
**位置**:
- `oblivious_expand` (1863-1878行)
- `group_by_common` 逐列AND (613-624行)

**方法**:
- 将多个独立的加密操作合并为批量操作
- 减少网络往返次数
- 使用向量化的加密操作

#### 3. 优化数据布局
**方法**:
- 使用列优先存储，减少转置操作
- 预分配内存，减少动态分配
- 使用内存池管理临时矩阵

### 5.2 中期优化（中优先级）

#### 1. 优化排列操作
**方法**:
- 合并多个排列操作
- 使用更高效的排列网络算法
- 减少排列操作的网络通信

#### 2. 优化bool2arith
**方法**:
- 批量处理多个矩阵
- 优化底层加密协议
- 使用SIMD指令优化本地计算

#### 3. 减少不必要的加密操作
**方法**:
- 识别可以本地计算的部分
- 延迟加密，批量加密
- 使用更轻量级的加密方案

### 5.3 长期优化（低优先级）

#### 1. 算法级优化
**方法**:
- 重新设计join和group算法，减少复杂度
- 使用分块处理，避免一次性处理所有数据
- 考虑使用更高效的MPC协议

#### 2. 并行化
**方法**:
- 识别可以并行执行的操作
- 使用多线程处理独立的数据块
- 优化网络通信的并行性

#### 3. 缓存优化
**方法**:
- 优化内存访问模式，提高缓存命中率
- 使用预取指令
- 重新组织数据结构

## 六、具体优化代码示例

### 示例1: 优化augment_table中的数据复制
```cpp
// 原代码（1086-1101行）: 逐元素复制
// 优化后: 使用整列赋值（如果T_1_key[j]是单列矩阵）
for(size_t j=0; j<key_num; j++){
    // 如果 T_1_key[j] 是单列矩阵（cols() == 1）
    T_c[0].col(j) = T_1_key[j].col(0);  // Eigen 会优化整列赋值
}

// 或者使用memcpy批量复制（如果数据布局允许）
for(size_t j=0; j<key_num; j++){
    // 复制T_1_key[j]的整列到T_c[0]的第j列
    // 注意：由于RowMajor布局，列数据不连续，需要逐行复制或使用Eigen的列操作
    for(size_t i=0; i<len1; i++){
        T_c[0].mShares[0](i, j) = T_1_key[j].mShares[0](i, 0);
        T_c[0].mShares[1](i, j) = T_1_key[j].mShares[1](i, 0);
    }
    // 或者使用 Eigen 的列操作（更高效）
    T_c[0].mShares[0].col(j).head(len1) = T_1_key[j].mShares[0].col(0).head(len1);
    T_c[0].mShares[1].col(j).head(len1) = T_1_key[j].mShares[1].col(0).head(len1);
}
```

### 示例2: 优化group_by_common中的列提取
```cpp
// 原代码（616-619行）: 逐元素提取
// 优化后: 使用Eigen的列操作
for(int j=0;j<key_cols;j++){
    sbMatrix f_vector_col(rows-1, 1);
    // 使用Eigen的列操作提取（注意f_vector是展平的）
    // 由于f_vector是展平的，需要特殊处理
    for(int i=0;i<rows-1;i++){
        f_vector_col.mShares[0](i, 0) = f_vector.mShares[0](i*key_cols+j, 0);
        f_vector_col.mShares[1](i, 0) = f_vector.mShares[1](i*key_cols+j, 0);
    }
    // 如果f_vector可以重新组织为矩阵，可以使用更高效的列提取
}
```

### 示例3: 优化join中的数据重组
```cpp
// 原代码（2747-2763行）: 三重嵌套循环
// 优化后: 使用Eigen的列操作
for(size_t l=0; l<key_num; l++){
    // 使用Eigen的列操作
    T_joined_sb_concat.mShares[0].col(0).segment(l*m, m) = 
        T_1_expanded[0].mShares[0].col(l);
    T_joined_sb_concat.mShares[1].col(0).segment(l*m, m) = 
        T_1_expanded[0].mShares[1].col(l);
}
```

## 七、性能预期

### 7.1 优化前（当前状态）
- 小数据（<10行）: ~80秒
- 大数据（1024行）: 预计数百秒到数千秒

### 7.2 优化后预期
- **使用整列赋值和memcpy优化**: 预计提升20-30%
- **批量处理加密操作**: 预计提升30-50%
- **优化排列操作**: 预计提升10-20%
- **综合优化**: 预计总体提升50-70%

### 7.3 注意事项
- 优化效果取决于具体的数据分布和网络条件
- 某些优化可能需要修改底层加密库
- 需要充分测试确保正确性
- Eigen的列操作在RowMajor布局下虽然不如行操作快，但仍比逐元素访问高效

## 八、监控和验证

### 8.1 性能监控
- 继续使用现有的监测脚本
- 添加更细粒度的时间测量
- 监控网络通信量

### 8.2 正确性验证
- 使用小数据验证优化后的正确性
- 逐步增加数据量进行测试
- 对比优化前后的结果

## 九、总结

q13（1024行）比小数据慢很多的主要原因是：

1. **大量的逐元素复制操作**，导致缓存不友好和CPU效率低
2. **频繁的网络通信**，每次加密操作都需要三方通信
3. **加密操作本身的开销**，在大数据量时成倍增加
4. **算法复杂度高**，许多操作是O(n²)或更高

**最有效的优化方向**:
1. 使用整列赋值和memcpy批量复制（立即实施，效果明显）
2. 批量处理加密操作（中期实施，效果显著）
3. 优化排列操作（中期实施，减少网络通信）
4. 算法级优化（长期实施，根本性改进）

建议优先实施整列赋值和批量处理加密操作，这两项优化相对简单但效果显著。

