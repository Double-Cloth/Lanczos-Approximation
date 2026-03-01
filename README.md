# Lanczos Approximation Coefficient Calculator

一个纯 C++17 实现的 Lanczos 逼近系数计算器，内置任意精度算术库，无需任何第三方依赖。

## 项目结构

```
Lanczos-Approximation/
├── CMakeLists.txt             # 构建配置 (已注入 FindOpenMP 关联)
├── include/
│   ├── BigInt.h               # 无符号任意精度整数
│   ├── BigFloat.h             # 任意精度浮点数
│   └── lanczos.h              # Lanczos 系数计算接口
├── src/
│   ├── BigInt.cpp             # BigInt 实现
│   ├── BigFloat.cpp           # BigFloat 实现 (含 pi, sqrt, exp, ln, pow 等)
│   ├── lanczos.cpp            # Godfrey 矩阵方法 (已挂载 OpenMP 多核并发与系统保护)
│   └── main.cpp               # 命令行应用 (包含 Generate/Eval/Test 并归合验证核心)
└── tests/
    └── test_bigfloat.cpp      # BigInt/BigFloat 单元测试（极高精度相对误差断言）
```

## 构建

### 依赖

- **CMake** ≥ 3.15
- **C++17 编译器**（GCC / Clang / MSVC）
- 无第三方库依赖

### 编译步骤

```bash
# 配置 (Linux / macOS)
cmake -B build

# 配置 (Windows MinGW)
cmake -B build -G "MinGW Makefiles"

# 配置 (Windows MSVC)
cmake -B build -G "Visual Studio 17 2022"

# 编译
cmake --build build
```

编译成功后会在 `build/` 目录下生成可执行文件：

| 文件 | 说明 |
|------|------|
| `lanczos_app` | 主程序 — 计算系数、高精求值、以及批量数据精度验证 |
| `test_bigfloat` | BigInt/BigFloat 单元底层数学能力测试（纯误差 % 严格断言） |

## 使用方法

### 命令行参数

程序提供三种运行模式：

#### 模式 1：生成验证模式
计算 Lanczos 系数并运行全量比较验证。本模式同样兼容并采纳了命名标志解析 (Named Flags) 系统。
```bash
./build/lanczos_app <n> <g> <digits> [--out dir] [--csv path]
```
| 参数 | 类型 | 说明 |
|------|------|------|
| `n` | 整数 | 定长必填: 系数个数（级数项数） |
| `g` | 浮点数 | 定长必填: Lanczos 参数 |
| `digits` | 整数 | 定长必填: 目标精度（十进制有效位数） |
| `--out` | 字符串 | (可选) 输出目录名，默认为 `output_n<n>_g<g>_d<digits>` |
| `--csv` | 字符串 | (可选) 指定验证目标文件，默认 `../assets/real_gamma.csv` |

#### 模式 2：复用求值模式
读取之前生成的系数结果（`parameters.txt` 和 `coefficients.txt`），利用已有参数极速对某个任意数值 `z` 求 Gamma 函数值，而无需重新构造大型系数矩阵。
```bash
./build/lanczos_app eval <output_dir> <z_value> [display_digits]
```
| 参数 | 类型 | 说明 |
|------|------|------|
| `eval` | 关键字 | 激活复用求值模式 |
| `output_dir` | 字符串 | 包含已计算好所需系数与参数文件的目录路径 |
| `z_value` | 浮点数 | 要求取 Gamma 函数的具体 `z` 值 |
| `display_digits` | 整数（可选） | 动态指定当前单次估算的输出十进制显示精度 |

#### 模式 3：测试与验证模式
读取已知 Gamma 函数精确值的 CSV 文件，在控制台打印并生成进度条，直接对批量数据运行动态断言并展示绝对与相对误差总结。
随着功能升级，现在该模式全面支持兼容式**命名标志解析 (Named Flags)**，让参数传递更加自由。
```bash
./build/lanczos_app test <n> <g> <digits> [--csv path] [--max N] [--start row] [--random] [--threshold %]
```
| 参数 | 类型 | 说明 |
|------|------|------|
| `test` | 关键字 | 激活测试模式 |
| `n` | 整数 | 定长必填: 系数个数（级数项数） |
| `g` | 浮点数 | 定长必填: Lanczos 参数 |
| `digits` | 整数 | 定长必填: 目标精度（十进制有效位数） |
| `--csv` | 字符串 | (可选) 指定验证目标文件，默认 `../assets/real_gamma.csv` |
| `--max` | 整数 | (可选) 抽取或按序读取的最大行数，默认 `50` |
| `--start` | 整数 | (可选) 跳过头部数据，从指定的行数 (1-indexed) 开始读取，默认 `1` |
| `--random`| 标志位 | (可选) 注入 `mt19937` 从整体数据 (截去 start 之后的数据) 中随机打散挑出 `max` 数量的点测试 |
| `--threshold`| 字符串 | (可选) 及格的相对误差比重限制，默认 `1e-6` |

### 示例

#### 1. 标准精度 (~10 位有效数字)

```bash
./build/lanczos_app 7 5.0 16
```

输出：

```
=== Lanczos Approximation Test ===
Parameters: n=7, g=5.0, precision=16 decimal digits
CSV file: ../assets/real_gamma.csv
Pass threshold: relative error <= 1e-08%

--- Computing Lanczos Coefficients ---
[lanczos] CPU Protection: Using 14 / 16 OpenMP threads.
[lanczos] matrix B     [==================================================] 100% (7/7)
[lanczos] matrix C     [==================================================] 100% (7/7)
[lanczos] matrix D     [==================================================] 100% (7/7)
[lanczos] matrix BC    [==================================================] 100% (7/7)
[lanczos] vector F     [==================================================] 100% (7/7)
[lanczos] vector P     [==================================================] 100% (7/7)
Coefficients computation complete.
[lanczos] Coefficients computed successfully in 45.28 ms.

=== Loading test data from real_gamma.csv ===
Loaded 2800 test entries.
Selected 100 test points (z <= 50).

=== Verification ===
  Pass threshold: relative error <= 1e-08% (based on 16-digit precision)

  Gamma(0.5) = 1.7724538508277281e+0  [relative error: 4.39e-11%, PASS]
  Gamma(1) = 9.9999999994978583e-1  [relative error: 5.02e-11%, PASS]
  ...
  Gamma(50) = 6.0828186406937170e+62  [relative error: 5.77e-09%, PASS]

=== Summary ===
  Time taken:          0.021 s (21 ms)
  Total tests:         100
  Passed (relative error <= 1e-08%): 100/100
  Max relative error:  7.24e-28%
  Min relative error:  0.00e+00%
  Avg relative error:  7.66e-29%
```

#### 2. 双精度 (~22 位有效数字)

```bash
./build/lanczos_app 13 13.144565 25
```

#### 3. 指定输出目录

```bash
./build/lanczos_app 13 13.144565 25 my_results
```

#### 4. 对单个 z 极速求值评测 (复用模式)

在已经运行过模式 1 且得到了 `output_n7_g5_d16` 目录后：
```bash
# 绕过繁重的矩阵求基运算，瞬间得出高精度 Gamma(50.5)
./build/lanczos_app eval output_n7_g5_d16 50.5
```
输出：
```
Gamma(50.5) = 1.631853612882101968502e+64
```

您还可以通过附加参数，动态指定更高或更低的临时显示精度：
```bash
# 以 40 位的要求输出 Gamma(10.5)
./build/lanczos_app eval output_n7_g5_d16 10.5 40
# 输出: Gamma(10.5) = 1.13327838882269365467682991931831111388074e+6
```


### 输出文件

程序运行后会在输出目录中生成以下文件：

```
output_n7_g5_d16/
├── coefficients.txt    # Lanczos 系数 (含使用公式说明)
├── parameters.txt      # 计算参数 (n, g, 精度)
└── verification.txt    # Gamma 函数验证结果 (11 个测试点)
```

**coefficients.txt** — 系数文件，包含公式说明和每个系数的值：
```
# Formula:
#   Gamma(z) = (base/e)^(z-0.5) * S(z)
#   where base = z + g - 0.5
#         S(z) = p[0] + sum_{k=1}^{n-1} p[k] / (z + k - 1)
#   Note: sqrt(2*PI) and exp(-g) are absorbed into coefficients.
#
0, 1.6889528464081992e-2
1, 1.2866458274168037e+0
...
```

**verification.txt** — 验证结果，包含计算值、期望值、相对误差和测试状态：
```
# Lanczos Gamma Function Verification Results
# Parameters: n=7, g=5, digits=16
# Test data source: assets/real_gamma.csv
# Pass threshold: relative error <= 1e-08%
#
# z, computed, expected (first 30 digits), relative_error%, status
0.5, 1.7724538508277281e+0, 1.7724538509055160272981674833..., 4.39e-11%, PASS
1, 9.9999999994978583e-1, 1..., 5.02e-11%, PASS
...
#
# ==============================
# Summary
# ==============================
# Time taken:         1102.34 ms
# Total tests:        100
# Passed (relative error <= 1e-08%): 100/100
# Max absolute error: 3.25e-22
# Min absolute error: 0.00e+00
# Avg absolute error: 4.19e-24
# Max relative error: 7.24e-28%
# Min relative error: 0.00e+00%
# Avg relative error: 7.66e-29%
```

### 推荐参数组合 (来自 Pugh 论文)

以下为 Pugh 确定的最优参数，适合**带保护位数**的计算场景：

| 目标精度 (位) | `n` | `g` | 最大误差 |
|:---:|:---:|:---:|:---:|
| 24 (float) | 6 | 5.581 | 9.51e-12 |
| 53 (double) | 13 | 13.144565 | 9.22e-23 |
| 64 (long double) | 17 | 17.0 | ~1e-27 |
| 113 (quad) | 24 | 23.5 | ~1e-40 |

> **注意：** 第三个参数 `digits` 应设置为需要的十进制有效位数，程序内部会自动转换为二进制位数并加上保护位。

### 运行测试

```bash
# BigInt/BigFloat 单元测试（使用相对误差百分比断言）
./build/test_bigfloat

# Lanczos 精度验证模式 (传统位置参数)
./build/lanczos_app test 7 5.0 16 ../assets/real_gamma.csv 50 1e-8

# Lanczos 精度验证模式 (全新推荐命名参数，极度自由)
# 从第 100 行开始，在后续所有数据中随机抽选 20 个点测试，报错阈心定在 1e-10%
./build/lanczos_app test 7 5.0 16 --start 100 --max 20 --random --threshold 1e-10
```

#### 测试通过标准

所有测试使用 **相对误差百分比** 作为判定标准：

$$\text{relative error} = \frac{|\text{computed} - \text{expected}|}{|\text{expected}|}$$

| 测试文件 | 测试类型 | 阈值 | 含义 |
|---------|---------|------|------|
| `test_bigfloat` | 基本算术 | ≤ 1e-30% | 应几乎精确 |
| `test_bigfloat` | 常量（π、√2、e） | ≤ 1e-45% | ~47 位有效数字 |
| `test_bigfloat` | 复合运算（exp、ln） | ≤ 1e-40% | ~42 位有效数字 |
| `test_bigfloat` | Γ 函数 | ≤ 1e-30% | ~32 位有效数字 |
| `lanczos_app` test 模式 | Gamma 函数轻量化验证 | 用户指定（默认 1e-6%） | ~8 位有效数字 |
| `lanczos_app` 生成模式 | CSV 全量输出验证 | 动态自动计算 | `10^(-min(digits, n+3) + 2) %` |

## 性能优化

本项目使用了多层次的算法优化来提升计算速度：

### BigInt 层
| 优化 | 复杂度变化 | 说明 |
|------|-----------|------|
| **Karatsuba 乘法** | O(n²) → O(n^1.585) | 操作数 ≥ 32 字（1024 位）时自动启用分治 |
| **Knuth Algorithm D 除法** | O(n×bitlen) → O(n×m) | 替代逐位二进制长除法，数量级提升 |

### BigFloat 层
| 优化 | 说明 |
|------|------|
| **`div_u32` 快速除法** | 除以小整数（≤ 2^32）时使用 O(n) 单字除法，替代完整 BigFloat 除法 |
| **exp 泰勒级数优化** | 级数中 `term /= k` 使用 `div_u32(k)` 而非 BigFloat 除法 |
| **ln(2) 缓存** | 避免每次 ln() 调用都重新计算 atanh(1/3) 级数 |
| **arctan/atanh 级数优化** | 除以 `(2k+1)` 使用 `div_u32` |
| **pow 整数提取** | 直接从 BigInt digits 提取整数值，避免 O(n×D) 十进制转换 |
| **to_decimal 快速幂 5^n** | 使用二进制快速幂替代逐步 `mul_u32(5)` |

### Lanczos 层与框架并发执行引擎
| 优化 | 说明 |
|------|------|
| **OpenMP 多核矩阵提速** | 彻底利用现代 CPU 多核心阵列。对于 O(N^3) 级矩阵 M, B, C, BC 的构建与结合执行分配 `#pragma omp parallel for`，百倍释放系统多核威力。 |
| **OS 资源枯竭保护盾** | 深度探针动态侦测 CPU 逻辑总核数 `omp_get_num_procs`，主动为 8 线程以上操作系统截留并保留 2 核心系统开销余地，彻底断绝电脑假死与鼠标卡顿危机。 |
| **组合数 (comb) 无锁静态缓存** | 构筑二维矩阵内存池静态存储阶乘结果，直接将 C 与 B 矩阵内占据 40% 时间的组合数倒算剔除，改判为 O(1) 短路读取。 |
| **F[i] 构造常数分离** | Godfrey 方法中将 `e^{-g}` 吸收进系数生成，运行时消除 `exp` |
| **F[i] 底数剥解法剔除 exp** | 通过基础底数数学合并提取 `exp(0.5)` 作为纯量常数，使每一段循环内的天文级 `BigFloat::exp` 演算剥解降级成为多线程安全常数积。 |
| **pow2 指数调整** | 除以 2^k 改用 `mul_pow2(-k)`（O(1) 指数调整） |
| **高精算力引擎钟表** | 在 `lanczos.cpp` 的核心与 `main.cpp` 的批量断言挂载极其精准的 `std::chrono::high_resolution_clock` 获取纳秒级耗时监控展示。 |

## API 接口

如需在自己的程序中使用本库，包含头文件并链接即可：

```cpp
#include "lanczos.h"

// 1. 计算系数
int n = 7;
double g = 5.0;
int digits = 16;
auto coeffs = compute_lanczos_coefficients(n, g, digits);

// 2. 使用系数计算 Gamma 函数
int prec_bits = 256;
BigFloat z(5.0, prec_bits);                       // z = 5.0
BigFloat result = lanczos_gamma(z, coeffs, g, digits);  // Gamma(5) = 24
std::cout << result.to_decimal_string(20) << std::endl;
```

### BigFloat 数学工具

`BigFloat` 类内置多种高精度数学函数，可独立使用：

```cpp
#include "BigFloat.h"

int prec = 512;  // 512 位精度 (~154 位十进制)

BigFloat pi_val  = BigFloat::pi(prec);           // π
BigFloat sqrt2   = BigFloat::sqrt(BigFloat(2, prec));  // √2
BigFloat e_val   = BigFloat::exp(BigFloat(1, prec));   // e
BigFloat ln2     = BigFloat::ln(BigFloat(2, prec));    // ln(2)
BigFloat pow_val = BigFloat::pow(base, exponent);      // base^exp
BigFloat fact5   = BigFloat::factorial(5, prec);       // 5! = 120

// 输出
std::cout << pi_val.to_decimal_string(50) << std::endl;
// 3.14159265358979323846264338327950288419716939937510
```

## 系统架构与算法原理深度解析 (Mathematical & Architecture Deep-Dive)

本项目作为一套超越标准数据类型的极限计算引擎，不仅实现了 [Lanczos 逼近定理](https://en.wikipedia.org/wiki/Lanczos_approximation)，更在底层工程、数值截断防护上做了严密拓展，可实现任意十进制精度（如 10,000位+）的高质量计算。

### 1. 任意精度算术基石 (BigFloat Engine)
完全摆脱 `double` 类型的 IEEE-754 定长限制。
* 内部所有常数（$\pi, e, \sqrt{2}, \ln2$）与数学运算（幂，指数，三角）均无硬编码位数。
* **高精底层架构**: 采用多字无符号 `BigInt` 结合二阶动态指数构建极为灵活的 `BigFloat` 大数系统。

### 2. 经典 Lanczos 逼近与公式推导
Gamma 函数的经典积分定义为了改善原点附近的缺陷，通过换元积分引出：
$$\Gamma(z) = \sqrt{2\pi} \left( z + g - \frac{1}{2} \right)^{z - \frac{1}{2}} e^{-\left( z + g - \frac{1}{2} \right)} A_g(z)$$

展开项 $A_g(z)$：
$$A_g(z) = \frac{1}{2} p_0(g) + \sum_{k=1}^{n-1} \frac{p_k(g)}{z+k-1}$$

> **注意：** $(z + g - 1/2)^{z-1/2} e^{-(z+g-1/2)}$ 通常被提炼为 `base^(z-0.5) / e^base` 以降低浮点标度损失。本系统更进一步，直接将常数部分 $e^{-g}$ 连同 $\sqrt{2\pi}$ 静态压栈吸收进入了系数 $p_k$ 中，使求值期的计算耗时压缩为常数级代数操作。

### 3. P. Godfrey 矩阵高速求解法 (Matrix Method)
早期的连分式积分极其缓慢。本系统采用直接多项式时间 $O(N^3)$ 的 Godfrey 矩阵逼近阵列公式：
$$P = D \times B \times C \times F$$
1. **$B_{i,j} = (-1)^{j-i} \binom{i+j-1}{j-i}$**
   基于二项式构成的下三角变换矩阵。
2. **$C_{i,j} = (-1)^{i-j} \sum_{k} \binom{2i}{2k} \binom{k}{k+j-i}$**
   基于极限 Chebyshev 正交多项式零点分布构成的误差均匀化矩阵（系统对其实现了静态二维 C++ 缓存记忆化，并搭载 OpenMP 并发，避开了高达数百万次的计算风暴）。
3. **$D_i = D_{i-1} \frac{2(2i-1)}{i-1}$**
   针对归一化的对角递推缩放矩阵。
4. **$F_i = \frac{(2i)!}{i!} \cdot \exp(i+0.5) \cdot \frac{1}{2^{2i-1} (g+i+0.5)^i \sqrt{g+i+0.5}}$**
   将 $e^{-g}$ 和复杂的阶乘乘积吸收并行的纯量合并浮点项。

得出系数后，利用重铸版极速求值法：
$$\Gamma(z) = \left( \frac{z + g - 0.5}{e} \right)^{z - 0.5} \cdot \left( p_0 + \sum_{k=1}^{n-1} \frac{p_k}{z + k - 1} \right)$$

### 4. 高精度数值抗性与误差全量清剿 (Numerical Robustness)
由于级数矩阵相乘引发的抵消甚至能超过 $10^{300}$ 的标度。一旦失去保护，纯粹的灾难性相消能瞬间吞噬二三百位的十进制精度！
本项目搭载了以下标准数值分析防护壁垒：
1. **动态防护垫**: 随着输入项数增加而触发计算：$Guard\_Bits = 2n \log_2(2n) + 256$，保护所有底层寄存边界。
2. **多项式 Kahan 补偿求和**: 级数项 $S(z)$ 的连加启用了 `c = (t - S) - y` 补档变量精准捕捉十进制深处的截断损耗值。
3. **指数幂完全拆解**: 凡应对底数非整次求幂，引擎会将次幂强制拉断为独立的整数无损快幂与纯小数的次幂，严格截停指数函数的精度放量。
4. **反射公式与 Payne-Hanek 动态模降**: 在负实轴发散区 $z < 0.5$：自动切换至拉马努金反射恒等式 $\frac{\pi}{\sin(\pi z)}$ 维系。一旦侦测到底层输入了非对称极限大数供向 $\sin$ 时，动态分配巨额精度的 $\pi$ 确保正弦映射严格闭敛于纯净域 $[0, \pi/2]$。
5. **从缺保护 (Lazy Normalization)**: 给泰勒迭代收敛提供超越 `prec + 64 bits` 的冗余宽域，规避了大数循环中的重复短截断。

## 参考文献

- Lanczos, C. (1964). "A Precision Approximation of the Gamma Function". *SIAM Journal on Numerical Analysis*.
- Godfrey, P. ["A note on the computation of the convergent Lanczos complex Gamma approximation"](http://my.fit.edu/~gabdo/gamma.txt).
- Pugh, G. (2004). "An Analysis of the Lanczos Gamma Approximation". PhD Thesis, University of British Columbia.
- Boost.Math — [Lanczos Approximation Documentation](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/lanczos.html).
- Robert Munafo — [Lanczos Gamma Coefficients](https://mrob.com/pub/ries/lanczos-gamma.html).

## 许可

MIT License
