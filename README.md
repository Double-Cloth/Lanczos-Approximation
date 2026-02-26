# Lanczos Approximation Coefficient Calculator

一个纯 C++17 实现的 Lanczos 逼近系数计算器，内置任意精度算术库，无需任何第三方依赖。

## 简介

[Lanczos 逼近](https://en.wikipedia.org/wiki/Lanczos_approximation) 是一种高效计算 Gamma 函数的数值方法。本程序通过 **Godfrey 矩阵方法** (`P = D × B × C × F`) 计算 Lanczos 系数，并使用这些系数来高精度求值 Gamma 函数。

矩阵公式严格遵循 [Boost.Math](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/lanczos.html) 的实现。

## 项目结构

```
Lanczos-Approximation/
├── CMakeLists.txt             # 构建配置
├── include/
│   ├── BigInt.h               # 无符号任意精度整数
│   ├── BigFloat.h             # 任意精度浮点数
│   └── lanczos.h              # Lanczos 系数计算接口
├── src/
│   ├── BigInt.cpp             # BigInt 实现
│   ├── BigFloat.cpp           # BigFloat 实现 (含 pi, sqrt, exp, ln, pow 等)
│   ├── lanczos.cpp            # Godfrey 矩阵方法 + Gamma 函数
│   └── main.cpp               # 命令行应用
└── tests/
    ├── test_bigfloat.cpp      # BigInt/BigFloat 单元测试（相对误差断言）
    └── test_lanczos.cpp       # Lanczos 精度验证测试（相对误差百分比）
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

编译成功后会在 `build/` 目录下生成三个可执行文件：

| 文件 | 说明 |
|------|------|
| `lanczos_app` | 主程序 — 计算系数并验证 |
| `test_bigfloat` | BigInt/BigFloat 单元测试（含相对误差 % 断言） |
| `test_lanczos` | Lanczos 精度验证（需 `real_gamma.csv`，相对误差 % 判定） |

## 使用方法

### 命令行参数

程序目前提供两种运行模式：

#### 模式 1：生成验证模式
计算 Lanczos 系数并运行全量比较验证。
```bash
./build/lanczos_app <n> <g> <digits> [output_dir] [csv_path]
```
| 参数 | 类型 | 说明 |
|------|------|------|
| `n` | 整数 | 系数个数（级数项数） |
| `g` | 浮点数 | Lanczos 参数 |
| `digits` | 整数 | 目标精度（十进制有效位数） |
| `output_dir` | 字符串（可选） | 输出目录名，默认为 `output_n<n>_g<g>_d<digits>` |
| `csv_path` | 字符串（可选） | 验证数据 CSV 文件路径，默认为 `../assets/real_gamma.csv` |

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

### 示例

#### 1. 标准精度 (~10 位有效数字)

```bash
./build/lanczos_app 7 5.0 16
```

输出：

```
=== Computed Coefficients ===
  p[0] = 2.5066282751072970e+0
  p[1] = 1.9095517189307640e+2
  p[2] = -2.1683668184372796e+2
  p[3] = 6.0194417640233329e+1
  p[4] = -3.0875132392854592e+0
  p[5] = 3.0296387052532575e-3
  p[6] = -1.3523859590726280e-5

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
#
0, 2.5066282751072970e+0
1, 1.9095517189307640e+2
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
# Total tests:        100
# Passed (relative error <= 1e-08%): 100/100
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

# Lanczos 精度验证 (需要 real_gamma.csv 文件在工作目录)
# 默认阈值: 相对误差 ≤ 1e-6% (约 8 位有效数字)
./build/test_lanczos

# 自定义参数: test_lanczos [n] [g] [precision] [csv_path] [max_tests] [threshold%]
./build/test_lanczos 7 5.0 16 real_gamma.csv 50 1e-8
```

#### 测试通过标准

所有测试使用 **相对误差百分比** 作为判定标准：

$$\text{relative\_error} = \frac{|\text{computed} - \text{expected}|}{|\text{expected}|} \times 100\%$$

| 测试文件 | 测试类型 | 阈值 | 含义 |
|---------|---------|------|------|
| `test_bigfloat` | 基本算术 | ≤ 1e-30% | 应几乎精确 |
| `test_bigfloat` | 常量（π、√2、e） | ≤ 1e-45% | ~47 位有效数字 |
| `test_bigfloat` | 复合运算（exp、ln） | ≤ 1e-40% | ~42 位有效数字 |
| `test_bigfloat` | Γ 函数 | ≤ 1e-30% | ~32 位有效数字 |
| `test_lanczos` | Gamma 函数验证 | ≤ 1e-6%（默认） | ~8 位有效数字 |
| `lanczos_app` | CSV 全量验证 | 动态计算 | `10^(-min(digits, n+3) + 2) %` |

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

### Lanczos 层与测试框架层
| 优化 | 说明 |
|------|------|
| **F[i] 构造常数分离** | Godfrey 方法中将 `e^{-g}` 吸收进系数生成，运行时消除 `exp` |
| **F[i] exp 复用** | `exp(power_{i+1}) = exp(power_i) × exp(1)`，将降维的高精指数计算进一步递推 |
| **运行时求值极速化** | 单个 Gamma 以一次大数幂运算 `pow(base/e, z-0.5)` 完成，避开高成本 `exp(-base)` |
| **pow2 指数调整** | 除以 2^k 改用 `mul_pow2(-k)`（O(1) 指数调整） |
| **高精验证动态截断** | 根据对比需求在转数字前动态截断 CSV 中 3 万位级别的极限验证数据，避免阻塞瓶颈 |

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

## 算法原理

本程序使用 **Godfrey 矩阵方法** 计算 Lanczos 系数：

$$P = D \times B \times C \times F$$

其中：
- **B** — 由二项式系数构成的矩阵，`B[i][j] = (-1)^{j-i} × C(i+j-1, j-i)`
- **D** — 对角矩阵，通过递推 `D[i] = D[i-1] × 2(2i-1)/(i-1)` 计算
- **C** — Chebyshev 系数矩阵，`C[i][j] = (-1)^{i-j} × Σ C(2i,2k)×C(k,k+j-i)`
- **F** — 包含纯量指数、阶乘和幂函数项的浮点向量，其通式已被优化吸收常数 $e^{-g}$：
  $$F_i = \frac{(2i)!}{i!} \cdot \exp(i+0.5) \cdot \frac{1}{2^{2i-1} (g+i+0.5)^i \sqrt{g+i+0.5}}$$

计算得到系数 `p_0, ..., p_{n-1}` 后，Gamma 函数通过以下公式求值：

$$\Gamma(z) = \left( \frac{z + g - 0.5}{e} \right)^{z - 0.5} \cdot S(z)$$

其中级数部分为：

$$S(z) = p_0 + \sum_{k=1}^{n-1} \frac{p_k}{z + k - 1}$$

> **注意：** 本实现中 `√(2π)` 因子已吸收在系数 `p_0` 中，无需额外乘以。

## 参考文献

- Lanczos, C. (1964). "A Precision Approximation of the Gamma Function". *SIAM Journal on Numerical Analysis*.
- Godfrey, P. ["A note on the computation of the convergent Lanczos complex Gamma approximation"](http://my.fit.edu/~gabdo/gamma.txt).
- Pugh, G. (2004). "An Analysis of the Lanczos Gamma Approximation". PhD Thesis, University of British Columbia.
- Boost.Math — [Lanczos Approximation Documentation](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/lanczos.html).
- Robert Munafo — [Lanczos Gamma Coefficients](https://mrob.com/pub/ries/lanczos-gamma.html).

## 许可

MIT License
