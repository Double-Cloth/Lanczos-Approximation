# Lanczos Approximation

[中文](README.md) | [English](README.en.md)

一个基于任意精度 BigInt/BigFloat 的 Lanczos Gamma 计算项目，包含：

- 系数生成与落盘
- 复用系数的单点 Gamma 求值
- 基于 CSV 的批量误差验证
- BigInt/BigFloat 单元测试

## 目录结构

```text
Lanczos-Approximation/
├─ include/
│  ├─ BigInt.h
│  ├─ BigFloat.h
│  └─ lanczos.h
├─ src/
│  ├─ BigInt.cpp
│  ├─ BigFloat.cpp
│  ├─ lanczos.cpp
│  └─ main.cpp
├─ tests/
│  └─ test_bigfloat.cpp
├─ assets/
│  └─ real_gamma.csv
└─ CMakeLists.txt
```

## 环境要求

- CMake >= 3.15
- 支持 C++17 的编译器（GCC/Clang/MSVC）
- 可选：OpenMP（若检测到会自动启用）

## 构建

### 通用（推荐）

```bash
cmake -B build
cmake --build build
```

### Windows MinGW

```bash
cmake -B build -G "MinGW Makefiles"
cmake --build build
```

### Windows MSVC

```bash
cmake -B build -G "Visual Studio 17 2022"
cmake --build build --config Release
```

构建后主要产物：

- lanczos_app：主程序
- test_bigfloat：单元测试

## 主程序用法

主程序有 3 种模式。

### 1) 生成模式（Generate）

计算系数，写出结果文件，并执行 CSV 验证（默认只验证 z <= 50，测试点数量不设上限）。

```bash
./build/lanczos_app <n> <g> <digits> [--out dir] [--csv path] [--auto-upgrade] [--hex]
```

参数说明：

- n：级数最高索引（正整数）
- g：Lanczos 参数（字符串形式读取）
- digits：十进制精度位数（正整数）
- --out：输出目录（默认 output_n<n>_g<g>_d<digits>）
- --csv：验证数据文件（默认 ../assets/real_gamma.csv）
- --auto-upgrade：根据 digits 自动升级 n/g 到推荐档位
- --hex：将 coefficients.txt 中的系数按十六进制科学计数法写出（`0x1.ffffp+N`）

示例：

```bash
./build/lanczos_app 20 21.5 80 --out output_n20_g21.5_d80

# 16进制系数输出
./build/lanczos_app 20 21.5 80 --hex
```

### 2) 求值模式（Eval）

从已有输出目录或单文件系数文本读取参数与系数，计算单点 Gamma(z)。

```bash
./build/lanczos_app eval <output_dir_or_file> <z_value> [display_digits] [--hex]
```

参数说明：

- output_dir_or_file：
  - 目录模式：目录下需有 parameters.txt 与 coefficients.txt
  - 单文件模式：支持包含 state / approx coef 段的文本
- z_value：要计算的 z
- display_digits：可选，显示位数；不传时程序自动选择
- --hex：可选，按十六进制科学计数法输出（格式类似 `0x1.ffffp+N`）

示例：

```bash
./build/lanczos_app eval output_n20_g21.5_d80 50.5 40

# 十六进制输出
./build/lanczos_app eval output_n20_g21.5_d80 50.5 40 --hex
```

### 3) 测试模式（Test）

计算系数后执行批量验证，可做顺序或随机抽样。

```bash
./build/lanczos_app test <n> <g> <digits> [--csv path] [--max N] [--start row] [--random] [--threshold %] [--auto-upgrade]
```

参数说明：

- --csv：CSV 路径（默认 ../assets/real_gamma.csv）
- --max：最大测试点数（默认 50）
- --start：起始行（1-based，默认 1）
- --random：随机抽样（否则顺序取样）
- --threshold：相对误差百分比阈值（默认 1e-6）
- --auto-upgrade：根据 digits 自动升级 n/g

示例：

```bash
./build/lanczos_app test 20 21.5 80 --start 100 --max 200 --random --threshold 1e-10
```

## 输出文件说明

生成模式会在输出目录写入：

- coefficients.txt：系数列表与公式注释
- parameters.txt：参数与精度元信息
- verification.txt：批量验证结果与统计摘要

parameters.txt 关键字段：

- n
- g
- coefficient_count
- precision_decimal_digits
- precision_binary_bits
- coefficient_dump_digits

## 单元测试

项目已接入 CTest。

```bash
cd build
ctest --output-on-failure
```

或直接运行：

```bash
./build/test_bigfloat
```

test_bigfloat 支持详细日志：

```bash
# 命令行开关
./build/test_bigfloat --verbose

# 或环境变量
TEST_VERBOSE=1 ./build/test_bigfloat
```

默认行为是只输出每个测试组结果与耗时；失败时自动打印详细误差上下文。

## 精度与参数建议

常用档位（与程序内推荐逻辑一致）：

- digits <= 7：n=6, g=5.581
- digits <= 16：n=13, g=13.144565
- digits <= 19：n=17, g=17.0
- 更高精度：n=24, g=23.5

说明：

- 你可以手动指定任意 n/g。
- 开启 --auto-upgrade 时，程序会在 n 偏小时自动提升到推荐档位。

## 误差判定

验证采用相对误差百分比：

$$
\text{relative error \%} = \frac{|\text{computed} - \text{expected}|}{|\text{expected}|}
$$

- 测试模式：阈值由 --threshold 指定。
- 生成模式：阈值由程序根据 digits 与 n 自动计算。

## 常见问题

1. 无参数运行直接退出

- 这是正常行为，程序会打印用法并返回非 0。

2. eval 报错 invalid or missing parameters found

- 检查输入目录是否包含 parameters.txt 与 coefficients.txt。
- 若是单文件模式，确认包含 state 与 approx coef 区段。

3. ctest 提示 DartConfiguration.tcl 缺失

- 一般不影响测试执行结果；属于 CDash 相关配置提示。

## 许可证

见 LICENSE。
