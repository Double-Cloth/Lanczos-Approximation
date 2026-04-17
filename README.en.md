# Lanczos Approximation

[中文](README.md) | [English](README.en.md)

A Lanczos Gamma computation project based on arbitrary-precision BigInt/BigFloat, including:

- Coefficient generation and persistence
- Single-point Gamma evaluation by reusing generated coefficients
- Batch error verification based on CSV data
- BigInt/BigFloat unit tests

## Project Structure

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

## Requirements

- CMake >= 3.15
- A C++17-compatible compiler (GCC/Clang/MSVC)
- Optional: OpenMP (enabled automatically if detected)

## Build

### Cross-platform (recommended)

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

Main artifacts after build:

- lanczos_app: main executable
- test_bigfloat: unit test executable

## Main Program Usage

The main program supports 3 modes.

### 1) Generate Mode

Computes coefficients, writes output files, and runs CSV verification
(by default only checks z <= 50, with no upper bound on the number of test points).

```bash
./build/lanczos_app <n> <g> <digits> [--out dir] [--csv path] [--auto-upgrade] [--hex]
```

Arguments:

- n: highest series index (positive integer)
- g: Lanczos parameter (read as string)
- digits: decimal precision digits (positive integer)
- --out: output directory (default: output_n<n>_g<g>_d<digits>)
- --csv: verification data file (default: ../assets/real_gamma.csv)
- --auto-upgrade: auto-upgrade n/g to recommended levels based on digits
- --hex: write coefficients in coefficients.txt using hexadecimal scientific notation (`0x1.ffffp+N`)

Example:

```bash
./build/lanczos_app 20 21.5 80 --out output_n20_g21.5_d80

# Hex coefficient output
./build/lanczos_app 20 21.5 80 --hex
```

### 2) Eval Mode

Loads parameters and coefficients from an existing output directory or a single coefficient text file, then evaluates Gamma(z).

```bash
./build/lanczos_app eval <output_dir_or_file> <z_value> [display_digits] [--hex]
```

Arguments:

- output_dir_or_file:
  - Directory mode: must contain parameters.txt and coefficients.txt
  - Single-file mode: supports text containing state / approx coef sections
- z_value: target z
- display_digits: optional display precision; if omitted, the program chooses automatically
- --hex: optional, print result in hexadecimal scientific notation (format like `0x1.ffffp+N`)

Example:

```bash
./build/lanczos_app eval output_n20_g21.5_d80 50.5 40

# Hex output
./build/lanczos_app eval output_n20_g21.5_d80 50.5 40 --hex
```

### 3) Test Mode

Generates coefficients and runs batch verification, either sequentially or randomly sampled.

```bash
./build/lanczos_app test <n> <g> <digits> [--csv path] [--max N] [--start row] [--random] [--threshold %] [--auto-upgrade]
```

Arguments:

- --csv: CSV path (default: ../assets/real_gamma.csv)
- --max: max number of test points (default: 50)
- --start: start row (1-based, default: 1)
- --random: random sampling (otherwise sequential)
- --threshold: relative error percentage threshold (default: 1e-6)
- --auto-upgrade: auto-upgrade n/g based on digits

Example:

```bash
./build/lanczos_app test 20 21.5 80 --start 100 --max 200 --random --threshold 1e-10
```

## Output Files

Generate mode writes the following files in the output directory:

- coefficients.txt: coefficient list and formula notes
- parameters.txt: parameter and precision metadata
- verification.txt: batch verification results and summary stats

Key fields in parameters.txt:

- n
- g
- coefficient_count
- precision_decimal_digits
- precision_binary_bits
- coefficient_dump_digits

## Unit Tests

CTest is integrated in this project.

```bash
cd build
ctest --output-on-failure
```

Or run directly:

```bash
./build/test_bigfloat
```

test_bigfloat supports verbose logs:

```bash
# CLI flag
./build/test_bigfloat --verbose

# Or environment variable
TEST_VERBOSE=1 ./build/test_bigfloat
```

Default behavior prints only each test group result and timing; detailed error context is printed automatically on failures.

## Precision and Parameter Suggestions

Common presets (aligned with the built-in recommendation logic):

- digits <= 7: n=6, g=5.581
- digits <= 16: n=13, g=13.144565
- digits <= 19: n=17, g=17.0
- higher precision: n=24, g=23.5

Notes:

- You can manually specify any n/g.
- With --auto-upgrade enabled, the program automatically raises small n values to recommended levels.

## Error Metric

Verification uses relative error percentage:

$$
\text{relative error \%} = \frac{|\text{computed} - \text{expected}|}{|\text{expected}|}
$$

- Test mode: threshold is set by --threshold.
- Generate mode: threshold is computed automatically from digits and n.

## FAQ

1. Program exits when run without arguments

- This is expected. The program prints usage and returns non-zero.

2. eval reports invalid or missing parameters found

- Check whether the input directory contains parameters.txt and coefficients.txt.
- For single-file mode, make sure state and approx coef sections exist.

3. ctest reports missing DartConfiguration.tcl

- Usually does not affect test execution; this is a CDash-related configuration notice.

## License

See LICENSE.
