/**
 * @file test_bigfloat.cpp
 * @brief BigInt 和 BigFloat 基本功能的单元测试
 *
 * @details 测试覆盖:
 *   - BigInt: 加法、乘法（含大数）、除法取模、位移
 *   - BigFloat: 加减乘除、分数表示精度
 *   - 数学函数: π、√2、e、exp、ln、sqrt(π)、Γ(0.5)、阶乘
 *
 *   所有测试通过 assert 断言验证，失败时程序会终止并报告错误位置。
 *   数学函数测试使用 **相对误差百分比** 作为通过标准：
 *     |计算值 − 期望值| / |期望值| × 100% ≤ 阈值%
 *   计算结果同时输出到控制台供人工比对已知值。
 */

#include "BigFloat.h"
#include "BigInt.h"
#include <cassert>

#include <iostream>
#include <string>

// ========================
// 辅助函数
// ========================

/**
 * @brief 计算 BigFloat 的相对误差百分比
 * @param computed  计算值
 * @param expected  期望值
 * @return 相对误差百分比（0.0 表示完全匹配）
 *
 * @details 计算公式: |computed - expected| / |expected| × 100%
 *   若 expected 为零，则 computed 也必须为零才算通过。
 *   通过将差值和期望值都转为高精度十进制字符串，
 *   再用 double 计算比值来得到百分比。
 */
double relative_error_percent(const BigFloat &computed,
                              const BigFloat &expected) {
  BigFloat diff = (computed - expected).abs();

  if (expected.is_zero()) {
    return computed.is_zero() ? 0.0 : 100.0;
  }

  // 使用足够多的有效数字来精确计算误差比值
  // 将 diff / |expected| 转换为十进制字符串再解析为 double
  BigFloat abs_expected = expected.abs();
  BigFloat ratio = diff / abs_expected;
  std::string ratio_str = ratio.to_decimal_string(20, true);
  double ratio_val = std::stod(ratio_str);

  return ratio_val * 100.0;
}

/**
 * @brief 检查计算值与期望值的相对误差是否在允许范围内
 * @param computed       计算值
 * @param expected       期望值
 * @param max_percent    允许的最大相对误差百分比
 * @param test_name      测试名称（用于输出）
 * @return true 如果通过，false 如果失败
 *
 * @details 输出格式:
 *   ✓ test_name: 相对误差 = x.xxe-yy% (≤ 阈值%) [PASS]
 *   ✗ test_name: 相对误差 = x.xxe-yy% (> 阈值%) [FAIL]
 */
bool check_relative_error(const BigFloat &computed, const BigFloat &expected,
                          double max_percent, const std::string &test_name) {
  double err = relative_error_percent(computed, expected);
  bool pass = (err <= max_percent);

  std::cout << "  " << (pass ? "[PASS]" : "[FAIL]") << " " << test_name
            << ": relative error = " << std::scientific << err << "%"
            << " (threshold: " << max_percent << "%)" << std::defaultfloat
            << std::endl;

  return pass;
}

/**
 * @brief BigInt 基本运算测试
 *
 * 测试项:
 *   - 加法: 12345 + 67890 = 80235
 *   - 乘法: 12345 × 67890 = 838102050
 *   - 大数乘法: 999999999999999999² = 999999999999999998000000000000000001
 *   - 除法: 100 ÷ 7 = 14 余 2
 *   - 左移: 1 << 100 的 bit_length 应为 101
 */
bool test_bigint_basic() {
  std::cout << "--- BigInt Basic Tests ---" << std::endl;

  BigInt a(12345);
  BigInt b(67890);

  // 加法测试
  BigInt c = a + b;
  assert(c.to_decimal_string() == "80235");
  std::cout << "  [PASS] 12345 + 67890 = " << c << std::endl;

  // 乘法测试
  BigInt d = a * b;
  assert(d.to_decimal_string() == "838102050");
  std::cout << "  [PASS] 12345 * 67890 = " << d << std::endl;

  // 大数乘法: 验证 18 位 × 18 位 = 36 位的正确性
  BigInt e("999999999999999999");
  BigInt f("999999999999999999");
  BigInt g = e * f;
  assert(g.to_decimal_string() == "999999999999999998000000000000000001");
  std::cout << "  [PASS] 999999999999999999^2 = " << g << std::endl;

  // 除法: 验证商和余数
  BigInt h(100);
  BigInt i(7);
  auto [q, r] = h.divmod(i);
  assert(q.to_decimal_string() == "14");
  assert(r.to_decimal_string() == "2");
  std::cout << "  [PASS] 100 / 7 = " << q << " remainder " << r << std::endl;

  // 位移: 1 << 100 应是 101 位
  BigInt j(1);
  BigInt k = j << 100;
  assert(k.bit_length() == 101);
  std::cout << "  [PASS] 1 << 100: bit_length = " << k.bit_length()
            << std::endl;

  return true;
}

/**
 * @brief BigFloat 算术运算测试
 *
 * 测试项:
 *   - 加法: 3 + 4 = 7
 *   - 减法: 7 - 4 = 3
 *   - 乘法: 3 × 4 = 12
 *   - 除法: 3 / 4 = 0.75, 1 / 3 = 0.333...（验证循环小数精度）
 *
 * 通过标准: 相对误差 ≤ 1e-30%（基本算术应几乎精确）
 */
bool test_bigfloat_arithmetic() {
  std::cout << std::endl << "--- BigFloat Arithmetic Tests ---" << std::endl;
  int prec = 128; // 128 位精度（约 38 位十进制）
  bool all_pass = true;

  // 加法: 3 + 4 = 7
  BigFloat a(3, prec);
  BigFloat b(4, prec);
  BigFloat c = a + b;
  BigFloat expected_sum(7, prec);
  std::cout << "  computed: 3 + 4 = " << c.to_decimal_string(10) << std::endl;
  all_pass &= check_relative_error(c, expected_sum, 1e-30, "3 + 4 = 7");

  // 减法: 7 - 4 = 3
  BigFloat sub_result = c - b;
  std::cout << "  computed: 7 - 4 = " << sub_result.to_decimal_string(10)
            << std::endl;
  all_pass &= check_relative_error(sub_result, a, 1e-30, "7 - 4 = 3");

  // 乘法: 3 * 4 = 12
  BigFloat d = a * b;
  BigFloat expected_prod(12, prec);
  std::cout << "  computed: 3 * 4 = " << d.to_decimal_string(10) << std::endl;
  all_pass &= check_relative_error(d, expected_prod, 1e-30, "3 * 4 = 12");

  // 除法: 3 / 4 = 0.75
  BigFloat e = a / b;
  BigFloat expected_div(0.75, prec);
  std::cout << "  computed: 3 / 4 = " << e.to_decimal_string(20) << std::endl;
  all_pass &= check_relative_error(e, expected_div, 1e-30, "3 / 4 = 0.75");

  // 除法: 1 / 3 = 0.333...（无限循环小数，验证 * 3 ≈ 1）
  BigFloat one(1, prec);
  BigFloat three(3, prec);
  BigFloat third = one / three;
  BigFloat reconstructed = third * three;
  std::cout << "  computed: 1 / 3 = " << third.to_decimal_string(40)
            << std::endl;
  std::cout << "  computed: (1/3) * 3 = " << reconstructed.to_decimal_string(40)
            << std::endl;
  // 循环小数允许更宽松的误差: 1e-25%（约 27 位有效数字）
  all_pass &= check_relative_error(reconstructed, one, 1e-25, "(1/3) * 3 ≈ 1");

  assert(all_pass);
  return all_pass;
}

/**
 * @brief BigFloat 数学函数测试
 *
 * 测试项（使用 256 位精度，约 77 位十进制）:
 *   - π  ≈ 3.14159265358979323846...
 *   - √2 ≈ 1.41421356237309504880...
 *   - e  ≈ 2.71828182845904523536...
 *   - exp(1) 应等于 e
 *   - ln(e) 应等于 1
 *   - √π ≈ 1.77245385090551602729...
 *   - Γ(0.5) = √π
 *   - 5! = 120
 *
 * 通过标准:
 *   - 常量级（π、√2、e）: 相对误差 ≤ 1e-45%（≈ 47 位有效数字）
 *     注: 参考值仅 50 位十进制数字，from_string 解析引入 ~1e-49% 舍入误差
 *   - 复合运算（exp、ln）: 相对误差 ≤ 1e-40%（≈ 42 位有效数字）
 *   - Γ 函数:             相对误差 ≤ 1e-30%（≈ 32 位有效数字）
 *   - 阶乘:               相对误差 ≤ 1e-60%（精确值）
 */
bool test_bigfloat_math() {
  std::cout << std::endl << "--- BigFloat Math Function Tests ---" << std::endl;
  int prec = 256; // 256 位精度
  bool all_pass = true;

  // 已知的高精度参考值（50+ 位有效数字）
  // 这些字符串来自公认的数学常量表

  // --- π ---
  std::cout << "  Computing pi..." << std::endl;
  BigFloat pi_val = BigFloat::pi(prec);
  std::string pi_str = pi_val.to_decimal_string(50);
  std::cout << "  pi       = " << pi_str << std::endl;
  std::cout
      << "  expected = 3.14159265358979323846264338327950288419716939937510"
      << std::endl;
  BigFloat pi_expected = BigFloat::from_string(
      "3.14159265358979323846264338327950288419716939937510", prec);
  all_pass &= check_relative_error(pi_val, pi_expected, 1e-45, "pi");

  // --- √2 ---
  std::cout << "  Computing sqrt(2)..." << std::endl;
  BigFloat two(2, prec);
  BigFloat sqrt2 = BigFloat::sqrt(two);
  std::cout << "  sqrt(2)  = " << sqrt2.to_decimal_string(50) << std::endl;
  std::cout
      << "  expected = 1.41421356237309504880168872420969807856967187537694"
      << std::endl;
  BigFloat sqrt2_expected = BigFloat::from_string(
      "1.41421356237309504880168872420969807856967187537694", prec);
  all_pass &= check_relative_error(sqrt2, sqrt2_expected, 1e-45, "sqrt(2)");

  // --- e ---
  std::cout << "  Computing e..." << std::endl;
  BigFloat e_val = BigFloat::e(prec);
  std::cout << "  e        = " << e_val.to_decimal_string(50) << std::endl;
  std::cout
      << "  expected = 2.71828182845904523536028747135266249775724709369995"
      << std::endl;
  BigFloat e_expected = BigFloat::from_string(
      "2.71828182845904523536028747135266249775724709369995", prec);
  all_pass &= check_relative_error(e_val, e_expected, 1e-45, "e");

  // --- exp(1) = e（交叉验证） ---
  BigFloat one(1, prec);
  BigFloat exp1 = BigFloat::exp(one);
  std::cout << "  exp(1)   = " << exp1.to_decimal_string(50) << std::endl;
  all_pass &= check_relative_error(exp1, e_expected, 1e-40, "exp(1) = e");

  // --- ln(e) = 1（对数验证） ---
  std::cout << "  Computing ln(e)..." << std::endl;
  BigFloat ln_e = BigFloat::ln(e_val);
  std::cout << "  ln(e)    = " << ln_e.to_decimal_string(50) << std::endl;
  all_pass &= check_relative_error(ln_e, one, 1e-40, "ln(e) = 1");

  // --- √π ---
  std::cout << "  Computing sqrt(pi)..." << std::endl;
  BigFloat sqrt_pi = BigFloat::sqrt(pi_val);
  std::cout << "  sqrt(pi) = " << sqrt_pi.to_decimal_string(50) << std::endl;
  std::cout
      << "  expected = 1.77245385090551602729816748334114518279754945612238"
      << std::endl;
  BigFloat sqrt_pi_expected = BigFloat::from_string(
      "1.77245385090551602729816748334114518279754945612238", prec);
  all_pass &=
      check_relative_error(sqrt_pi, sqrt_pi_expected, 1e-45, "sqrt(pi)");

  // --- Γ(0.5) = √π ---
  // 通过 half_factorial(0) 计算 (−0.5)! = Γ(0.5)
  std::cout << "  Computing Gamma(0.5) = (-0.5)! ..." << std::endl;
  BigFloat gamma_half = BigFloat::half_factorial(0, prec);
  std::cout << "  Gamma(0.5) = " << gamma_half.to_decimal_string(50)
            << std::endl;
  all_pass &= check_relative_error(gamma_half, sqrt_pi_expected, 1e-30,
                                   "Gamma(0.5) = sqrt(pi)");

  // --- 5! = 120 ---
  BigFloat f5 = BigFloat::factorial(5, prec);
  BigFloat expected_f5(120, prec);
  std::cout << "  5!       = " << f5.to_decimal_string(5) << std::endl;
  all_pass &= check_relative_error(f5, expected_f5, 1e-45, "5! = 120");

  assert(all_pass);
  return all_pass;
}

/**
 * @brief 测试入口
 */
int main() {
  std::cout << "=============================" << std::endl;
  std::cout << " BigFloat Unit Tests" << std::endl;
  std::cout << "=============================" << std::endl;

  bool all_pass = true;
  all_pass &= test_bigint_basic();
  all_pass &= test_bigfloat_arithmetic();
  all_pass &= test_bigfloat_math();

  std::cout << std::endl << "=============================" << std::endl;
  if (all_pass) {
    std::cout << " All tests PASSED." << std::endl;
  } else {
    std::cout << " Some tests FAILED!" << std::endl;
  }
  std::cout << "=============================" << std::endl;

  return all_pass ? 0 : 1;
}
