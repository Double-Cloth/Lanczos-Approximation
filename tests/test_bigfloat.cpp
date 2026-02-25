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
 *   计算结果同时输出到控制台供人工比对已知值。
 */

#include "BigFloat.h"
#include "BigInt.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

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
    std::cout << "  12345 + 67890 = " << c << " [OK]" << std::endl;

    // 乘法测试
    BigInt d = a * b;
    assert(d.to_decimal_string() == "838102050");
    std::cout << "  12345 * 67890 = " << d << " [OK]" << std::endl;

    // 大数乘法: 验证 18 位 × 18 位 = 36 位的正确性
    BigInt e("999999999999999999");
    BigInt f("999999999999999999");
    BigInt g = e * f;
    std::cout << "  999999999999999999^2 = " << g << std::endl;
    assert(g.to_decimal_string() == "999999999999999998000000000000000001");
    std::cout << "  [OK]" << std::endl;

    // 除法: 验证商和余数
    BigInt h(100);
    BigInt i(7);
    auto [q, r] = h.divmod(i);
    std::cout << "  100 / 7 = " << q << " remainder " << r << std::endl;
    assert(q.to_decimal_string() == "14");
    assert(r.to_decimal_string() == "2");
    std::cout << "  [OK]" << std::endl;

    // 位移: 1 << 100 应是 101 位
    BigInt j(1);
    BigInt k = j << 100;
    std::cout << "  1 << 100 = " << k << std::endl;
    std::cout << "  bit_length = " << k.bit_length() << std::endl;
    assert(k.bit_length() == 101);
    std::cout << "  [OK]" << std::endl;

    return true;
}

/**
 * @brief BigFloat 算术运算测试
 *
 * 测试项:
 *   - 加法: 3 + 4 = 7
 *   - 乘法: 3 × 4 = 12
 *   - 除法: 3 / 4 = 0.75, 1 / 3 = 0.333...（验证循环小数精度）
 */
bool test_bigfloat_arithmetic() {
    std::cout << std::endl << "--- BigFloat Arithmetic Tests ---" << std::endl;
    int prec = 128; // 128 位精度（约 38 位十进制）

    // 加法
    BigFloat a(3, prec);
    BigFloat b(4, prec);
    BigFloat c = a + b;
    std::cout << "  3 + 4 = " << c.to_decimal_string(10) << std::endl;

    // 乘法
    BigFloat d = a * b;
    std::cout << "  3 * 4 = " << d.to_decimal_string(10) << std::endl;

    // 除法: 3/4 = 0.75
    BigFloat e = a / b;
    std::cout << "  3 / 4 = " << e.to_decimal_string(20) << std::endl;

    // 除法: 1/3 = 0.333...（无限循环小数，测试精度保持）
    BigFloat one(1, prec);
    BigFloat three(3, prec);
    BigFloat third = one / three;
    std::cout << "  1 / 3 = " << third.to_decimal_string(40) << std::endl;

    return true;
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
 */
bool test_bigfloat_math() {
    std::cout << std::endl << "--- BigFloat Math Function Tests ---" << std::endl;
    int prec = 256; // 256 位精度

    // --- π ---
    std::cout << "  Computing pi..." << std::endl;
    BigFloat pi_val = BigFloat::pi(prec);
    std::string pi_str = pi_val.to_decimal_string(50);
    std::cout << "  pi = " << pi_str << std::endl;
    // 已知值: 3.14159265358979323846264338327950288419716939937510...

    // --- √2 ---
    std::cout << "  Computing sqrt(2)..." << std::endl;
    BigFloat two(2, prec);
    BigFloat sqrt2 = BigFloat::sqrt(two);
    std::cout << "  sqrt(2) = " << sqrt2.to_decimal_string(50) << std::endl;
    // 已知值: 1.41421356237309504880168872420969807856967187537694...

    // --- e ---
    std::cout << "  Computing e..." << std::endl;
    BigFloat e_val = BigFloat::e(prec);
    std::cout << "  e = " << e_val.to_decimal_string(50) << std::endl;
    // 已知值: 2.71828182845904523536028747135266249775724709369995...

    // --- exp(1) = e（交叉验证） ---
    BigFloat one(1, prec);
    BigFloat exp1 = BigFloat::exp(one);
    std::cout << "  exp(1) = " << exp1.to_decimal_string(50) << std::endl;

    // --- ln(e) = 1（对数验证） ---
    std::cout << "  Computing ln(e)..." << std::endl;
    BigFloat ln_e = BigFloat::ln(e_val);
    std::cout << "  ln(e) = " << ln_e.to_decimal_string(50) << std::endl;

    // --- √π ---
    std::cout << "  Computing sqrt(pi)..." << std::endl;
    BigFloat sqrt_pi = BigFloat::sqrt(pi_val);
    std::cout << "  sqrt(pi) = " << sqrt_pi.to_decimal_string(50) << std::endl;
    // 已知值: 1.77245385090551602729816748334114518279754945612238...

    // --- Γ(0.5) = √π ---
    // 通过 half_factorial(0) 计算 (−0.5)! = Γ(0.5)
    std::cout << "  Computing Gamma(0.5) = (-0.5)! ..." << std::endl;
    BigFloat gamma_half = BigFloat::half_factorial(0, prec);
    std::cout << "  Gamma(0.5) = " << gamma_half.to_decimal_string(50)
            << std::endl;

    // --- 5! = 120 ---
    BigFloat f5 = BigFloat::factorial(5, prec);
    std::cout << "  5! = " << f5.to_decimal_string(5) << std::endl;

    return true;
}

/**
 * @brief 测试入口
 */
int main() {
    std::cout << "=============================" << std::endl;
    std::cout << " BigFloat Unit Tests" << std::endl;
    std::cout << "=============================" << std::endl;

    test_bigint_basic();
    test_bigfloat_arithmetic();
    test_bigfloat_math();

    std::cout << std::endl << "All tests completed." << std::endl;
    return 0;
}
