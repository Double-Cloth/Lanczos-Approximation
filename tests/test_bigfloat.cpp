#include "BigFloat.h"
#include "BigInt.h"
#include <iostream>
#include <cmath>
#include <string>
#include <cassert>

/**
 * BigFloat 基本功能测试
 */

bool test_bigint_basic() {
    std::cout << "--- BigInt Basic Tests ---" << std::endl;
    
    BigInt a(12345);
    BigInt b(67890);
    
    // 加法
    BigInt c = a + b;
    assert(c.to_decimal_string() == "80235");
    std::cout << "  12345 + 67890 = " << c << " [OK]" << std::endl;
    
    // 乘法
    BigInt d = a * b;
    assert(d.to_decimal_string() == "838102050");
    std::cout << "  12345 * 67890 = " << d << " [OK]" << std::endl;
    
    // 大数乘法
    BigInt e("999999999999999999");
    BigInt f("999999999999999999");
    BigInt g = e * f;
    std::cout << "  999999999999999999^2 = " << g << std::endl;
    assert(g.to_decimal_string() == "999999999999999998000000000000000001");
    std::cout << "  [OK]" << std::endl;
    
    // 除法
    BigInt h(100);
    BigInt i(7);
    auto [q, r] = h.divmod(i);
    std::cout << "  100 / 7 = " << q << " remainder " << r << std::endl;
    assert(q.to_decimal_string() == "14");
    assert(r.to_decimal_string() == "2");
    std::cout << "  [OK]" << std::endl;
    
    // 移位
    BigInt j(1);
    BigInt k = j << 100;
    std::cout << "  1 << 100 = " << k << std::endl;
    std::cout << "  bit_length = " << k.bit_length() << std::endl;
    assert(k.bit_length() == 101);
    std::cout << "  [OK]" << std::endl;
    
    return true;
}

bool test_bigfloat_arithmetic() {
    std::cout << std::endl << "--- BigFloat Arithmetic Tests ---" << std::endl;
    int prec = 128;
    
    // 加法
    BigFloat a(3, prec);
    BigFloat b(4, prec);
    BigFloat c = a + b;
    std::cout << "  3 + 4 = " << c.to_decimal_string(10) << std::endl;
    
    // 乘法
    BigFloat d = a * b;
    std::cout << "  3 * 4 = " << d.to_decimal_string(10) << std::endl;
    
    // 除法
    BigFloat e = a / b;
    std::cout << "  3 / 4 = " << e.to_decimal_string(20) << std::endl;
    
    // 1/3
    BigFloat one(1, prec);
    BigFloat three(3, prec);
    BigFloat third = one / three;
    std::cout << "  1 / 3 = " << third.to_decimal_string(40) << std::endl;
    
    return true;
}

bool test_bigfloat_math() {
    std::cout << std::endl << "--- BigFloat Math Function Tests ---" << std::endl;
    int prec = 256;

    // pi
    std::cout << "  Computing pi..." << std::endl;
    BigFloat pi_val = BigFloat::pi(prec);
    std::string pi_str = pi_val.to_decimal_string(50);
    std::cout << "  pi = " << pi_str << std::endl;
    // 验证前几位: 3.14159265358979323846...
    
    // sqrt(2)
    std::cout << "  Computing sqrt(2)..." << std::endl;
    BigFloat two(2, prec);
    BigFloat sqrt2 = BigFloat::sqrt(two);
    std::cout << "  sqrt(2) = " << sqrt2.to_decimal_string(50) << std::endl;
    // 1.41421356237309504880168872420969807856967187537694...
    
    // e
    std::cout << "  Computing e..." << std::endl;
    BigFloat e_val = BigFloat::e(prec);
    std::cout << "  e = " << e_val.to_decimal_string(50) << std::endl;
    // 2.71828182845904523536028747135266249775724709369995...
    
    // exp(1) = e
    BigFloat one(1, prec);
    BigFloat exp1 = BigFloat::exp(one);
    std::cout << "  exp(1) = " << exp1.to_decimal_string(50) << std::endl;
    
    // ln(e) = 1
    std::cout << "  Computing ln(e)..." << std::endl;
    BigFloat ln_e = BigFloat::ln(e_val);
    std::cout << "  ln(e) = " << ln_e.to_decimal_string(50) << std::endl;
    
    // sqrt(pi)
    std::cout << "  Computing sqrt(pi)..." << std::endl;
    BigFloat sqrt_pi = BigFloat::sqrt(pi_val);
    std::cout << "  sqrt(pi) = " << sqrt_pi.to_decimal_string(50) << std::endl;
    // 1.77245385090551602729816748334114518279754945612238...
    
    // Gamma(0.5) = sqrt(pi), 通过 half_factorial
    std::cout << "  Computing Gamma(0.5) = (-0.5)! ..." << std::endl;
    BigFloat gamma_half = BigFloat::half_factorial(0, prec);
    std::cout << "  Gamma(0.5) = " << gamma_half.to_decimal_string(50) << std::endl;

    // factorial
    BigFloat f5 = BigFloat::factorial(5, prec);
    std::cout << "  5! = " << f5.to_decimal_string(5) << std::endl;
    
    return true;
}

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
