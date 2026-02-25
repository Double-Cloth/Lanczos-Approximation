#pragma once
#include "BigInt.h"
#include <string>
#include <cstdint>

/**
 * 任意精度浮点数类
 * 
 * 表示: value = (-1)^sign * mantissa * 2^exponent
 * 
 * mantissa 是一个 BigInt (无符号)
 * mantissa 的位长被限制在 precision_bits 左右（运算后截断）
 * 
 * 约定: 零值时 mantissa=0, exponent=0, sign=false
 */
class BigFloat {
public:
    // === 构造 ===
    BigFloat();
    BigFloat(int val, int prec_bits = 256);
    BigFloat(int64_t val, int prec_bits = 256);
    BigFloat(double val, int prec_bits = 256);

    // 从分数 num/den 构造
    static BigFloat from_fraction(const BigInt& num, const BigInt& den, int prec_bits);
    
    // 从 BigInt 构造
    static BigFloat from_bigint(const BigInt& val, int prec_bits);
    
    // 从十进制字符串构造
    static BigFloat from_string(const std::string& str, int prec_bits);

    // 复制/移动
    BigFloat(const BigFloat&) = default;
    BigFloat(BigFloat&&) noexcept = default;
    BigFloat& operator=(const BigFloat&) = default;
    BigFloat& operator=(BigFloat&&) noexcept = default;

    // === 精度 ===
    int precision() const { return precision_bits_; }
    void set_precision(int bits) { precision_bits_ = bits; normalize(); }

    // === 符号 ===
    bool is_negative() const { return sign_; }
    bool is_zero() const { return mantissa_.is_zero(); }
    BigFloat operator-() const;
    BigFloat abs() const;

    // === 比较 ===
    int compare(const BigFloat& other) const;
    bool operator==(const BigFloat& other) const { return compare(other) == 0; }
    bool operator!=(const BigFloat& other) const { return compare(other) != 0; }
    bool operator<(const BigFloat& other) const  { return compare(other) < 0; }
    bool operator<=(const BigFloat& other) const { return compare(other) <= 0; }
    bool operator>(const BigFloat& other) const  { return compare(other) > 0; }
    bool operator>=(const BigFloat& other) const { return compare(other) >= 0; }

    // === 算术 ===
    BigFloat operator+(const BigFloat& other) const;
    BigFloat operator-(const BigFloat& other) const;
    BigFloat operator*(const BigFloat& other) const;
    BigFloat operator/(const BigFloat& other) const;
    BigFloat& operator+=(const BigFloat& other);
    BigFloat& operator-=(const BigFloat& other);
    BigFloat& operator*=(const BigFloat& other);
    BigFloat& operator/=(const BigFloat& other);

    // 乘以 2^n
    BigFloat mul_pow2(int n) const;

    // === 数学函数（静态）===
    static BigFloat pi(int prec_bits);
    static BigFloat e(int prec_bits);
    static BigFloat sqrt(const BigFloat& x);
    static BigFloat exp(const BigFloat& x);
    static BigFloat ln(const BigFloat& x);
    static BigFloat pow(const BigFloat& base, const BigFloat& exp);
    
    // 阶乘和半整数阶乘:  (k - 0.5)! = Gamma(k + 0.5)
    static BigFloat factorial(int n, int prec_bits);
    static BigFloat half_factorial(int k, int prec_bits); // (k-0.5)! 

    // === I/O ===
    // 转换为十进制字符串, decimal_digits 个小数有效数字
    // 默认输出科学记数法
    std::string to_decimal_string(int decimal_digits, bool scientific_notation = true) const;
    friend std::ostream& operator<<(std::ostream& os, const BigFloat& f);

    // === 内部访问 ===
    const BigInt& mantissa() const { return mantissa_; }
    int64_t exponent() const { return exponent_; }

private:
    bool sign_;            // true = negative
    BigInt mantissa_;      // 无符号尾数
    int64_t exponent_;     // 二进制指数
    int precision_bits_;   // 目标精度

    // 截断尾数到 precision_bits_ 位
    void normalize();
    
    // 对齐两个 BigFloat 的指数，使尾数可以直接加减
    static void align(const BigFloat& a, const BigFloat& b,
                      BigInt& a_mant, BigInt& b_mant, int64_t& common_exp);
};
