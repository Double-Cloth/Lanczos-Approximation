#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <iostream>

/**
 * 无符号大整数类
 * 内部使用 vector<uint32_t> 存储，低位在前 (little-endian)
 * 基数 = 2^32
 */
class BigInt {
public:
    // === 构造 ===
    BigInt();
    BigInt(uint64_t val);
    BigInt(const std::string& decimal_str);
    BigInt(const BigInt& other) = default;
    BigInt(BigInt&& other) noexcept = default;
    BigInt& operator=(const BigInt& other) = default;
    BigInt& operator=(BigInt&& other) noexcept = default;

    // === 比较 ===
    int compare(const BigInt& other) const;
    bool operator==(const BigInt& other) const { return compare(other) == 0; }
    bool operator!=(const BigInt& other) const { return compare(other) != 0; }
    bool operator<(const BigInt& other) const  { return compare(other) < 0; }
    bool operator<=(const BigInt& other) const { return compare(other) <= 0; }
    bool operator>(const BigInt& other) const  { return compare(other) > 0; }
    bool operator>=(const BigInt& other) const { return compare(other) >= 0; }

    // === 算术（返回新对象）===
    BigInt operator+(const BigInt& other) const;
    BigInt operator-(const BigInt& other) const; // 要求 *this >= other
    BigInt operator*(const BigInt& other) const;
    
    // 除法: 返回 {商, 余数}
    std::pair<BigInt, BigInt> divmod(const BigInt& divisor) const;
    BigInt operator/(const BigInt& other) const;
    BigInt operator%(const BigInt& other) const;

    // 与小数相乘/相除（便捷方法）
    BigInt mul_u32(uint32_t val) const;
    std::pair<BigInt, uint32_t> divmod_u32(uint32_t val) const;

    // === 复合赋值 ===
    BigInt& operator+=(const BigInt& other);
    BigInt& operator-=(const BigInt& other);
    BigInt& operator*=(const BigInt& other);

    // === 位移 ===
    BigInt operator<<(int bits) const;    // 左移
    BigInt operator>>(int bits) const;    // 右移
    BigInt& operator<<=(int bits);
    BigInt& operator>>=(int bits);

    // === 属性 ===
    bool is_zero() const;
    int bit_length() const;     // 最高位的位置 (0 表示值为 0)
    bool get_bit(int pos) const;

    // === I/O ===
    std::string to_decimal_string() const;
    friend std::ostream& operator<<(std::ostream& os, const BigInt& b);

    // === 内部访问 ===
    const std::vector<uint32_t>& digits() const { return digits_; }
    std::vector<uint32_t>& digits_mut() { return digits_; }

    // 去除高位的零
    void trim();

private:
    std::vector<uint32_t> digits_;  // 低位在前, digits_[0] 是最低 32-bit
};
