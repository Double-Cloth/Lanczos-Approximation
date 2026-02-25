#include "BigInt.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>

// ============================================
// 构造函数
// ============================================

BigInt::BigInt() : digits_(1, 0) {}

BigInt::BigInt(uint64_t val) {
    if (val == 0) {
        digits_.push_back(0);
    } else {
        digits_.push_back(static_cast<uint32_t>(val & 0xFFFFFFFFULL));
        uint32_t hi = static_cast<uint32_t>(val >> 32);
        if (hi != 0) {
            digits_.push_back(hi);
        }
    }
}

BigInt::BigInt(const std::string& decimal_str) {
    if (decimal_str.empty()) {
        digits_.push_back(0);
        return;
    }
    // 逐字符处理：result = result * 10 + digit
    digits_.push_back(0);
    for (char c : decimal_str) {
        if (c < '0' || c > '9') continue;
        // *this = *this * 10 + (c - '0')
        *this = this->mul_u32(10);
        uint64_t carry = static_cast<uint64_t>(c - '0');
        for (size_t i = 0; i < digits_.size() && carry > 0; i++) {
            carry += digits_[i];
            digits_[i] = static_cast<uint32_t>(carry & 0xFFFFFFFFULL);
            carry >>= 32;
        }
        if (carry > 0) {
            digits_.push_back(static_cast<uint32_t>(carry));
        }
    }
    trim();
}

// ============================================
// 辅助
// ============================================

void BigInt::trim() {
    while (digits_.size() > 1 && digits_.back() == 0) {
        digits_.pop_back();
    }
}

bool BigInt::is_zero() const {
    return digits_.size() == 1 && digits_[0] == 0;
}

int BigInt::bit_length() const {
    if (is_zero()) return 0;
    int top = static_cast<int>(digits_.size()) - 1;
    uint32_t hi = digits_[top];
    int bits = top * 32;
    while (hi > 0) {
        bits++;
        hi >>= 1;
    }
    return bits;
}

bool BigInt::get_bit(int pos) const {
    int word = pos / 32;
    int bit = pos % 32;
    if (word >= static_cast<int>(digits_.size())) return false;
    return (digits_[word] >> bit) & 1;
}

// ============================================
// 比较
// ============================================

int BigInt::compare(const BigInt& other) const {
    if (digits_.size() != other.digits_.size()) {
        return digits_.size() < other.digits_.size() ? -1 : 1;
    }
    for (int i = static_cast<int>(digits_.size()) - 1; i >= 0; i--) {
        if (digits_[i] != other.digits_[i]) {
            return digits_[i] < other.digits_[i] ? -1 : 1;
        }
    }
    return 0;
}

// ============================================
// 加法
// ============================================

BigInt BigInt::operator+(const BigInt& other) const {
    BigInt result;
    size_t n = std::max(digits_.size(), other.digits_.size());
    result.digits_.resize(n, 0);
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t a = (i < digits_.size()) ? digits_[i] : 0;
        uint64_t b = (i < other.digits_.size()) ? other.digits_[i] : 0;
        carry += a + b;
        result.digits_[i] = static_cast<uint32_t>(carry & 0xFFFFFFFFULL);
        carry >>= 32;
    }
    if (carry > 0) {
        result.digits_.push_back(static_cast<uint32_t>(carry));
    }
    return result;
}

BigInt& BigInt::operator+=(const BigInt& other) {
    *this = *this + other;
    return *this;
}

// ============================================
// 减法 (要求 *this >= other)
// ============================================

BigInt BigInt::operator-(const BigInt& other) const {
    assert(*this >= other);
    BigInt result;
    result.digits_.resize(digits_.size(), 0);
    int64_t borrow = 0;
    for (size_t i = 0; i < digits_.size(); i++) {
        int64_t a = digits_[i];
        int64_t b = (i < other.digits_.size()) ? other.digits_[i] : 0;
        int64_t diff = a - b - borrow;
        if (diff < 0) {
            diff += (1LL << 32);
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.digits_[i] = static_cast<uint32_t>(diff);
    }
    result.trim();
    return result;
}

BigInt& BigInt::operator-=(const BigInt& other) {
    *this = *this - other;
    return *this;
}

// ============================================
// 乘法: 经典 O(n*m) 算法
// ============================================

BigInt BigInt::operator*(const BigInt& other) const {
    if (is_zero() || other.is_zero()) return BigInt(0);
    
    size_t n = digits_.size();
    size_t m = other.digits_.size();
    BigInt result;
    result.digits_.resize(n + m, 0);
    
    for (size_t i = 0; i < n; i++) {
        uint64_t carry = 0;
        for (size_t j = 0; j < m; j++) {
            uint64_t prod = static_cast<uint64_t>(digits_[i]) * other.digits_[j]
                          + result.digits_[i + j] + carry;
            result.digits_[i + j] = static_cast<uint32_t>(prod & 0xFFFFFFFFULL);
            carry = prod >> 32;
        }
        if (carry > 0) {
            result.digits_[i + m] += static_cast<uint32_t>(carry);
        }
    }
    result.trim();
    return result;
}

BigInt& BigInt::operator*=(const BigInt& other) {
    *this = *this * other;
    return *this;
}

// 乘以单个 uint32_t
BigInt BigInt::mul_u32(uint32_t val) const {
    if (val == 0 || is_zero()) return BigInt(0);
    BigInt result;
    result.digits_.resize(digits_.size(), 0);
    uint64_t carry = 0;
    for (size_t i = 0; i < digits_.size(); i++) {
        uint64_t prod = static_cast<uint64_t>(digits_[i]) * val + carry;
        result.digits_[i] = static_cast<uint32_t>(prod & 0xFFFFFFFFULL);
        carry = prod >> 32;
    }
    if (carry > 0) {
        result.digits_.push_back(static_cast<uint32_t>(carry));
    }
    return result;
}

// ============================================
// 除法: 经典长除法
// ============================================

// 除以 uint32_t, 返回 {商, 余数}
std::pair<BigInt, uint32_t> BigInt::divmod_u32(uint32_t val) const {
    if (val == 0) throw std::runtime_error("Division by zero");
    BigInt quotient;
    quotient.digits_.resize(digits_.size(), 0);
    uint64_t rem = 0;
    for (int i = static_cast<int>(digits_.size()) - 1; i >= 0; i--) {
        rem = (rem << 32) | digits_[i];
        quotient.digits_[i] = static_cast<uint32_t>(rem / val);
        rem %= val;
    }
    quotient.trim();
    return {quotient, static_cast<uint32_t>(rem)};
}

// 大数除大数: Knuth Algorithm D 的简化实现
std::pair<BigInt, BigInt> BigInt::divmod(const BigInt& divisor) const {
    if (divisor.is_zero()) throw std::runtime_error("Division by zero");
    
    int cmp = compare(divisor);
    if (cmp < 0) return {BigInt(0), *this};
    if (cmp == 0) return {BigInt(1), BigInt(0)};
    
    // 对于单字除数，使用快速路径
    if (divisor.digits_.size() == 1) {
        auto [q, r] = divmod_u32(divisor.digits_[0]);
        return {q, BigInt(r)};
    }
    
    // 二进制长除法 (简单但可靠)
    BigInt quotient(0);
    BigInt remainder(0);
    
    int total_bits = bit_length();
    for (int i = total_bits - 1; i >= 0; i--) {
        remainder <<= 1;
        if (get_bit(i)) {
            remainder.digits_[0] |= 1;
        }
        if (remainder >= divisor) {
            remainder -= divisor;
            // 设置 quotient 的第 i 位
            int word = i / 32;
            int bit = i % 32;
            while (static_cast<int>(quotient.digits_.size()) <= word) {
                quotient.digits_.push_back(0);
            }
            quotient.digits_[word] |= (1U << bit);
        }
    }
    quotient.trim();
    remainder.trim();
    return {quotient, remainder};
}

BigInt BigInt::operator/(const BigInt& other) const {
    return divmod(other).first;
}

BigInt BigInt::operator%(const BigInt& other) const {
    return divmod(other).second;
}

// ============================================
// 位移
// ============================================

BigInt BigInt::operator<<(int bits) const {
    if (bits == 0 || is_zero()) return *this;
    
    int word_shift = bits / 32;
    int bit_shift = bits % 32;
    
    BigInt result;
    result.digits_.resize(digits_.size() + word_shift + 1, 0);
    
    uint64_t carry = 0;
    for (size_t i = 0; i < digits_.size(); i++) {
        uint64_t val = (static_cast<uint64_t>(digits_[i]) << bit_shift) | carry;
        result.digits_[i + word_shift] = static_cast<uint32_t>(val & 0xFFFFFFFFULL);
        carry = val >> 32;
    }
    if (carry > 0) {
        result.digits_[digits_.size() + word_shift] = static_cast<uint32_t>(carry);
    }
    result.trim();
    return result;
}

BigInt& BigInt::operator<<=(int bits) {
    *this = *this << bits;
    return *this;
}

BigInt BigInt::operator>>(int bits) const {
    if (bits == 0) return *this;
    
    int word_shift = bits / 32;
    int bit_shift = bits % 32;
    
    if (word_shift >= static_cast<int>(digits_.size())) return BigInt(0);
    
    BigInt result;
    size_t new_size = digits_.size() - word_shift;
    result.digits_.resize(new_size, 0);
    
    for (size_t i = 0; i < new_size; i++) {
        result.digits_[i] = digits_[i + word_shift] >> bit_shift;
        if (bit_shift > 0 && i + word_shift + 1 < digits_.size()) {
            result.digits_[i] |= digits_[i + word_shift + 1] << (32 - bit_shift);
        }
    }
    result.trim();
    return result;
}

BigInt& BigInt::operator>>=(int bits) {
    *this = *this >> bits;
    return *this;
}

// ============================================
// 十进制输出
// ============================================

std::string BigInt::to_decimal_string() const {
    if (is_zero()) return "0";
    
    std::string result;
    BigInt temp = *this;
    while (!temp.is_zero()) {
        auto [q, r] = temp.divmod_u32(10);
        result += static_cast<char>('0' + r);
        temp = q;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

std::ostream& operator<<(std::ostream& os, const BigInt& b) {
    os << b.to_decimal_string();
    return os;
}
