/**
 * @file BigInt.cpp
 * @brief BigInt 类的实现 —— 无符号任意精度大整数
 *
 * @details
 * 内部使用 vector<uint32_t> 存储，以 2^32 为基数，小端序排列。
 * 所有运算均为无符号运算，减法要求被减数 ≥ 减数。
 */

#include "BigInt.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>

// ============================================
// 构造函数
// ============================================

/**
 * @brief 默认构造: 值为 0（一个元素，值为 0）
 */
BigInt::BigInt() : digits_(1, 0) {
}

/**
 * @brief 从 uint64_t 构造
 * @details 将 64 位整数拆分为低 32 位和高 32 位存入 digits_
 */
BigInt::BigInt(const uint64_t val) {
    if (val == 0) {
        digits_.push_back(0);
    } else {
        // 低 32 位
        digits_.push_back(static_cast<uint32_t>(val & 0xFFFFFFFFULL));
        // 高 32 位（仅在非零时添加）
        if (const auto hi = static_cast<uint32_t>(val >> 32); hi != 0) {
            digits_.push_back(hi);
        }
    }
}

/**
 * @brief 从十进制字符串构造
 * @details 逐字符处理: result = result × 10 + digit
 *          非数字字符（如空格、逗号）被静默跳过
 */
BigInt::BigInt(const std::string &decimal_str) {
    if (decimal_str.empty()) {
        digits_.push_back(0);
        return;
    }
    digits_.push_back(0);
    for (const char c: decimal_str) {
        if (c < '0' || c > '9')
            continue; // 跳过非数字字符
        // *this = *this * 10 + (c - '0')
        *this = this->mul_u32(10);
        auto carry = static_cast<uint64_t>(c - '0');
        for (size_t i = 0; i < digits_.size() && carry > 0; i++) {
            carry += digits_[i];
            digits_[i] = static_cast<uint32_t>(carry & 0xFFFFFFFFULL);
            carry >>= 32;
        }
        if (carry > 0) {
            digits_.push_back(static_cast<uint32_t>(carry));
        }
    }
    trim(); // 去除可能的前导零
}

// ============================================
// 辅助函数
// ============================================

/**
 * @brief 去除高位多余的零，保证至少保留一个元素
 * @details 维护不变量: digits_.back() != 0（除非数值为 0）
 */
void BigInt::trim() {
    while (digits_.size() > 1 && digits_.back() == 0) {
        digits_.pop_back();
    }
}

/**
 * @brief 判断是否为 0
 * @details 由于 trim() 的不变量，值为 0 等价于只有一个元素且为 0
 */
bool BigInt::is_zero() const { return digits_.size() == 1 && digits_[0] == 0; }

/**
 * @brief 计算二进制位长（最高有效位的位置）
 * @return 0 表示值为 0; 否则返回 floor(log2(val)) + 1
 * @details 先计算最高字的位长，再加上低位字的总位数
 */
int BigInt::bit_length() const {
    if (is_zero())
        return 0;
    // 最高字的索引
    const int top = static_cast<int>(digits_.size()) - 1;
    uint32_t hi = digits_[top];
    // 低位字贡献 top * 32 位
    int bits = top * 32;
    // 计算最高字中有效位数
    while (hi > 0) {
        bits++;
        hi >>= 1;
    }
    return bits;
}

/**
 * @brief 获取第 pos 位的值（从 0 开始计数，0 = 最低位）
 * @return true 表示该位为 1
 */
bool BigInt::get_bit(const int pos) const {
    const int word = pos / 32; // 所在的 uint32_t 索引
    const int bit = pos % 32; // 在该字内的位偏移
    if (word >= static_cast<int>(digits_.size()))
        return false;
    return (digits_[word] >> bit) & 1;
}

// ============================================
// 比较
// ============================================

/**
 * @brief 三路比较: 先比较字长，相同则从最高字逐字比较
 * @return -1 / 0 / 1
 */
int BigInt::compare(const BigInt &other) const {
    // 字长不同，字长大者值更大
    if (digits_.size() != other.digits_.size()) {
        return digits_.size() < other.digits_.size() ? -1 : 1;
    }
    // 字长相同，从最高字向最低字逐字比较
    for (int i = static_cast<int>(digits_.size()) - 1; i >= 0; i--) {
        if (digits_[i] != other.digits_[i]) {
            return digits_[i] < other.digits_[i] ? -1 : 1;
        }
    }
    return 0; // 完全相同
}

// ============================================
// 加法: 逐字相加并传播进位
// ============================================

/**
 * @brief 无符号加法
 * @details 从低位到高位逐字相加，carry 在每步传播
 *          时间复杂度: O(max(n, m))
 */
BigInt BigInt::operator+(const BigInt &other) const {
    BigInt result;
    const size_t n = std::max(digits_.size(), other.digits_.size());
    result.digits_.resize(n, 0);
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        const uint64_t a = i < digits_.size() ? digits_[i] : 0;
        const uint64_t b = i < other.digits_.size() ? other.digits_[i] : 0;
        carry += a + b;
        result.digits_[i] = static_cast<uint32_t>(carry & 0xFFFFFFFFULL);
        carry >>= 32;
    }
    if (carry > 0) {
        result.digits_.push_back(static_cast<uint32_t>(carry));
    }
    return result;
}

BigInt &BigInt::operator+=(const BigInt &other) {
    *this = *this + other;
    return *this;
}

// ============================================
// 减法: 逐字相减并传播借位（要求 *this >= other）
// ============================================

/**
 * @brief 无符号减法
 * @pre *this >= other
 * @details 从低位到高位逐字相减，借位通过 borrow 传播
 *          时间复杂度: O(n)
 */
BigInt BigInt::operator-(const BigInt &other) const {
    assert(*this >= other); // 前提条件: 被减数 ≥ 减数
    BigInt result;
    result.digits_.resize(digits_.size(), 0);
    int64_t borrow = 0;
    for (size_t i = 0; i < digits_.size(); i++) {
        const int64_t a = digits_[i];
        const int64_t b = i < other.digits_.size() ? other.digits_[i] : 0;
        int64_t diff = a - b - borrow;
        if (diff < 0) {
            diff += 1LL << 32; // 向高位借 2^32
            borrow = 1;
        } else {
            borrow = 0;
        }
        result.digits_[i] = static_cast<uint32_t>(diff);
    }
    result.trim();
    return result;
}

BigInt &BigInt::operator-=(const BigInt &other) {
    *this = *this - other;
    return *this;
}

// ============================================
// 乘法: 经典 O(n×m) 长乘法
// ============================================

/**
 * @brief 经典长乘法算法
 * @details 对于 n 字 × m 字的乘法:
 *          result[i+j] += digits_[i] × other.digits_[j]
 *          进位逐步向高位传播。
 *          结果最多 n+m 字。时间复杂度 O(n×m)。
 *
 * @note 对于超大数，可考虑 Karatsuba 或 FFT 乘法优化，
 *       但当前场景下经典算法已足够。
 */
/**
 * @brief 乘法: 自动选择经典算法或 Karatsuba 分治
 * @details 当两操作数均 > KARATSUBA_THRESHOLD 字时启用 Karatsuba，
 *          否则使用经典 O(n×m) 长乘法。
 */
BigInt BigInt::operator*(const BigInt &other) const {
    if (is_zero() || other.is_zero())
        return BigInt(0);

    const size_t n = digits_.size();
    const size_t m = other.digits_.size();

    // 当两操作数都超过阈值时启用 Karatsuba
    if (n >= KARATSUBA_THRESHOLD && m >= KARATSUBA_THRESHOLD) {
        return karatsuba_mul(*this, other);
    }

    // 经典长乘法
    BigInt result;
    result.digits_.resize(n + m, 0);

    for (size_t i = 0; i < n; i++) {
        uint64_t carry = 0;
        for (size_t j = 0; j < m; j++) {
            uint64_t prod = static_cast<uint64_t>(digits_[i]) * other.digits_[j] +
                            result.digits_[i + j] + carry;
            result.digits_[i + j] = static_cast<uint32_t>(prod & 0xFFFFFFFFULL);
            carry = prod >> 32;
        }
        for (size_t k = i + m; carry > 0 && k < result.digits_.size(); k++) {
            uint64_t sum = static_cast<uint64_t>(result.digits_[k]) + carry;
            result.digits_[k] = static_cast<uint32_t>(sum & 0xFFFFFFFFULL);
            carry = sum >> 32;
        }
    }
    result.trim();
    return result;
}

/**
 * @brief Karatsuba 分治乘法
 * @details 将两个大数各分为高低两半，使用 3 次递归乘法代替 4 次：
 *   a = a_hi × B + a_lo,  b = b_hi × B + b_lo  (B = 2^(half*32))
 *   z0 = a_lo × b_lo
 *   z2 = a_hi × b_hi
 *   z1 = (a_lo + a_hi) × (b_lo + b_hi) - z0 - z2
 *   result = z2 × B² + z1 × B + z0
 */
BigInt BigInt::karatsuba_mul(const BigInt &a, const BigInt &b) {
    size_t n = a.digits_.size();
    size_t m = b.digits_.size();

    // 小规模回退到经典乘法
    if (n < KARATSUBA_THRESHOLD || m < KARATSUBA_THRESHOLD) {
        // 经典乘法（内联）
        BigInt result;
        result.digits_.resize(n + m, 0);
        for (size_t i = 0; i < n; i++) {
            uint64_t carry = 0;
            for (size_t j = 0; j < m; j++) {
                uint64_t prod = static_cast<uint64_t>(a.digits_[i]) * b.digits_[j] +
                                result.digits_[i + j] + carry;
                result.digits_[i + j] = static_cast<uint32_t>(prod & 0xFFFFFFFFULL);
                carry = prod >> 32;
            }
            for (size_t k = i + m; carry > 0 && k < result.digits_.size(); k++) {
                uint64_t sum = static_cast<uint64_t>(result.digits_[k]) + carry;
                result.digits_[k] = static_cast<uint32_t>(sum & 0xFFFFFFFFULL);
                carry = sum >> 32;
            }
        }
        result.trim();
        return result;
    }

    // 分割点: 取两操作数较大者的一半
    size_t half = std::max(n, m) / 2;

    // a = a_hi * B^half + a_lo
    BigInt a_lo, a_hi;
    if (half < n) {
        a_lo.digits_.assign(a.digits_.begin(), a.digits_.begin() + half);
        a_hi.digits_.assign(a.digits_.begin() + half, a.digits_.end());
    } else {
        a_lo = a;
        a_hi = BigInt(0);
    }
    a_lo.trim();
    a_hi.trim();

    // b = b_hi * B^half + b_lo
    BigInt b_lo, b_hi;
    if (half < m) {
        b_lo.digits_.assign(b.digits_.begin(), b.digits_.begin() + half);
        b_hi.digits_.assign(b.digits_.begin() + half, b.digits_.end());
    } else {
        b_lo = b;
        b_hi = BigInt(0);
    }
    b_lo.trim();
    b_hi.trim();

    // 3 次递归乘法
    BigInt z0 = a_lo * b_lo; // a_lo * b_lo
    BigInt z2 = a_hi * b_hi; // a_hi * b_hi
    BigInt z1 = (a_lo + a_hi) * (b_lo + b_hi); // (a_lo+a_hi) * (b_lo+b_hi)
    z1 = z1 - z0 - z2; // z1 = cross term

    // result = z2 * B^(2*half) + z1 * B^half + z0
    int shift_bits = static_cast<int>(half) * 32;
    BigInt result = (z2 << (2 * shift_bits)) + (z1 << shift_bits) + z0;
    result.trim();
    return result;
}

BigInt &BigInt::operator*=(const BigInt &other) {
    *this = *this * other;
    return *this;
}

/**
 * @brief 乘以单个 uint32_t：O(n) 快速路径
 * @details 避免了通用乘法的双重循环开销
 */
BigInt BigInt::mul_u32(uint32_t val) const {
    if (val == 0 || is_zero())
        return BigInt(0);
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
// 除法
// ============================================

/**
 * @brief 除以单个 uint32_t，返回 {商, 余数}
 * @throws std::runtime_error 若 val 为 0
 * @details 从最高字到最低字逐字处理:
 *          rem = (rem << 32) | digit[i]
 *          quotient[i] = rem / val
 *          rem = rem % val
 *          时间复杂度: O(n)
 */
std::pair<BigInt, uint32_t> BigInt::divmod_u32(uint32_t val) const {
    if (val == 0)
        throw std::runtime_error("Division by zero");
    BigInt quotient;
    quotient.digits_.resize(digits_.size(), 0);
    uint64_t rem = 0;
    // 从最高字向最低字处理
    for (int i = static_cast<int>(digits_.size()) - 1; i >= 0; i--) {
        rem = (rem << 32) | digits_[i];
        quotient.digits_[i] = static_cast<uint32_t>(rem / val);
        rem %= val;
    }
    quotient.trim();
    return {quotient, static_cast<uint32_t>(rem)};
}

/**
 * @brief 大数除大数: 二进制长除法
 * @throws std::runtime_error 若 divisor 为 0
 * @details 算法流程:
 *          1. 快速路径: 单字除数使用 divmod_u32
 *          2. 通用路径: 从被除数的最高位到最低位逐位处理：
 *             - remainder 左移 1 位
 *             - 设置最低位为被除数的当前位
 *             - 若 remainder ≥ divisor，则 remainder -= divisor，
 *               并在商中设置对应位
 *          时间复杂度: O(n × bit_length)，其中 n 为除数的字长
 *
 * @note 这是简化版的二进制长除法，效率不如 Knuth Algorithm D，
 *       但实现更简单可靠。对于本项目的使用场景已足够。
 */
/**
 * @brief 大数除大数: Knuth Algorithm D
 * @throws std::runtime_error 若 divisor 为 0
 * @details 算法流程:
 *   1. 快速路径: 单字除数使用 divmod_u32
 *   2. Knuth Algorithm D:
 *      a. 规范化: 左移使除数最高字 >= BASE/2
 *      b. 对每个商字: 估算 q̂ = (u[j+n]*B + u[j+n-1]) / v[n-1]
 *      c. 修正: 最多修正 2 次
 *      d. u -= q̂ * v, 若借位则 u += v 并 q̂--
 *   复杂度: O(n×m)，远快于二进制逐位法的 O(n × bitlen)
 */
std::pair<BigInt, BigInt> BigInt::divmod(const BigInt &divisor) const {
    if (divisor.is_zero())
        throw std::runtime_error("Division by zero");

    int cmp = compare(divisor);
    if (cmp < 0)
        return {BigInt(0), *this};
    if (cmp == 0)
        return {BigInt(1), BigInt(0)};

    // 单字除数: 使用 O(n) 快速路径
    if (divisor.digits_.size() == 1) {
        auto [q, r] = divmod_u32(divisor.digits_[0]);
        return {q, BigInt(r)};
    }

    // === Knuth Algorithm D ===
    size_t n = divisor.digits_.size();
    size_t m = digits_.size() - n;

    // D1: 规范化 — 左移使除数最高字的最高位为 1
    uint32_t d_top = divisor.digits_[n - 1];
    int s = 0; // 规范化移位量
    while ((d_top & 0x80000000U) == 0) {
        d_top <<= 1;
        s++;
    }

    // 对除数和被除数都左移 s 位
    BigInt v = divisor << s;
    BigInt u = *this << s;

    // 确保 u 有 n+m+1 个字（最高位可能为 0）
    while (u.digits_.size() <= n + m) {
        u.digits_.push_back(0);
    }

    BigInt quotient;
    quotient.digits_.resize(m + 1, 0);

    // D2-D7: 主循环，对每个商字 j
    for (int j = static_cast<int>(m); j >= 0; j--) {
        // D3: 估算商字 q̂
        uint64_t u_hi =
                (static_cast<uint64_t>(u.digits_[j + n]) << 32) | u.digits_[j + n - 1];
        uint64_t qhat = u_hi / v.digits_[n - 1];
        uint64_t rhat = u_hi % v.digits_[n - 1];

        // 修正 qhat（Knuth 的精确调整）
        while (qhat >= 0x100000000ULL ||
               (n >= 2 &&
                qhat * v.digits_[n - 2] > (rhat << 32) + u.digits_[j + n - 2])) {
            qhat--;
            rhat += v.digits_[n - 1];
            if (rhat >= 0x100000000ULL)
                break;
        }

        // D4: u[j..j+n] -= qhat * v
        int64_t borrow = 0;
        for (size_t i = 0; i < n; i++) {
            uint64_t prod = qhat * v.digits_[i];
            int64_t diff = static_cast<int64_t>(u.digits_[j + i]) -
                           static_cast<int64_t>(prod & 0xFFFFFFFFULL) - borrow;
            u.digits_[j + i] = static_cast<uint32_t>(diff & 0xFFFFFFFFULL);
            borrow = static_cast<int64_t>(prod >> 32) - (diff >> 32);
        }
        int64_t diff_top = static_cast<int64_t>(u.digits_[j + n]) - borrow;
        u.digits_[j + n] = static_cast<uint32_t>(diff_top & 0xFFFFFFFFULL);

        // D5: 设置商字
        quotient.digits_[j] = static_cast<uint32_t>(qhat);

        // D6: 若借位（diff_top < 0），补充加回
        if (diff_top < 0) {
            quotient.digits_[j]--;
            uint64_t carry = 0;
            for (size_t i = 0; i < n; i++) {
                carry += static_cast<uint64_t>(u.digits_[j + i]) + v.digits_[i];
                u.digits_[j + i] = static_cast<uint32_t>(carry & 0xFFFFFFFFULL);
                carry >>= 32;
            }
            u.digits_[j + n] += static_cast<uint32_t>(carry);
        }
    }

    // D8: 反规范化余数
    quotient.trim();
    BigInt remainder = u >> s;
    // 余数只取前 n 个字
    if (remainder.digits_.size() > n) {
        remainder.digits_.resize(n);
    }
    remainder.trim();
    return {quotient, remainder};
}

BigInt BigInt::operator/(const BigInt &other) const {
    return divmod(other).first;
}

BigInt BigInt::operator%(const BigInt &other) const {
    return divmod(other).second;
}

// ============================================
// 位移运算
// ============================================

/**
 * @brief 逻辑左移 bits 位
 * @details 拆分为整字移位和比特移位两步:
 *          1. 整字移位: 新数组前 word_shift 个字为 0
 *          2. 比特移位: 每个字左移 bit_shift 位，溢出部分进入高一字
 */
BigInt BigInt::operator<<(int bits) const {
    if (bits == 0 || is_zero())
        return *this;

    int word_shift = bits / 32; // 整字偏移量
    int bit_shift = bits % 32; // 字内比特偏移量

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

BigInt &BigInt::operator<<=(int bits) {
    *this = *this << bits;
    return *this;
}

/**
 * @brief 逻辑右移 bits 位
 * @details 拆分为整字移位和比特移位两步:
 *          1. 整字移位: 跳过低 word_shift 个字
 *          2. 比特移位: 每个字右移 bit_shift 位，
 *             高一字的低位补入当前字的高位
 */
BigInt BigInt::operator>>(int bits) const {
    if (bits == 0)
        return *this;

    int word_shift = bits / 32;
    int bit_shift = bits % 32;

    // 若移位量超过总位长，结果为 0
    if (word_shift >= static_cast<int>(digits_.size()))
        return BigInt(0);

    BigInt result;
    size_t new_size = digits_.size() - word_shift;
    result.digits_.resize(new_size, 0);

    for (size_t i = 0; i < new_size; i++) {
        result.digits_[i] = digits_[i + word_shift] >> bit_shift;
        // 从高一字借入低位（注意 bit_shift 为 0 时不需要借位）
        if (bit_shift > 0 && i + word_shift + 1 < digits_.size()) {
            result.digits_[i] |= digits_[i + word_shift + 1] << (32 - bit_shift);
        }
    }
    result.trim();
    return result;
}

BigInt &BigInt::operator>>=(int bits) {
    *this = *this >> bits;
    return *this;
}

// ============================================
// 十进制输出
// ============================================

/**
 * @brief 转换为十进制字符串
 * @details 反复除以 10 取余数，收集各位数字后反转
 *          时间复杂度: O(n × D)，D 为十进制位数
 */
std::string BigInt::to_decimal_string() const {
    if (is_zero())
        return "0";

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

/** @brief 流输出运算符 */
std::ostream &operator<<(std::ostream &os, const BigInt &b) {
    os << b.to_decimal_string();
    return os;
}
