/**
 * @file BigInt.h
 * @brief 无符号任意精度大整数类
 *
 * @details
 * BigInt 使用 vector<uint32_t> 作为内部存储，以 2^32 为基数，
 * 采用小端序排列（digits_[0] 存储最低 32 位）。
 *
 * 主要功能:
 *   - 支持从 uint64_t、十进制字符串构造
 *   - 支持加、减、乘、除、取模运算
 *   - 支持左移/右移位运算
 *   - 支持十进制字符串输出
 *
 * 注意: 本类仅处理**无符号**整数，减法要求 *this >= other，
 *       否则会触发 assert 失败。如需处理有符号整数，
 *       请在外部维护符号位（参见 lanczos.cpp 中的 SignedBigInt）。
 *
 * 时间复杂度:
 *   - 加/减法: O(n)，n 为数字长度（uint32_t 字数）
 *   - 乘法: O(n*m)，经典长乘法
 *   - 除法: O(n * bit_length)，二进制长除法
 */
#pragma once
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

class BigInt {
public:
  // ========================
  // 构造函数
  // ========================

  /** @brief 默认构造，值为 0 */
  BigInt();

  /**
   * @brief 从 uint64_t 构造
   * @param val 无符号 64 位整数值
   */
  BigInt(uint64_t val);

  /**
   * @brief 从十进制字符串构造
   * @param decimal_str 仅包含 '0'-'9' 的字符串，非数字字符会被跳过
   * @details 内部逐字符执行 result = result * 10 + digit
   */
  BigInt(const std::string &decimal_str);

  /** @brief 复制/移动构造与赋值（使用编译器默认实现） */
  BigInt(const BigInt &other) = default;

  BigInt(BigInt &&other) noexcept = default;

  BigInt &operator=(const BigInt &other) = default;

  BigInt &operator=(BigInt &&other) noexcept = default;

  // ========================
  // 比较运算
  // ========================

  /**
   * @brief 三路比较
   * @return -1 若 *this < other，0 若相等，1 若 *this > other
   * @details 先比较字长（digits_.size()），相同则从高位到低位逐字比较
   */
  int compare(const BigInt &other) const;

  bool operator==(const BigInt &other) const { return compare(other) == 0; }
  bool operator!=(const BigInt &other) const { return compare(other) != 0; }
  bool operator<(const BigInt &other) const { return compare(other) < 0; }
  bool operator<=(const BigInt &other) const { return compare(other) <= 0; }
  bool operator>(const BigInt &other) const { return compare(other) > 0; }
  bool operator>=(const BigInt &other) const { return compare(other) >= 0; }

  // ========================
  // 算术运算（返回新对象）
  // ========================

  /** @brief 无符号加法，O(n) */
  BigInt operator+(const BigInt &other) const;

  /**
   * @brief 无符号减法
   * @pre *this >= other，否则 assert 失败
   */
  BigInt operator-(const BigInt &other) const;

  /** @brief 经典长乘法，O(n*m) */
  BigInt operator*(const BigInt &other) const;

  /**
   * @brief 大数除法，返回 {商, 余数}
   * @throws std::runtime_error 若 divisor 为 0
   * @details 单字除数走快速路径 divmod_u32，否则使用二进制长除法
   */
  std::pair<BigInt, BigInt> divmod(const BigInt &divisor) const;

  BigInt operator/(const BigInt &other) const;

  BigInt operator%(const BigInt &other) const;

  /**
   * @brief 乘以单个 uint32_t，O(n) 快速路径
   * @param val 乘数（32 位无符号）
   */
  BigInt mul_u32(uint32_t val) const;

  /**
   * @brief 除以单个 uint32_t，返回 {商, 余数}
   * @throws std::runtime_error 若 val 为 0
   */
  std::pair<BigInt, uint32_t> divmod_u32(uint32_t val) const;

  // ========================
  // 复合赋值运算
  // ========================

  BigInt &operator+=(const BigInt &other);

  BigInt &operator-=(const BigInt &other);

  BigInt &operator*=(const BigInt &other);

  // ========================
  // 位移运算
  // ========================

  /**
   * @brief 逻辑左移
   * @param bits 左移位数，必须 >= 0
   * @details 分为整字移位（word_shift = bits/32）和比特移位（bit_shift =
   * bits%32）
   */
  BigInt operator<<(int bits) const;

  BigInt &operator<<=(int bits);

  /**
   * @brief 逻辑右移
   * @param bits 右移位数，必须 >= 0
   * @details 若 bits >= 总位长，返回 0
   */
  BigInt operator>>(int bits) const;

  BigInt &operator>>=(int bits);

  // ========================
  // 属性查询
  // ========================

  /** @brief 判断是否为 0 */
  bool is_zero() const;

  /**
   * @brief 返回最高有效位的位置（从 1 开始计数）
   * @return 0 表示值为 0，否则返回 floor(log2(val)) + 1
   */
  int bit_length() const;

  /**
   * @brief 获取指定位置的比特值
   * @param pos 位位置（从 0 开始，0 = 最低位）
   * @return true 表示该位为 1
   */
  bool get_bit(int pos) const;

  // ========================
  // I/O
  // ========================

  /**
   * @brief 转换为十进制字符串
   * @details 反复除以 10 取余数，最后反转字符串
   */
  std::string to_decimal_string() const;

  friend std::ostream &operator<<(std::ostream &os, const BigInt &b);

  // ========================
  // 内部访问
  // ========================

  /** @brief 只读访问内部数据（小端序 uint32_t 数组） */
  const std::vector<uint32_t> &digits() const { return digits_; }

  /** @brief 可写访问内部数据（慎用，调用者需自行维护不变量） */
  std::vector<uint32_t> &digits_mut() { return digits_; }

  /** @brief 去除高位多余的零（维护内部不变量：至少保留一个元素） */
  void trim();

private:
  /**
   * @brief 内部存储：小端序 uint32_t 数组
   * @details digits_[0] 是最低 32 位，digits_.back() 是最高 32 位
   *          不变量：digits_.size() >= 1，且除 digits_[0] 外不应有前导零
   */
  std::vector<uint32_t> digits_;

  /**
   * @brief Karatsuba 分治乘法（内部实现）
   * @details 当操作数超过 KARATSUBA_THRESHOLD 字时启用。
   *          复杂度 O(n^1.585)，相比经典 O(n²) 显著提升。
   */
  static BigInt karatsuba_mul(const BigInt &a, const BigInt &b);

  /// Karatsuba 启用阈值（字数）
  static constexpr size_t KARATSUBA_THRESHOLD = 32;
};
