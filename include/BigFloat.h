/**
 * @file BigFloat.h
 * @brief 任意精度浮点数类
 *
 * @details
 * BigFloat 采用 **二进制浮点** 表示:
 *
 *   value = (-1)^sign × mantissa × 2^exponent
 *
 * 其中:
 *   - sign_          : bool，true 表示负数
 *   - mantissa_      : BigInt（无符号），存储二进制尾数
 *   - exponent_      : int64_t，二进制指数
 *   - precision_bits_ : int，目标精度（二进制位数）
 *
 * 规范化不变量:
 *   - mantissa_ 的位长始终 ≈ precision_bits_（normalize() 负责截断或扩展）
 *   - 零值时 mantissa_=0, exponent_=0, sign_=false
 *   - normalize() 使用四舍五入（round-to-nearest）
 *
 * 主要功能:
 *   - 构造: 从 int/int64_t/double/字符串/BigInt/分数构造
 *   - 算术: +, -, ×, ÷, 乘以 2^n
 *   - 数学函数: pi, e, sqrt, exp, ln, pow, factorial, half_factorial
 *   - 输出: 十进制字符串（定点或科学记数法）
 */
#pragma once
#include "BigInt.h"
#include <cstdint>
#include <string>

class BigFloat {
public:
  // ========================
  // 构造函数
  // ========================

  /** @brief 默认构造，值为 0，精度 256 位 */
  BigFloat();

  /**
   * @brief 从 int 构造
   * @param val    整数值
   * @param prec_bits 目标精度（二进制位数），默认 256
   */
  BigFloat(int val, int prec_bits = 256);

  /**
   * @brief 从 int64_t 构造
   * @param val    64 位有符号整数值
   * @param prec_bits 目标精度
   */
  BigFloat(int64_t val, int prec_bits = 256);

  /**
   * @brief 从 double 构造
   * @param val    双精度浮点数
   * @param prec_bits 目标精度
   * @details 使用 frexp 分解 double 为尾数和指数，保留 53 位精度
   */
  BigFloat(double val, int prec_bits = 256);

  /**
   * @brief 从分数 num/den 构造高精度浮点数
   * @param num      分子（BigInt，无符号）
   * @param den      分母（BigInt，无符号）
   * @param prec_bits 目标精度
   * @details 通过 (num << (prec_bits + 32)) / den 实现整数除法获取高精度商
   */
  static BigFloat from_fraction(const BigInt &num, const BigInt &den,
                                int prec_bits);

  /**
   * @brief 从 BigInt 构造（整数值，无符号）
   * @param val      大整数值
   * @param prec_bits 目标精度
   */
  static BigFloat from_bigint(const BigInt &val, int prec_bits);

  /**
   * @brief 从十进制字符串构造
   * @param str      十进制数字字符串，支持 "[-]123.456[e±789]" 格式
   * @param prec_bits 目标精度
   * @details 解析整数部分、小数部分和指数部分，
   *          利用 10^n = 5^n × 2^n 分解来保持二进制精度
   */
  static BigFloat from_string(const std::string &str, int prec_bits);

  /** @brief 复制/移动（编译器默认实现） */
  BigFloat(const BigFloat &) = default;

  BigFloat(BigFloat &&) noexcept = default;

  BigFloat &operator=(const BigFloat &) = default;

  BigFloat &operator=(BigFloat &&) noexcept = default;

  // ========================
  // 精度控制
  // ========================

  /** @brief 获取当前精度（二进制位数） */
  int precision() const { return precision_bits_; }

  /**
   * @brief 设置精度并立即 normalize
   * @param bits 新的精度（二进制位数）
   * @details 如果新精度小于当前尾数位长，会截断并四舍五入；
   *          如果新精度大于当前尾数位长，会左移扩展
   */
  void set_precision(int bits) {
    precision_bits_ = bits;
    normalize();
  }

  // ========================
  // 符号与零判断
  // ========================

  /** @brief 判断是否为负数 */
  bool is_negative() const { return sign_; }

  /** @brief 判断是否为零（尾数为零即为零值） */
  bool is_zero() const { return mantissa_.is_zero(); }

  /** @brief 取相反数 */
  BigFloat operator-() const;

  /** @brief 取绝对值 */
  BigFloat abs() const;

  // ========================
  // 取整与奇偶判断
  // ========================

  /**
   * @brief 向下取整
   * @return 返回不超过当前值的最大整数对应的 BigFloat
   */
  BigFloat floor() const;

  /**
   * @brief 判断该浮点数所代表的值是否为奇整数
   * @return 如果值是精确的整数且为奇数，则返回 true；否则 false
   */
  bool is_odd_integer() const;

  // ========================
  // 比较运算
  // ========================

  /**
   * @brief 三路比较
   * @return -1 / 0 / 1
   * @details 先比较符号，同号时对齐指数后比较尾数
   */
  int compare(const BigFloat &other) const;

  bool operator==(const BigFloat &other) const { return compare(other) == 0; }
  bool operator!=(const BigFloat &other) const { return compare(other) != 0; }
  bool operator<(const BigFloat &other) const { return compare(other) < 0; }
  bool operator<=(const BigFloat &other) const { return compare(other) <= 0; }
  bool operator>(const BigFloat &other) const { return compare(other) > 0; }
  bool operator>=(const BigFloat &other) const { return compare(other) >= 0; }

  // ========================
  // 算术运算
  // ========================

  BigFloat operator+(const BigFloat &other) const;

  BigFloat operator-(const BigFloat &other) const;

  BigFloat operator*(const BigFloat &other) const;

  /**
   * @brief 除法
   * @throws std::runtime_error 若 other 为 0
   * @details 通过 (mantissa << (prec + 32)) / other.mantissa 实现
   */
  BigFloat operator/(const BigFloat &other) const;

  BigFloat &operator+=(const BigFloat &other);

  BigFloat &operator-=(const BigFloat &other);

  BigFloat &operator*=(const BigFloat &other);

  BigFloat &operator/=(const BigFloat &other);

  /**
   * @brief 高效地乘以 2^n（只调整指数，不改变尾数）
   * @param n 指数偏移量（可为负数）
   */
  BigFloat mul_pow2(int n) const;

  /**
   * @brief 高效除以小整数（O(n) 单字除法）
   * @param val 除数（32 位无符号）
   * @details 直接对尾数做 BigInt::divmod_u32，避免构造 BigFloat 的 O(n²) 除法
   */
  BigFloat div_u32(uint32_t val) const;

  // ========================
  // 数学函数（静态方法）
  // ========================

  /**
   * @brief 计算圆周率 π
   * @param prec_bits 目标精度
   * @details 使用 Machin 公式: π/4 = 4·arctan(1/5) − arctan(1/239)
   */
  static BigFloat pi(int prec_bits);

  /**
   * @brief 计算自然常数 e
   * @param prec_bits 目标精度
   * @details 等价于 exp(1)
   */
  static BigFloat e(int prec_bits);

  /**
   * @brief 计算平方根
   * @param x 被开方数（必须 ≥ 0）
   * @throws std::runtime_error 若 x < 0
   * @details 牛顿迭代法: x_{n+1} = (x_n + a/x_n) / 2
   */
  static BigFloat sqrt(const BigFloat &x);

  /**
   * @brief 计算正弦函数 sin(x)
   * @param x 弧度值
   * @details 基于严格的 2π 模取余 + 三倍角区间缩减 +
   * 泰勒级数，保证收敛的绝对稳定性。
   */
  static BigFloat sin(const BigFloat &x);

  /**
   * @brief 计算指数函数 e^x
   * @param x 指数
   * @details 范围缩减 + 泰勒级数:
   *          1. 先将 x 缩小为 x/(2^r) 使其足够小
   *          2. 用泰勒级数计算 e^(x/2^r)
   *          3. 反向平方 r 次得到 e^x
   */
  static BigFloat exp(const BigFloat &x);

  /**
   * @brief 计算自然对数 ln(x)
   * @param x 真数（必须 > 0）
   * @throws std::runtime_error 若 x ≤ 0
   * @details 将 x 分解为 m × 2^e（m ∈ [1,2)），然后:
   *          ln(x) = ln(m) + e × ln(2)
   *          ln(m) = 2 × atanh((m-1)/(m+1))
   */
  static BigFloat ln(const BigFloat &x);

  /**
   * @brief 计算幂函数 base^exp
   * @details 整数指数 ≤ 1000 时使用二进制快速幂（避免 ln 的精度损失），
   *          否则使用通用路径 exp(exp × ln(base))
   */
  static BigFloat pow(const BigFloat &base, const BigFloat &exp);

  /**
   * @brief 计算阶乘 n!
   * @param n    非负整数
   * @param prec_bits 目标精度
   * @details 直接用 BigInt 逐次乘法计算，再转为 BigFloat
   */
  static BigFloat factorial(int n, int prec_bits);

  /**
   * @brief 计算半整数阶乘 (k − 0.5)! = Γ(k + 0.5)
   * @param k    非负整数
   * @param prec_bits 目标精度
   * @details 使用公式: Γ(k+0.5) = (2k−1)!! / 2^k × √π
   *          其中 (2k−1)!! = 1 × 3 × 5 × ⋯ × (2k−1)
   *          特别地，k=0 时 Γ(0.5) = √π
   */
  static BigFloat half_factorial(int k, int prec_bits);

  // ========================
  // I/O
  // ========================

  /**
   * @brief 转换为十进制字符串
   * @param decimal_digits 有效数字位数
   * @param scientific_notation 是否使用科学记数法（默认 true）
   * @details 通过将二进制浮点值乘以适当的 10 的幂转换为十进制整数，
   *          再格式化为字符串输出
   */
  std::string to_decimal_string(int decimal_digits,
                                bool scientific_notation = true) const;

  friend std::ostream &operator<<(std::ostream &os, const BigFloat &f);

  // ========================
  // 内部访问
  // ========================

  /** @brief 只读访问尾数 */
  const BigInt &mantissa() const { return mantissa_; }

  /** @brief 获取二进制指数 */
  int64_t exponent() const { return exponent_; }

private:
  bool sign_;          ///< 符号位: true = 负数, false = 非负
  BigInt mantissa_;    ///< 无符号尾数（二进制）
  int64_t exponent_;   ///< 二进制指数: value = mantissa × 2^exponent
  int precision_bits_; ///< 目标精度（尾数的二进制位数）

  /**
   * @brief 规范化: 将尾数截断或扩展到 precision_bits_ 位
   * @details 若尾数位长 > precision_bits_: 右移并四舍五入（round-to-nearest）
   *          若尾数位长 < precision_bits_: 左移补零
   *          零值时重置 exponent_ 和 sign_
   */
  void normalize();

  /**
   * @brief 对齐两个 BigFloat 的指数，使尾数可以直接加减
   * @param a, b       输入的两个 BigFloat
   * @param a_mant, b_mant 输出的对齐后尾数
   * @param common_exp 输出的公共指数
   * @details 将指数较大的数的尾数左移，使两个数共享同一指数。
   *          若指数差过大（超出精度范围），则小数视为零。
   */
  static void align(const BigFloat &a, const BigFloat &b, BigInt &a_mant,
                    BigInt &b_mant, int64_t &common_exp);
};
