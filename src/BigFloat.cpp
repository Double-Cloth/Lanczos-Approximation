/**
 * @file BigFloat.cpp
 * @brief BigFloat 类的实现 —— 任意精度浮点数
 *
 * @details
 * 表示: value = (-1)^sign × mantissa × 2^exponent
 *
 * 核心设计:
 *   - 所有运算使用额外保护位（guard bits）确保精度
 *   - normalize() 在每次运算后截断尾数并四舍五入
 *   - 数学函数（exp, ln, sqrt, pow）都使用纯 BigFloat 运算，
 *     不依赖 double 中间结果（除 sqrt 的初始估计外）
 */

#include "BigFloat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

// ============================================
// 构造函数
// ============================================

/**
 * @brief 默认构造: 值为 0，精度 256 位
 */
BigFloat::BigFloat()
    : sign_(false), mantissa_(0), exponent_(0), precision_bits_(256) {}

/**
 * @brief 从 int 构造
 * @details 取绝对值存入 mantissa_，符号单独记录。
 *          使用 int64_t 中转以安全处理 INT_MIN。
 */
BigFloat::BigFloat(int val, int prec_bits)
    : sign_(val < 0), mantissa_(static_cast<uint64_t>(
                          val < 0 ? -static_cast<int64_t>(val) : val)),
      exponent_(0), precision_bits_(prec_bits) {
  normalize();
}

/**
 * @brief 从 int64_t 构造
 */
BigFloat::BigFloat(int64_t val, int prec_bits)
    : sign_(val < 0), mantissa_(static_cast<uint64_t>(val < 0 ? -val : val)),
      exponent_(0), precision_bits_(prec_bits) {
  normalize();
}

/**
 * @brief 从 double 构造
 * @details 使用 std::frexp 将 double 分解为 frac × 2^exp2
 *          其中 0.5 ≤ frac < 1，则 frac × 2^53 是精确的 53 位整数尾数。
 *          exponent_ = exp2 - 53
 */
BigFloat::BigFloat(double val, int prec_bits)
    : sign_(val < 0), exponent_(0), precision_bits_(prec_bits) {
  if (val == 0.0) {
    mantissa_ = BigInt(0);
    return;
  }
  double abs_val = std::abs(val);
  int exp2;
  double frac = std::frexp(abs_val, &exp2); // abs_val = frac × 2^exp2
  // frac ∈ [0.5, 1)，乘以 2^53 得到 53 位整数尾数
  uint64_t mant = static_cast<uint64_t>(frac * (1ULL << 53));
  mantissa_ = BigInt(mant);
  exponent_ = exp2 - 53;
  normalize();
}

/**
 * @brief 从 BigInt 构造（结果为非负整数值）
 * @details exponent_ 初始为 0，normalize() 会调整尾数位长
 */
BigFloat BigFloat::from_bigint(const BigInt &val, int prec_bits) {
  BigFloat result;
  result.sign_ = false;
  result.mantissa_ = val;
  result.exponent_ = 0;
  result.precision_bits_ = prec_bits;
  result.normalize();
  return result;
}

/**
 * @brief 从分数 num/den 构造
 * @details 为保持精度，先将分子左移 (prec_bits + 32) 位再做整数除法:
 *          result = (num << (prec_bits + 32)) / den
 *          exponent = -(prec_bits + 32)
 *          额外 32 位作为保护位
 */
BigFloat BigFloat::from_fraction(const BigInt &num, const BigInt &den,
                                 int prec_bits) {
  int guard = std::max(32, prec_bits / 10);
  BigInt shifted = num << (prec_bits + guard); // 动态保护位
  auto [q, r] = shifted.divmod(den);
  BigFloat result;
  result.sign_ = false;
  result.mantissa_ = q;
  result.exponent_ = -(prec_bits + guard);
  result.precision_bits_ = prec_bits;
  result.normalize();
  return result;
}

// ============================================
// 规范化: 截断尾数到目标精度
// ============================================

/**
 * @brief 将尾数调整为恰好 precision_bits_ 位
 * @details
 *   - 零值: 重置 exponent_ 和 sign_
 *   - 尾数过长: 右移并四舍五入（检查被丢弃的最高位）
 *   - 尾数过短: 左移补零
 *
 *   四舍五入策略: round-to-nearest（检查被截断部分的最高位）
 *   注意: 四舍五入后可能导致尾数位长增加 1 位，需再次调整
 */
void BigFloat::normalize() {
  if (mantissa_.is_zero()) {
    exponent_ = 0;
    sign_ = false;
    return;
  }
  int bl = mantissa_.bit_length();
  if (bl > precision_bits_) {
    int shift = bl - precision_bits_;
    // 四舍五入: 检查被丢弃部分的最高位
    bool round_up = (shift > 0) && mantissa_.get_bit(shift - 1);
    mantissa_ >>= shift;
    exponent_ += shift;
    if (round_up) {
      mantissa_ += BigInt(1);
      // 进位可能使位长增加一位，需再截断
      if (mantissa_.bit_length() > precision_bits_) {
        mantissa_ >>= 1;
        exponent_++;
      }
    }
  } else if (bl < precision_bits_) {
    // 尾数过短: 左移扩展到目标精度
    int shift = precision_bits_ - bl;
    mantissa_ <<= shift;
    exponent_ -= shift;
  }
}

// ============================================
// 符号操作
// ============================================

/** @brief 取相反数（零值不改变符号） */
BigFloat BigFloat::operator-() const {
  BigFloat result = *this;
  if (!result.is_zero())
    result.sign_ = !result.sign_;
  return result;
}

/** @brief 取绝对值 */
BigFloat BigFloat::abs() const {
  BigFloat result = *this;
  result.sign_ = false;
  return result;
}

// ============================================
// 取整与奇偶判断
// ============================================

BigFloat BigFloat::floor() const {
  if (is_zero())
    return *this;
  int64_t shift = -exponent_;
  if (shift <= 0)
    return *this; // 已是整数

  int64_t int_bits = static_cast<int64_t>(mantissa_.bit_length()) + exponent_;
  if (int_bits <= 0) {
    if (sign_) {
      return BigFloat(-1, precision_bits_);
    } else {
      return BigFloat(0, precision_bits_);
    }
  }

  BigInt m = mantissa_ >> static_cast<int>(shift);
  if (sign_ && ((m << static_cast<int>(shift)) != mantissa_)) {
    m = m + BigInt(1);
  }

  BigFloat res = BigFloat::from_bigint(m, precision_bits_);
  res.sign_ = sign_;
  res.normalize();
  return res;
}

bool BigFloat::is_odd_integer() const {
  if (is_zero())
    return false;
  if (exponent_ > 0)
    return false;
  if (exponent_ == 0) {
    if (mantissa_.digits().empty())
      return false;
    return (mantissa_.digits()[0] & 1) != 0;
  }

  int64_t shift = -exponent_;
  int64_t int_bits = static_cast<int64_t>(mantissa_.bit_length()) + exponent_;
  if (int_bits <= 0)
    return false; // purely fractional

  BigInt m = mantissa_ >> static_cast<int>(shift);
  if (m.digits().empty())
    return false;
  return (m.digits()[0] & 1) != 0;
}

// ============================================
// 指数对齐: 使两个 BigFloat 可以直接做尾数加减
// ============================================

/**
 * @brief 对齐两个 BigFloat 的指数
 * @param a, b          两个输入 BigFloat
 * @param a_mant, b_mant 输出: 对齐后的尾数
 * @param common_exp     输出: 公共指数
 *
 * @details 原理:
 *   若 a.exp < b.exp，则 diff = b.exp - a.exp
 *   将 b 的尾数左移 diff 位，使指数统一为 a.exp:
 *     a = a.mant × 2^(a.exp)
 *     b = (b.mant << diff) × 2^(a.exp)
 *
 *   若 diff 超过最大精度 + 32 位保护位，说明一方远大于另一方，
 *   可以安全地将小数视为零（避免超大位移导致内存问题）。
 */
void BigFloat::align(const BigFloat &a, const BigFloat &b, BigInt &a_mant,
                     BigInt &b_mant, int64_t &common_exp) {
  // 处理零值: 防止无意义的位移
  if (a.is_zero()) {
    a_mant = BigInt(0);
    b_mant = b.mantissa_;
    common_exp = b.exponent_;
    return;
  }
  if (b.is_zero()) {
    a_mant = a.mantissa_;
    b_mant = BigInt(0);
    common_exp = a.exponent_;
    return;
  }

  if (a.exponent_ < b.exponent_) {
    int64_t diff = b.exponent_ - a.exponent_;
    int max_prec =
        std::max(a.precision_bits_, b.precision_bits_) + 32; // 32 位保护位

    if (diff > max_prec) {
      // b 远大于 a，a 相对可忽略
      a_mant = BigInt(0);
      b_mant = b.mantissa_;
      common_exp = b.exponent_;
    } else {
      // 正常对齐: 左移 b 的尾数，指数统一为较小的 a.exponent_
      a_mant = a.mantissa_;
      b_mant = b.mantissa_ << static_cast<int>(diff);
      common_exp = a.exponent_;
    }
  } else {
    int64_t diff = a.exponent_ - b.exponent_;
    int max_prec = std::max(a.precision_bits_, b.precision_bits_) + 32;

    if (diff > max_prec) {
      // a 远大于 b，b 相对可忽略
      a_mant = a.mantissa_;
      b_mant = BigInt(0);
      common_exp = a.exponent_;
    } else {
      // 正常对齐: 左移 a 的尾数，指数统一为较小的 b.exponent_
      b_mant = b.mantissa_;
      a_mant = a.mantissa_ << static_cast<int>(diff);
      common_exp = b.exponent_;
    }
  }
}

// ============================================
// 从十进制字符串构造
// ============================================

/**
 * @brief 解析十进制字符串为 BigFloat
 * @param str       支持格式: "[-]123.456[e±789]"
 * @param prec_bits 目标精度
 *
 * @details 解析步骤:
 *   1. 解析符号、整数部分、小数部分、指数部分
 *   2. combined = int_part + frac_part（合并为无小数点的整数字符串）
 *   3. dec_exp = 原始指数 − 小数部分长度（十进制位移量）
 *      即 value = combined × 10^dec_exp
 *   4. 利用 10^n = 5^n × 2^n 将十进制转为二进制:
 *      - dec_exp > 0: mantissa *= 5^dec_exp, exponent += dec_exp
 *      - dec_exp < 0: mantissa = mantissa × 2^K / 5^|dec_exp|,
 *                     exponent -= (|dec_exp| + K)
 *        其中 K 选择足够大以保持精度
 */
BigFloat BigFloat::from_string(const std::string &str, int prec_bits) {
  if (str.empty())
    return BigFloat(0, prec_bits);

  // --- 解析符号 ---
  size_t i = 0; // 使用 size_t 避免与 str.size() 比较时的符号警告
  bool sign = false;
  if (str[i] == '-') {
    sign = true;
    i++;
  } else if (str[i] == '+') {
    i++;
  }

  // --- 解析整数部分和小数部分 ---
  std::string int_part;
  std::string frac_part;
  int exp_val = 0;

  // 整数部分: 连续数字
  while (i < str.size() && std::isdigit(str[i])) {
    int_part += str[i++];
  }
  // 小数部分: '.' 后的连续数字
  if (i < str.size() && str[i] == '.') {
    i++;
    while (i < str.size() && std::isdigit(str[i])) {
      frac_part += str[i++];
    }
  }
  // 指数部分: 'e' 或 'E' 后的整数
  if (i < str.size() && (str[i] == 'e' || str[i] == 'E')) {
    i++;
    std::string exp_str;
    while (i < str.size()) {
      exp_str += str[i++];
    }
    if (!exp_str.empty()) {
      exp_val = std::stoi(exp_str);
    }
  }

  // --- 合并整数和小数部分 ---
  std::string combined = int_part + frac_part;

  // 去除前导零（保留至少一个字符）
  size_t start = 0;
  while (start < combined.size() - 1 && combined[start] == '0')
    start++;
  combined = combined.substr(start);

  if (combined.empty() || combined == "0")
    return BigFloat(0, prec_bits);

  // --- 计算十进制指数偏移 ---
  // value = combined × 10^dec_exp
  int dec_exp = exp_val - static_cast<int>(frac_part.size());
  BigInt num(combined);

  // --- 转换为二进制浮点 ---
  int guard = std::max(64, prec_bits / 10);
  int work_prec = prec_bits + guard; // 动态保护位
  BigFloat result = from_bigint(num, work_prec);

  // 快速幂计算 5^p（利用二进制分解）
  auto power5 = [](int p) {
    BigInt res(1);
    BigInt base(5);
    while (p > 0) {
      if (p & 1)
        res = res * base;
      base = base * base;
      p >>= 1;
    }
    return res;
  };

  if (dec_exp > 0) {
    // value = combined × 5^dec_exp × 2^dec_exp
    BigInt five_pow = power5(dec_exp);
    result.mantissa_ = result.mantissa_ * five_pow;
    result.exponent_ += dec_exp;
  } else if (dec_exp < 0) {
    // value = combined / 5^|dec_exp| / 2^|dec_exp|
    // 为保持精度: 先乘以 2^K 再除以 5^|dec_exp|
    int abs_exp = -dec_exp;
    BigInt five_pow = power5(abs_exp);

    // 选择 K 使得尾数除以 5^|dec_exp| 后仍有足够精度
    int bl_num = result.mantissa_.bit_length();
    int bl_five = five_pow.bit_length();
    int K = work_prec + bl_five - bl_num + 32;
    if (K < 0)
      K = 0;

    result.mantissa_ = (result.mantissa_ << K) / five_pow;
    result.exponent_ -= (abs_exp + K);
  }

  result.sign_ = sign;
  result.set_precision(prec_bits);
  return result;
}

// ============================================
// 比较
// ============================================

/**
 * @brief 三路比较
 * @details 优先级: 零值判断 → 符号比较 → 同号时对齐指数后比较尾数
 *          负数时结果取反
 */
int BigFloat::compare(const BigFloat &other) const {
  // 两个零相等
  if (is_zero() && other.is_zero())
    return 0;
  // 一方为零: 非零方的符号决定结果
  if (is_zero())
    return other.sign_ ? 1 : -1;
  if (other.is_zero())
    return sign_ ? -1 : 1;
  // 符号不同: 正数 > 负数
  if (sign_ != other.sign_)
    return sign_ ? -1 : 1;

  // 同号: 对齐指数后比较尾数
  BigInt a_mant, b_mant;
  int64_t common_exp;
  align(*this, other, a_mant, b_mant, common_exp);
  int cmp = a_mant.compare(b_mant);
  // 若为负数，比较结果取反
  return sign_ ? -cmp : cmp;
}

// ============================================
// 加法与减法
// ============================================

/**
 * @brief 加法
 * @details 对齐指数后:
 *   - 同号: 尾数直接相加
 *   - 异号: 尾数相减（大减小），结果符号跟随绝对值较大者
 *   结果精度取两者的较大值
 */
BigFloat BigFloat::operator+(const BigFloat &other) const {
  // 快速路径: 一方为零
  if (is_zero()) {
    BigFloat r = other;
    r.precision_bits_ = std::max(precision_bits_, other.precision_bits_);
    return r;
  }
  if (other.is_zero()) {
    BigFloat r = *this;
    r.precision_bits_ = std::max(precision_bits_, other.precision_bits_);
    return r;
  }

  BigFloat result;
  result.precision_bits_ = std::max(precision_bits_, other.precision_bits_);

  // 对齐指数
  BigInt a_mant, b_mant;
  int64_t common_exp;
  align(*this, other, a_mant, b_mant, common_exp);

  if (sign_ == other.sign_) {
    // 同号: 直接相加
    result.sign_ = sign_;
    result.mantissa_ = a_mant + b_mant;
  } else {
    // 异号: 做减法
    int cmp = a_mant.compare(b_mant);
    if (cmp == 0)
      return BigFloat(0, result.precision_bits_); // 完全抵消
    if (cmp > 0) {
      result.sign_ = sign_; // |a| > |b|, 符号跟随 a
      result.mantissa_ = a_mant - b_mant;
    } else {
      result.sign_ = other.sign_; // |b| > |a|, 符号跟随 b
      result.mantissa_ = b_mant - a_mant;
    }
  }
  result.exponent_ = common_exp;
  result.normalize();
  return result;
}

/**
 * @brief 减法: 通过 a - b = a + (-b) 实现
 */
BigFloat BigFloat::operator-(const BigFloat &other) const {
  BigFloat neg = -other;
  return *this + neg;
}

BigFloat &BigFloat::operator+=(const BigFloat &other) {
  *this = *this + other;
  return *this;
}

BigFloat &BigFloat::operator-=(const BigFloat &other) {
  *this = *this - other;
  return *this;
}

// ============================================
// 乘法
// ============================================

/**
 * @brief 乘法: 尾数相乘，指数相加，符号异或
 * @details result.mantissa = a.mantissa × b.mantissa
 *          result.exponent = a.exponent + b.exponent
 *          result.sign = a.sign ⊕ b.sign
 */
BigFloat BigFloat::operator*(const BigFloat &other) const {
  if (is_zero() || other.is_zero()) {
    return BigFloat(0, std::max(precision_bits_, other.precision_bits_));
  }
  BigFloat result;
  result.precision_bits_ = std::max(precision_bits_, other.precision_bits_);
  result.sign_ = sign_ != other.sign_;            // 异号为负
  result.mantissa_ = mantissa_ * other.mantissa_; // BigInt 乘法
  result.exponent_ = exponent_ + other.exponent_; // 指数相加
  result.normalize();
  return result;
}

BigFloat &BigFloat::operator*=(const BigFloat &other) {
  *this = *this * other;
  return *this;
}

// ============================================
// 除法: 通过整数除法实现
// ============================================

/**
 * @brief 除法: a / b
 * @throws std::runtime_error 若 b 为 0
 *
 * @details 为保持精度，将被除数尾数左移 (prec + 32) 位后做整数除法:
 *   result.mantissa = (a.mantissa << shift) / b.mantissa
 *   result.exponent = a.exponent - b.exponent - shift
 *   result.sign = a.sign ⊕ b.sign
 */
BigFloat BigFloat::operator/(const BigFloat &other) const {
  if (other.is_zero())
    throw std::runtime_error("BigFloat division by zero");
  if (is_zero())
    return BigFloat(0, std::max(precision_bits_, other.precision_bits_));

  int prec = std::max(precision_bits_, other.precision_bits_);

  // 为保持精度，先将被除数尾数左移足够多的位
  int guard = std::max(32, prec / 10);
  int shift = prec + guard;
  BigInt num = mantissa_ << shift;
  auto [q, r] = num.divmod(other.mantissa_);

  BigFloat result;
  result.precision_bits_ = prec;
  result.sign_ = sign_ != other.sign_;
  result.mantissa_ = q;
  result.exponent_ = exponent_ - other.exponent_ - shift;
  result.normalize();
  return result;
}

BigFloat &BigFloat::operator/=(const BigFloat &other) {
  *this = *this / other;
  return *this;
}

/**
 * @brief 高效地乘以 2^n（只调整指数，O(1) 操作）
 * @param n 可为正数（乘以 2^n）或负数（除以 2^|n|）
 */
BigFloat BigFloat::mul_pow2(int n) const {
  BigFloat result = *this;
  result.exponent_ += n;
  return result;
}

/**
 * @brief 高效除以小整数（O(n) 单字除法）
 * @details 直接对尾数做 BigInt::divmod_u32，避免构造 BigFloat 的 O(n²) 除法
 *          这比 operator/(BigFloat(k)) 快得多，因为后者需要：
 *          1. 构造 BigFloat(k)
 *          2. 尾数左移 (prec+32) 位
 *          3. BigInt divmod（即使是 Knuth D 也是 O(n)）
 *          4. normalize
 *          而 div_u32 只需要 BigInt 的单字除法 + normalize
 */
BigFloat BigFloat::div_u32(uint32_t val) const {
  if (val == 0)
    throw std::runtime_error("BigFloat div_u32: division by zero");
  if (is_zero())
    return *this;
  if (val == 1)
    return *this;

  BigFloat result;
  result.sign_ = sign_;
  result.precision_bits_ = precision_bits_;

  // 为保持精度，先将尾数左移一定保护位再除
  int guard = std::max(32, precision_bits_ / 10);
  BigInt shifted = mantissa_ << guard;
  auto [q, r] = shifted.divmod_u32(val);
  result.mantissa_ = q;
  result.exponent_ = exponent_ - guard;
  result.normalize();
  return result;
}

// ============================================
// 数学常数与函数
// ============================================

// --- arctan(1/x): 用于 Machin 公式计算 π ---

/**
 * @brief 计算 arctan(1/x) 的级数展开
 *
 * @details 使用 arctan(1/x) 的 Leibniz 级数:
 *   arctan(1/x) = 1/x − 1/(3x³) + 1/(5x⁵) − ...
 *               = Σ_{k=0}^∞ (−1)^k / ((2k+1) × x^(2k+1))
 *
 *   每项 = power / (2k+1)，每次 power /= x²
 *   收敛条件: term 的有效位低于精度阈值时停止
 */
static BigFloat arctan_recip(int x, int prec_bits) {
  BigFloat result(0, prec_bits);
  BigFloat x_bf(x, prec_bits);
  BigFloat x2 = x_bf * x_bf;                      // x²
  BigFloat power = BigFloat(1, prec_bits) / x_bf; // 初始 power = 1/x

  for (int k = 0;; k++) {
    BigFloat term = power.div_u32(static_cast<uint32_t>(2 * k + 1));

    // 收敛检查: 当 term 的有效大小远小于精度要求时停止
    int guard = std::max(16, prec_bits / 10);
    if (term.mantissa().bit_length() + term.exponent() <
        -static_cast<int64_t>(prec_bits + guard)) {
      break;
    }

    // 交替加减
    if (k % 2 == 0)
      result += term; // 偶数项: 加
    else
      result -= term; // 奇数项: 减

    power /= x2; // power = power / x²，即 power = 1/x^(2k+3)
  }
  return result;
}

/**
 * @brief 计算圆周率 π
 * @details 使用 Machin 公式:
 *   π/4 = 4 × arctan(1/5) − arctan(1/239)
 *
 *   此公式收敛快速（arctan(1/5) 和 arctan(1/239) 的级数收敛都很快），
 *   适合中等精度计算。对于超高精度（>10000 位），
 *   可考虑 Chudnovsky 算法。
 */
BigFloat BigFloat::pi(int prec_bits) {
  // 用动态保护位计算
  int guard = std::max(32, prec_bits / 10);
  int work_prec = prec_bits + guard;
  BigFloat a = arctan_recip(5, work_prec);
  BigFloat b = arctan_recip(239, work_prec);
  BigFloat four(4, work_prec);
  // π = 4 × (4 × arctan(1/5) − arctan(1/239))
  BigFloat result = (four * a - b) * four;
  result.set_precision(prec_bits);
  return result;
}

/**
 * @brief 计算自然常数 e = exp(1)
 */
BigFloat BigFloat::e(int prec_bits) { return exp(BigFloat(1, prec_bits)); }

// ============================================
// sqrt: 牛顿迭代法
// ============================================

/**
 * @brief 计算 √x
 * @throws std::runtime_error 若 x < 0
 *
 * @details 牛顿迭代法:
 *   x_{n+1} = (x_n + a/x_n) / 2
 *
 *   收敛性: 二次收敛，每次迭代有效位数翻倍
 *   初始估计: 使用 double 精度的 ldexp 给出粗略值
 *   终止条件: |x_{n+1} − x_n| 的有效位远小于精度要求
 */
BigFloat BigFloat::sqrt(const BigFloat &x) {
  if (x.is_zero())
    return BigFloat(0, x.precision());
  if (x.is_negative())
    throw std::runtime_error("sqrt of negative number");

  int prec = x.precision();
  int guard = std::max(64, prec / 10);
  int work_prec = prec + guard; // 动态保护位

  // --- 初始估计: 使用 double 给出粗略的量级 ---
  int bl = x.mantissa().bit_length();
  int total_exp = bl + static_cast<int>(x.exponent());
  // √(2^total_exp) ≈ 2^(total_exp/2)
  double approx = std::ldexp(1.0, total_exp / 2);

  BigFloat guess(approx, work_prec);
  // 预计算 0.5
  BigFloat half(0, work_prec);
  {
    BigFloat one(1, work_prec);
    BigFloat two(2, work_prec);
    half = one / two;
  }
  BigFloat val = x;
  val.set_precision(work_prec);

  // --- 牛顿迭代 ---
  int iters = 0;
  int max_iters = prec; // 最坏情况（实际远少于此）
  while (iters < max_iters) {
    // x_{n+1} = (x_n + val/x_n) × 0.5
    BigFloat next = (guess + val / guess) * half;

    // 检查收敛: |next − guess| 足够小
    BigFloat diff = (next - guess).abs();
    guess = next;
    iters++;

    if (diff.is_zero())
      break;

    // 若 diff 相对于 guess 已经小于精度要求，则收敛
    int diff_bl =
        diff.mantissa().bit_length() + static_cast<int>(diff.exponent());
    int guess_bl =
        guess.mantissa().bit_length() + static_cast<int>(guess.exponent());
    int guard = std::max(8, prec / 32);
    if (guess_bl - diff_bl > prec + guard)
      break;
  }
  guess.set_precision(prec);
  return guess;
}

// ============================================
// sin: 周期缩减 + 三倍角降阶 + 泰勒级数
// ============================================

/**
 * @brief 计算 sin(x)
 *
 * @details 算法分四步:
 *
 *   1. **周期缩减**: 利用 2π 周期性，将任意大的输入 x 压缩至 [-π/2, π/2] 内，
 *      这有效避免了直接吃大坐标导致的系统性爆级和误差。
 *   2. **范围缩减**: 对于压缩后的结果再次除以 3^r，并利用 O(n) 单字运算。
 *   3. **泰勒级数**: 对极小数展开泰勒累加。
 *   4. **反向补偿**: 利用三次公式回退 r 轮得到最终精准结果。
 */
BigFloat BigFloat::sin(const BigFloat &x) {
  if (x.is_zero())
    return BigFloat(0, x.precision());

  int prec = x.precision();
  int guard = std::max(128, prec / 8);
  int work_prec = prec + guard; // 动态保护位

  BigFloat val = x;
  val.set_precision(work_prec);

  // 符号处理
  bool negate = val.is_negative();
  val = val.abs();

  // --- 步骤 1: 周期缩减 (Range Reduction Module) ---
  // 将 val 根据 pi 取模到 [-pi/2, pi/2]
  BigFloat pi_val = BigFloat::pi(work_prec);
  BigFloat half_pi = pi_val.mul_pow2(-1); // pi / 2

  if (val > half_pi) {
    // n = floor(val/pi + 0.5) 取最近整数
    BigFloat ratio = val / pi_val;
    BigFloat half(0, work_prec);
    {
      BigFloat one(1, work_prec);
      half = one.mul_pow2(-1);
    }
    BigFloat n_float = (ratio + half).floor();

    // 取余: val = val - n * pi
    val = val - n_float * pi_val;

    // 若 n 为奇数，sin(x + n*pi) = -sin(x)
    if (n_float.is_odd_integer()) {
      negate = !negate;
    }
  }

  // 若存在因漂浮计算变负的极其微小的误差
  if (val.is_negative()) {
    val = val.abs();
    negate = !negate;
  }
  if (val.is_zero())
    return BigFloat(0, prec);

  // --- 步骤 2: 三倍角降阶缩减 ---
  // 获取缩小倍率 r，使 t 减至足够小加速泰勒级数收敛（无限拓展）
  int r = std::max(4, prec / 32);

  BigFloat reduced = val;
  for (int i = 0; i < r; i++) {
    reduced = reduced.div_u32(3);
  }

  // --- 步骤 3: 泰勒级数计算 ---
  // sin(t) = t - t^3/3! + t^5/5! - ...
  BigFloat sum = reduced;
  BigFloat term = reduced;
  BigFloat reduced_sq = reduced * reduced;

  for (int k = 1;; k++) {
    term *= reduced_sq;
    // 除以 2k * (2k+1)
    uint32_t divisor = static_cast<uint32_t>(2 * k * (2 * k + 1));
    term = term.div_u32(divisor);

    if (term.is_zero())
      break;

    if (k % 2 == 1)
      sum -= term;
    else
      sum += term;

    int term_bl =
        term.mantissa().bit_length() + static_cast<int>(term.exponent());
    int sum_bl = sum.mantissa().bit_length() + static_cast<int>(sum.exponent());
    if (sum_bl - term_bl > work_prec)
      break;
  }

  // --- 步骤 4: 补偿还原 ---
  // 依靠三倍角公式: sin(3t) = 3sin(t) - 4sin^3(t) 还原 r 次
  for (int i = 0; i < r; i++) {
    BigFloat sum_cubed = sum * sum * sum;
    // 3*sum - 4*sum^3 = (sum + 2*sum) - 4*sum^3 = (sum + sum<<1) -
    // (sum_cubed<<2)
    BigFloat term1 = sum + sum.mul_pow2(1);
    BigFloat term2 = sum_cubed.mul_pow2(2);
    sum = term1 - term2;
  }

  if (negate)
    sum = -sum;
  sum.set_precision(prec);
  return sum;
}

// ============================================
// exp: 范围缩减 + 泰勒级数
// ============================================

/**
 * @brief 计算 e^x
 *
 * @details 算法分三步:
 *
 *   1. **范围缩减**: 将 x 除以 2^r，使 x/(2^r) 足够小（接近 0）
 *      这样泰勒级数可以快速收敛。r 的选择基于 x 的量级。
 *
 *   2. **泰勒级数**: 计算 e^(x/2^r)
 *      e^t = Σ_{k=0}^∞ t^k / k! = 1 + t + t²/2! + t³/3! + ...
 *      每项 term_k = term_{k-1} × t / k（递推避免重复计算幂和阶乘）
 *
 *   3. **反向平方**: 将结果平方 r 次
 *      e^x = (e^(x/2^r))^(2^r)
 *
 *   对于负指数，先计算 e^|x|，再取倒数: e^(-x) = 1/e^x
 *
 *   保护位: 额外 128 位以抵消 r 次平方累积的误差
 */
BigFloat BigFloat::exp(const BigFloat &x) {
  if (x.is_zero())
    return BigFloat(1, x.precision());

  int prec = x.precision();
  int guard = std::max(128, prec / 8);
  int work_prec = prec + guard; // 保护位抵消平方累积误差

  BigFloat val = x;
  val.set_precision(work_prec);
  bool negate = val.is_negative(); // 记录是否为负指数
  if (negate)
    val = val.abs();

  // --- 步骤 1: 范围缩减 ---
  // 将 val 除以 2^r 使其足够小，以加速泰勒级数收敛
  int r = 0;
  {
    int bl = val.mantissa().bit_length();
    int total = bl + static_cast<int>(val.exponent());
    if (total > 0) {
      r = total + 8; // val > 1 时需要较大的缩减
    } else {
      r = std::max(0, 8 + total); // val < 1 时适度缩减
    }
  }
  BigFloat reduced = val.mul_pow2(-r); // reduced = val / 2^r

  // --- 步骤 2: 泰勒级数 ---
  // e^reduced = 1 + reduced + reduced²/2! + reduced³/3! + ...
  // 优化: 用 div_u32(k) 替代 operator/(BigFloat(k))，避免每次迭代都做完整
  // BigFloat 除法
  BigFloat sum(1, work_prec);
  BigFloat term(1, work_prec);
  for (int k = 1;; k++) {
    term *= reduced;
    term = term.div_u32(static_cast<uint32_t>(k)); // O(n) 单字除法
    sum += term;

    if (term.is_zero())
      break;
    int term_bl =
        term.mantissa().bit_length() + static_cast<int>(term.exponent());
    int sum_bl = sum.mantissa().bit_length() + static_cast<int>(sum.exponent());
    if (sum_bl - term_bl > work_prec)
      break;
  }

  // --- 步骤 3: 反向平方 ---
  // result = sum^(2^r)，通过连续 r 次平方实现
  for (int i = 0; i < r; i++) {
    sum = sum * sum;
  }

  // 负指数取倒数
  if (negate) {
    sum = BigFloat(1, work_prec) / sum;
  }

  sum.set_precision(prec);
  return sum;
}

// ============================================
// ln: range reduction + atanh 级数
// ============================================

/**
 * @brief 计算自然对数 ln(x)
 * @throws std::runtime_error 若 x ≤ 0
 *
 * @details 算法:
 *
 *   1. **范围缩减**: 分解 x = m × 2^e，其中 m ∈ [1, 2)
 *      ln(x) = ln(m) + e × ln(2)
 *
 *   2. **计算 ln(m)**: 使用恒等式
 *      ln(m) = 2 × atanh((m−1)/(m+1))
 *      其中 m ∈ [1,2) 时 (m−1)/(m+1) ∈ [0, 1/3)
 *
 *   3. **atanh 级数**:
 *      atanh(t) = t + t³/3 + t⁵/5 + ...
 *      由于 t ∈ [0, 1/3)，级数收敛很快
 *
 *   4. **计算 ln(2)**: 同样使用 atanh
 *      ln(2) = 2 × atanh(1/3)
 *      （仅在 e ≠ 0 时计算，避免不必要的开销）
 */
BigFloat BigFloat::ln(const BigFloat &x) {
  if (x.is_zero() || x.is_negative())
    throw std::runtime_error("ln of non-positive number");

  int prec = x.precision();
  int guard = std::max(128, prec / 8);
  int work_prec = prec + guard; // 动态保护位

  // --- 步骤 1: 范围缩减 ---
  // 分解 x = m × 2^e，使 m ∈ [1, 2)
  BigFloat val = x;
  val.set_precision(work_prec);
  int bl = val.mantissa().bit_length();
  int64_t total_exp_bits = static_cast<int64_t>(bl) + val.exponent();

  // 构造 m: 将 val 缩放到 [1, 2) 范围
  // m = val × 2^(−total_exp_bits + 1)
  BigFloat m = val.mul_pow2(static_cast<int>(-total_exp_bits + 1));
  int64_t e_val = total_exp_bits - 1; // val ≈ m × 2^e_val

  // --- 步骤 2: 计算 ln(m) ---
  // ln(m) = 2 × atanh((m−1)/(m+1))
  BigFloat one(1, work_prec);
  BigFloat t = (m - one) / (m + one); // t ∈ [0, 1/3)
  BigFloat t2 = t * t;                // t² 用于递推

  // atanh(t) = t + t³/3 + t⁵/5 + ...
  BigFloat sum = t;
  BigFloat power = t; // power = t^(2k+1)
  for (int k = 1;; k++) {
    power *= t2; // power = t^(2k+1)
    BigFloat term = power.div_u32(static_cast<uint32_t>(2 * k + 1));
    sum += term;

    // 收敛检查
    if (term.is_zero())
      break;
    int term_bl =
        term.mantissa().bit_length() + static_cast<int>(term.exponent());
    int sum_bl = sum.mantissa().bit_length() + static_cast<int>(sum.exponent());
    if (sum_bl - term_bl > work_prec)
      break;
  }

  BigFloat ln_m = sum * BigFloat(2, work_prec); // ln(m) = 2 × atanh(t)

  // --- 步骤 3: 计算 ln(2)（缓存以避免重复计算） ---
  BigFloat ln2(0, work_prec);
  if (e_val != 0) {
    // 缓存 ln(2) — 同一精度下复用，避免每次 ln() 都重算 atanh(1/3) 级数
    static BigFloat cached_ln2;
    static int cached_ln2_prec = 0;
    if (cached_ln2_prec < work_prec) {
      // ln(2) = 2 × atanh(1/3)
      BigFloat one_wp(1, work_prec);
      BigFloat third = one_wp / BigFloat(3, work_prec);
      BigFloat third2 = third * third;
      BigFloat ln2_sum = third;
      BigFloat ln2_power = third;
      for (int k = 1;; k++) {
        ln2_power *= third2;
        BigFloat term = ln2_power.div_u32(static_cast<uint32_t>(2 * k + 1));
        ln2_sum += term;

        if (term.is_zero())
          break;
        int term_bl =
            term.mantissa().bit_length() + static_cast<int>(term.exponent());
        int sum_bl = ln2_sum.mantissa().bit_length() +
                     static_cast<int>(ln2_sum.exponent());
        if (sum_bl - term_bl > work_prec)
          break;
      }
      cached_ln2 = ln2_sum * BigFloat(2, work_prec);
      cached_ln2_prec = work_prec;
    }
    ln2 = cached_ln2;
    ln2.set_precision(work_prec);
  }

  // --- 步骤 4: 合并结果 ---
  // ln(x) = ln(m) + e × ln(2)
  BigFloat result =
      ln_m + BigFloat(static_cast<int64_t>(e_val), work_prec) * ln2;
  result.set_precision(prec);
  return result;
}

// ============================================
// pow: 整数快速幂 / 通用路径
// ============================================

/**
 * @brief 计算 base^exponent
 *
 * @details 两条路径:
 *
 *   1. **整数快速幂**（exponent 为非负整数且 ≤ 1000 时）:
 *      使用二进制分解: result = base^(b_0 + 2×b_1 + 4×b_2 + ...)
 *      每步 result *= base（当对应位为 1）, base *= base
 *      优点: 无需 ln 运算，精度无损失
 *
 *   2. **通用路径**:
 *      base^exp = exp(exp × ln(base))
 *      适用于非整数指数或大指数
 *      注意: ln 和 exp 会引入额外误差，保护位很重要
 */
BigFloat BigFloat::pow(const BigFloat &base, const BigFloat &exponent) {
  if (exponent.is_zero())
    return BigFloat(1, base.precision());
  if (base.is_zero())
    return BigFloat(0, base.precision());

  int prec = std::max(base.precision(), exponent.precision());
  int guard = std::max(128, prec / 8);
  int work_prec = prec + guard;

  // --- 检查指数是否为非负整数 ---
  // 条件: exponent ≥ 0 且 exponent_ ≥ 0（没有小数部分）
  if (!exponent.is_negative() && exponent.exponent() >= 0) {
    // 直接从 BigInt 提取整数值，避免 to_decimal_string + stoll 的 O(n×D) 开销
    BigInt mant = exponent.mantissa();
    int64_t shift = exponent.exponent();
    if (shift > 0 && shift <= 64) {
      mant <<= static_cast<int>(shift);
    }
    int64_t total_bits = static_cast<int64_t>(mant.bit_length());
    if (total_bits <= 31) {
      // 直接从 digits_ 提取整数值，无需十进制转换
      uint64_t n_val = mant.digits()[0];
      if (n_val <= 1000) {
        // 使用二进制快速幂: O(log(n)) 次乘法
        BigFloat result(1, work_prec);
        BigFloat b = base;
        b.set_precision(work_prec);
        int64_t p = n_val;
        while (p > 0) {
          if (p & 1)
            result *= b; // 当前位为 1 时乘入结果
          b *= b;        // base 平方
          p >>= 1;
        }
        result.set_precision(prec);
        return result;
      }
    }
  }

  // --- 通用路径: exp(exponent × ln(base)) ---
  BigFloat ln_base = ln(base.abs());
  ln_base.set_precision(work_prec);
  BigFloat exp_val = exponent;
  exp_val.set_precision(work_prec);

  return exp(exp_val * ln_base);
}

// ============================================
// 阶乘与半整数阶乘
// ============================================

/**
 * @brief 计算整数阶乘 n!
 * @details 用 BigInt 逐次乘法计算（精确无损），再转为 BigFloat
 */
BigFloat BigFloat::factorial(int n, int prec_bits) {
  BigInt result(1);
  for (int i = 2; i <= n; i++) {
    result = result.mul_u32(static_cast<uint32_t>(i));
  }
  return from_bigint(result, prec_bits);
}

/**
 * @brief 计算半整数阶乘 (k − 0.5)! = Γ(k + 0.5)
 *
 * @details 公式:
 *   Γ(k + 0.5) = (2k−1)!! / 2^k × √π
 *
 *   其中 (2k−1)!! = 1 × 3 × 5 × ... × (2k−1) 为双阶乘
 *
 *   特殊情况:
 *     k = 0: Γ(0.5) = √π ≈ 1.7724538509...
 */
BigFloat BigFloat::half_factorial(int k, int prec_bits) {
  int guard = std::max(64, prec_bits / 10);
  int work_prec = prec_bits + guard;

  if (k == 0) {
    // Γ(0.5) = √π
    return sqrt(pi(work_prec));
  }

  // 计算 (2k−1)!! = 1 × 3 × 5 × ... × (2k−1)
  BigInt double_fact(1);
  for (int i = 1; i <= 2 * k - 1; i += 2) {
    double_fact = double_fact.mul_u32(static_cast<uint32_t>(i));
  }

  // result = (2k−1)!! / 2^k × √π
  BigFloat result = from_bigint(double_fact, work_prec);
  result = result.mul_pow2(-k); // 除以 2^k（高效指数调整）
  result *= sqrt(pi(work_prec));
  result.set_precision(prec_bits);
  return result;
}

// ============================================
// 十进制输出
// ============================================

/**
 * @brief 将 BigFloat 转换为十进制字符串
 * @param decimal_digits       有效数字位数
 * @param scientific_notation  是否使用科学记数法（默认 true）
 *
 * @details 算法:
 *
 *   1. **估算十进制指数**: e10 ≈ exponent × log₁₀(2) + bit_length × log₁₀(2)
 *
 *   2. **缩放**: 将值乘以 10^shift_10 使其成为整数:
 *      - 科学记数法: shift_10 = decimal_digits − e10
 *      - 定点表示: shift_10 = decimal_digits
 *      利用 10^n = 5^n × 2^n 分解来避免大数乘法的精度损失
 *
 *   3. **格式化**: 将整数形式的数字插入小数点和指数标记
 *
 *   注意: 由于 e10 是估算的，实际位数可能差 1，
 *         代码中有修正逻辑来处理这种情况。
 */
std::string BigFloat::to_decimal_string(int decimal_digits,
                                        bool scientific_notation) const {
  // 零值特殊处理
  if (is_zero()) {
    std::string s = "0.";
    s += std::string(decimal_digits, '0');
    if (scientific_notation)
      s += "e+00";
    return s;
  }

  std::string result;
  if (sign_)
    result = "-";

  // --- 步骤 1: 估算十进制指数 ---
  // e10 ≈ floor(value 的十进制量级)
  double log10_2 = 0.30102999566398119521;
  double log10_val = exponent_ * log10_2 + mantissa_.bit_length() * log10_2;
  int64_t e10 = static_cast<int64_t>(std::floor(log10_val));

  // 超大或超小数字强制使用科学记数法
  if (e10 >= 40 || e10 <= -10) {
    scientific_notation = true;
  }

  // --- 步骤 2: 确定缩放因子 ---
  int shift_10;
  if (scientific_notation) {
    // 科学记数法: 小数点前只有 1 位数字
    shift_10 = decimal_digits - static_cast<int>(e10);
  } else {
    // 定点数表示
    shift_10 = decimal_digits;
  }

  // --- 步骤 3: 缩放 ---
  // scaled = value × 10^shift_10
  // 利用 10^n = 5^n × 2^n 分解
  BigInt scaled = mantissa_;
  int64_t net_exp = exponent_;

  if (shift_10 > 0) {
    // 乘以 10^shift_10 = 5^shift_10 × 2^shift_10
    // 使用二进制快速幂计算 5^shift_10，替代逐步 mul_u32(5)
    BigInt five_power(1);
    BigInt base5(5);
    int p = shift_10;
    while (p > 0) {
      if (p & 1)
        five_power = five_power * base5;
      base5 = base5 * base5;
      p >>= 1;
    }
    scaled = scaled * five_power;
    net_exp += shift_10; // 记录 2^shift_10 的贡献到指数
  } else if (shift_10 < 0) {
    // 除以 10^|shift_10| = 除以 5^|shift_10| / 2^|shift_10|
    BigInt five_power(1);
    BigInt base5(5);
    int abs_shift = -shift_10;
    int p = abs_shift;
    while (p > 0) {
      if (p & 1)
        five_power = five_power * base5;
      base5 = base5 * base5;
      p >>= 1;
    }
    // 先乘以 2^extra_bits 保持精度，再做整数除法
    int extra_bits = abs_shift * 3 + 64;
    scaled = (scaled << extra_bits) / five_power;
    net_exp = net_exp - abs_shift - extra_bits;
  }

  // 应用剩余的二进制指数
  if (net_exp > 0) {
    scaled = scaled << static_cast<int>(net_exp);
  } else if (net_exp < 0) {
    scaled = scaled >> static_cast<int>(-net_exp);
  }

  std::string digits = scaled.to_decimal_string();

  // --- 步骤 4: 格式化输出 ---
  if (scientific_notation) {
    // 期望总位数为 decimal_digits + 1
    // 修正估算偏差（e10 可能差 1）
    while (static_cast<int>(digits.size()) < decimal_digits + 1) {
      digits += "0";
      e10--;
    }
    if (static_cast<int>(digits.size()) > decimal_digits + 1) {
      int extra = static_cast<int>(digits.size()) - (decimal_digits + 1);
      digits = digits.substr(0, decimal_digits + 1);
      e10 += extra;
    }

    // 格式: X.YYYYe+ZZ
    result += digits.substr(0, 1);
    result += ".";
    result += digits.substr(1);
    result += (e10 >= 0) ? "e+" : "e";
    result += std::to_string(e10);
  } else {
    // 定点数格式
    if (static_cast<int>(digits.size()) <= decimal_digits) {
      // 需要在前面补零
      std::string zeros(decimal_digits - digits.size() + 1, '0');
      digits = zeros + digits;
    }

    int int_len = static_cast<int>(digits.size()) - decimal_digits;
    result += digits.substr(0, int_len);
    result += ".";
    result += digits.substr(int_len);
  }

  return result;
}

/** @brief 流输出运算符（默认 20 位有效数字） */
std::ostream &operator<<(std::ostream &os, const BigFloat &f) {
  os << f.to_decimal_string(20);
  return os;
}
