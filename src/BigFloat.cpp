#include "BigFloat.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <sstream>

// ============================================
// 构造函数
// ============================================

BigFloat::BigFloat()
    : sign_(false), mantissa_(0), exponent_(0), precision_bits_(256) {}

BigFloat::BigFloat(int val, int prec_bits)
    : sign_(val < 0), mantissa_(static_cast<uint64_t>(val < 0 ? -static_cast<int64_t>(val) : val)),
      exponent_(0), precision_bits_(prec_bits) {
    normalize();
}

BigFloat::BigFloat(int64_t val, int prec_bits)
    : sign_(val < 0), mantissa_(static_cast<uint64_t>(val < 0 ? -val : val)),
      exponent_(0), precision_bits_(prec_bits) {
    normalize();
}

BigFloat::BigFloat(double val, int prec_bits)
    : sign_(val < 0), exponent_(0), precision_bits_(prec_bits) {
    if (val == 0.0) {
        mantissa_ = BigInt(0);
        return;
    }
    double abs_val = std::abs(val);
    // 将 double 转化为整数尾数: 乘以 2^53 然后截断
    int exp2;
    double frac = std::frexp(abs_val, &exp2); // abs_val = frac * 2^exp2, 0.5 <= frac < 1
    // frac * 2^53 是整数
    uint64_t mant = static_cast<uint64_t>(frac * (1ULL << 53));
    mantissa_ = BigInt(mant);
    exponent_ = exp2 - 53;
    normalize();
}

BigFloat BigFloat::from_bigint(const BigInt& val, int prec_bits) {
    BigFloat result;
    result.sign_ = false;
    result.mantissa_ = val;
    result.exponent_ = 0;
    result.precision_bits_ = prec_bits;
    result.normalize();
    return result;
}

BigFloat BigFloat::from_fraction(const BigInt& num, const BigInt& den, int prec_bits) {
    // 计算 num / den 到 prec_bits 精度
    // result = (num << prec_bits) / den, exponent = -prec_bits
    BigInt shifted = num << (prec_bits + 32); // 多加 32 位保护位
    auto [q, r] = shifted.divmod(den);
    BigFloat result;
    result.sign_ = false;
    result.mantissa_ = q;
    result.exponent_ = -(prec_bits + 32);
    result.precision_bits_ = prec_bits;
    result.normalize();
    return result;
}

// ============================================
// 截断/规范化
// ============================================

void BigFloat::normalize() {
    if (mantissa_.is_zero()) {
        exponent_ = 0;
        sign_ = false;
        return;
    }
    int bl = mantissa_.bit_length();
    if (bl > precision_bits_) {
        int shift = bl - precision_bits_;
        // 四舍五入: 检查被丢弃的最高位 (round-to-nearest)
        bool round_up = (shift > 0) && mantissa_.get_bit(shift - 1);
        mantissa_ >>= shift;
        exponent_ += shift;
        if (round_up) {
            mantissa_ += BigInt(1);
            // 进位可能使位长增加一位
            if (mantissa_.bit_length() > precision_bits_) {
                mantissa_ >>= 1;
                exponent_++;
            }
        }
    } else if (bl < precision_bits_) {
        int shift = precision_bits_ - bl;
        mantissa_ <<= shift;
        exponent_ -= shift;
    }
}

// ============================================
// 符号
// ============================================

BigFloat BigFloat::operator-() const {
    BigFloat result = *this;
    if (!result.is_zero()) result.sign_ = !result.sign_;
    return result;
}

BigFloat BigFloat::abs() const {
    BigFloat result = *this;
    result.sign_ = false;
    return result;
}

// ============================================
// 对齐指数
// ============================================

void BigFloat::align(const BigFloat& a, const BigFloat& b,
                     BigInt& a_mant, BigInt& b_mant, int64_t& common_exp) {
    // 提早处理 0 的情况，防止无意义的位移
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
        int max_prec = std::max(a.precision_bits_, b.precision_bits_) + 32; // 32位保护位
        
        if (diff > max_prec) {
            // b 的指数远大于 a，说明 b 的绝对值远大于 a
            // a 完全可以忽略不计（相当于 a 为 0）
            a_mant = BigInt(0);
            b_mant = b.mantissa_;
            common_exp = b.exponent_; // 以大数 b 的指数为准
        } else {
            // 正常对齐：将大数 b 的尾数左移，指数统一为较小的 a.exponent_
            a_mant = a.mantissa_;
            b_mant = b.mantissa_ << static_cast<int>(diff);
            common_exp = a.exponent_;
        }
    } else {
        int64_t diff = a.exponent_ - b.exponent_;
        int max_prec = std::max(a.precision_bits_, b.precision_bits_) + 32;
        
        if (diff > max_prec) {
            // a 的指数远大于 b，说明 a 的绝对值远大于 b
            // b 完全可以忽略不计
            a_mant = a.mantissa_;
            b_mant = BigInt(0);
            common_exp = a.exponent_; // 以大数 a 的指数为准
        } else {
            // 正常对齐：将大数 a 的尾数左移，指数统一为较小的 b.exponent_
            b_mant = b.mantissa_;
            a_mant = a.mantissa_ << static_cast<int>(diff);
            common_exp = b.exponent_;
        }
    }
}

// --- 从十进制字符串构造 ---
BigFloat BigFloat::from_string(const std::string& str, int prec_bits) {
    if (str.empty()) return BigFloat(0, prec_bits);
    int i = 0;
    bool sign = false;
    if (str[i] == '-') { sign = true; i++; }
    else if (str[i] == '+') { i++; }

    std::string int_part = "";
    std::string frac_part = "";
    int exp_val = 0;

    // parse int part
    while (i < str.size() && std::isdigit(str[i])) {
        int_part += str[i++];
    }
    // parse frac part
    if (i < str.size() && str[i] == '.') {
        i++;
        while (i < str.size() && std::isdigit(str[i])) {
            frac_part += str[i++];
        }
    }
    // parse exponent
    if (i < str.size() && (str[i] == 'e' || str[i] == 'E')) {
        i++;
        std::string exp_str = "";
        while (i < str.size()) {
            exp_str += str[i++];
        }
        if (!exp_str.empty()) {
            exp_val = std::stoi(exp_str);
        }
    }

    std::string combined = int_part + frac_part;
    
    // remove leading zeros
    int start = 0;
    while (start < static_cast<int>(combined.size()) - 1 && combined[start] == '0') start++;
    combined = combined.substr(start);

    if (combined.empty() || combined == "0") return BigFloat(0, prec_bits);

    int dec_exp = exp_val - static_cast<int>(frac_part.size());
    BigInt num(combined);

    int work_prec = prec_bits + 64;
    BigFloat result = from_bigint(num, work_prec);
    
    auto power5 = [](int p) {
        BigInt res(1);
        BigInt base(5);
        while (p > 0) {
            if (p & 1) res = res * base;
            base = base * base;
            p >>= 1;
        }
        return res;
    };
    
    if (dec_exp > 0) {
        BigInt five_pow = power5(dec_exp);
        result.mantissa_ = result.mantissa_ * five_pow;
        result.exponent_ += dec_exp;
    } else if (dec_exp < 0) {
        int abs_exp = -dec_exp;
        BigInt five_pow = power5(abs_exp);
        
        int bl_num = result.mantissa_.bit_length();
        int bl_five = five_pow.bit_length();
        int K = work_prec + bl_five - bl_num + 32;
        if (K < 0) K = 0;

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

int BigFloat::compare(const BigFloat& other) const {
    if (is_zero() && other.is_zero()) return 0;
    if (is_zero()) return other.sign_ ? 1 : -1;
    if (other.is_zero()) return sign_ ? -1 : 1;
    if (sign_ != other.sign_) return sign_ ? -1 : 1;

    // 同号: 比较绝对值
    BigInt a_mant, b_mant;
    int64_t common_exp;
    align(*this, other, a_mant, b_mant, common_exp);
    int cmp = a_mant.compare(b_mant);
    return sign_ ? -cmp : cmp;
}

// ============================================
// 加减法
// ============================================

BigFloat BigFloat::operator+(const BigFloat& other) const {
    if (is_zero()) { BigFloat r = other; r.precision_bits_ = std::max(precision_bits_, other.precision_bits_); return r; }
    if (other.is_zero()) { BigFloat r = *this; r.precision_bits_ = std::max(precision_bits_, other.precision_bits_); return r; }

    BigFloat result;
    result.precision_bits_ = std::max(precision_bits_, other.precision_bits_);

    BigInt a_mant, b_mant;
    int64_t common_exp;
    align(*this, other, a_mant, b_mant, common_exp);

    if (sign_ == other.sign_) {
        result.sign_ = sign_;
        result.mantissa_ = a_mant + b_mant;
    } else {
        int cmp = a_mant.compare(b_mant);
        if (cmp == 0) return BigFloat(0, result.precision_bits_);
        if (cmp > 0) {
            result.sign_ = sign_;
            result.mantissa_ = a_mant - b_mant;
        } else {
            result.sign_ = other.sign_;
            result.mantissa_ = b_mant - a_mant;
        }
    }
    result.exponent_ = common_exp;
    result.normalize();
    return result;
}

BigFloat BigFloat::operator-(const BigFloat& other) const {
    BigFloat neg = -other;
    return *this + neg;
}

BigFloat& BigFloat::operator+=(const BigFloat& other) { *this = *this + other; return *this; }
BigFloat& BigFloat::operator-=(const BigFloat& other) { *this = *this - other; return *this; }

// ============================================
// 乘法
// ============================================

BigFloat BigFloat::operator*(const BigFloat& other) const {
    if (is_zero() || other.is_zero()) {
        return BigFloat(0, std::max(precision_bits_, other.precision_bits_));
    }
    BigFloat result;
    result.precision_bits_ = std::max(precision_bits_, other.precision_bits_);
    result.sign_ = sign_ != other.sign_;
    result.mantissa_ = mantissa_ * other.mantissa_;
    result.exponent_ = exponent_ + other.exponent_;
    result.normalize();
    return result;
}

BigFloat& BigFloat::operator*=(const BigFloat& other) { *this = *this * other; return *this; }

// ============================================
// 除法: a / b = a * (1/b), 使用牛顿迭代求倒数
// ============================================

BigFloat BigFloat::operator/(const BigFloat& other) const {
    if (other.is_zero()) throw std::runtime_error("BigFloat division by zero");
    if (is_zero()) return BigFloat(0, std::max(precision_bits_, other.precision_bits_));

    int prec = std::max(precision_bits_, other.precision_bits_);
    
    // 直接用 BigInt 整数除法
    // a/b = (a.mant * 2^a.exp) / (b.mant * 2^b.exp)
    //     = (a.mant / b.mant) * 2^(a.exp - b.exp)
    // 为了保证精度，先将 a.mant 左移足够多的位
    int shift = prec + 32;
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

BigFloat& BigFloat::operator/=(const BigFloat& other) { *this = *this / other; return *this; }

BigFloat BigFloat::mul_pow2(int n) const {
    BigFloat result = *this;
    result.exponent_ += n;
    return result;
}

// ============================================
// 数学常数与函数
// ============================================

// --- pi: Machin 公式 pi/4 = 4*arctan(1/5) - arctan(1/239) ---
// arctan(1/x) = 1/x - 1/(3*x^3) + 1/(5*x^5) - ...

static BigFloat arctan_recip(int x, int prec_bits) {
    // 计算 arctan(1/x) 使用级数
    BigFloat result(0, prec_bits);
    BigFloat x_bf(x, prec_bits);
    BigFloat x2 = x_bf * x_bf;
    BigFloat power = BigFloat(1, prec_bits) / x_bf; // 1/x
    
    for (int k = 0; ; k++) {
        BigFloat denom(2 * k + 1, prec_bits);
        BigFloat term = power / denom;
        if (term.mantissa().bit_length() + term.exponent() < -static_cast<int64_t>(prec_bits + 16)) {
            break; // 项已经足够小
        }
        if (k % 2 == 0)
            result += term;
        else
            result -= term;
        power /= x2;
    }
    return result;
}

BigFloat BigFloat::pi(int prec_bits) {
    // Machin: pi/4 = 4*arctan(1/5) - arctan(1/239)
    BigFloat a = arctan_recip(5, prec_bits + 32);
    BigFloat b = arctan_recip(239, prec_bits + 32);
    BigFloat four(4, prec_bits + 32);
    BigFloat result = (four * a - b) * four;
    result.set_precision(prec_bits);
    return result;
}

BigFloat BigFloat::e(int prec_bits) {
    return exp(BigFloat(1, prec_bits));
}

// --- sqrt: 牛顿迭代 x_{n+1} = (x_n + a/x_n) / 2 ---
BigFloat BigFloat::sqrt(const BigFloat& x) {
    if (x.is_zero()) return BigFloat(0, x.precision());
    if (x.is_negative()) throw std::runtime_error("sqrt of negative number");
    
    int prec = x.precision();
    int work_prec = prec + 64; // 增加保护位
    
    // 初始估计: 使用 double
    double approx = 0.5; // 粗略
    int bl = x.mantissa().bit_length();
    int total_exp = bl + static_cast<int>(x.exponent());
    approx = std::ldexp(1.0, total_exp / 2);
    
    BigFloat guess(approx, work_prec);
    BigFloat half(0, work_prec);
    {
        BigFloat one(1, work_prec);
        BigFloat two(2, work_prec);
        half = one / two;
    }
    BigFloat val = x;
    val.set_precision(work_prec);
    
    // 牛顿迭代
    int iters = 0;
    int max_iters = prec; // 最坏情况
    while (iters < max_iters) {
        BigFloat next = (guess + val / guess) * half;
        // 检查收敛: |next - guess| 足够小
        BigFloat diff = (next - guess).abs();
        guess = next;
        iters++;
        if (diff.is_zero()) break;
        // diff 的有效位比 guess 小至少 prec 位则收敛
        int diff_bl = diff.mantissa().bit_length() + static_cast<int>(diff.exponent());
        int guess_bl = guess.mantissa().bit_length() + static_cast<int>(guess.exponent());
        if (guess_bl - diff_bl > prec + 8) break;
    }
    guess.set_precision(prec);
    return guess;
}

// --- exp: 泰勒级数 e^x = sum(x^k / k!) ---
// 对大参数，先做 range reduction: e^x = (e^(x/2^r))^(2^r)
BigFloat BigFloat::exp(const BigFloat& x) {
    if (x.is_zero()) return BigFloat(1, x.precision());
    
    int prec = x.precision();
    int work_prec = prec + 128; // 增加保护位以抵消平方累积误差
    
    BigFloat val = x;
    val.set_precision(work_prec);
    bool negate = val.is_negative();
    if (negate) val = val.abs();
    
    // Range reduction: 将 val 除以 2^r 使得 val 足够小
    int r = 0;
    {
        int bl = val.mantissa().bit_length();
        int total = bl + static_cast<int>(val.exponent());
        if (total > 0) {
            r = total + 8; // 适度缩小，减少平方次数以降低误差累积
        } else {
            r = std::max(0, 8 + total); // 避免不必要的过度缩减
        }
    }
    BigFloat reduced = val.mul_pow2(-r);
    
    // 泰勒级数: e^reduced = 1 + reduced + reduced^2/2! + ...
    BigFloat sum(1, work_prec);
    BigFloat term(1, work_prec);
    for (int k = 1; ; k++) {
        term *= reduced;
        term /= BigFloat(k, work_prec);
        sum += term;
        // 检查 term 是否足够小
        if (term.is_zero()) break;
        int term_bl = term.mantissa().bit_length() + static_cast<int>(term.exponent());
        int sum_bl = sum.mantissa().bit_length() + static_cast<int>(sum.exponent());
        if (sum_bl - term_bl > work_prec) break;
    }
    
    // 反向 squaring: result = sum^(2^r)
    for (int i = 0; i < r; i++) {
        sum = sum * sum;
    }
    
    if (negate) {
        sum = BigFloat(1, work_prec) / sum;
    }
    
    sum.set_precision(prec);
    return sum;
}

// --- ln: AGM method 不够简单, 用 Taylor + range reduction ---
// ln(x) = ln(m * 2^e) = ln(m) + e * ln(2)
// 对 m in [0.5, 1), 用 ln(1+t) 的级数, t = m - 1
// 实际做法: 利用 ln(x) = 2 * atanh( (x-1)/(x+1) )
// atanh(t) = t + t^3/3 + t^5/5 + ...
BigFloat BigFloat::ln(const BigFloat& x) {
    if (x.is_zero() || x.is_negative()) throw std::runtime_error("ln of non-positive number");
    
    int prec = x.precision();
    int work_prec = prec + 128; // 增加保护位

    // 分解: x = mant * 2^exp
    // 提取使 mant 在 [1, 2) 附近
    BigFloat val = x;
    val.set_precision(work_prec);
    int bl = val.mantissa().bit_length();
    int64_t total_exp_bits = static_cast<int64_t>(bl) + val.exponent();
    
    // 构造 m: 将 val 缩放到 [1, 2)
    // m = val * 2^(-total_exp_bits + 1)
    BigFloat m = val.mul_pow2(static_cast<int>(-total_exp_bits + 1));
    int64_t e_val = total_exp_bits - 1; // val ≈ m * 2^e_val, m ∈ [1, 2)

    // ln(val) = ln(m) + e_val * ln(2)  
    // ln(m) = 2 * atanh( (m-1)/(m+1) ), 其中 m ∈ [1,2) 所以 (m-1)/(m+1) ∈ [0, 1/3)
    
    BigFloat one(1, work_prec);
    BigFloat t = (m - one) / (m + one); // t ∈ [0, 1/3)
    BigFloat t2 = t * t;
    
    // atanh(t) = t + t^3/3 + t^5/5 + ...
    BigFloat sum = t;
    BigFloat power = t;
    for (int k = 1; ; k++) {
        power *= t2;
        BigFloat term = power / BigFloat(2 * k + 1, work_prec);
        sum += term;
        if (term.is_zero()) break;
        int term_bl = term.mantissa().bit_length() + static_cast<int>(term.exponent());
        int sum_bl = sum.mantissa().bit_length() + static_cast<int>(sum.exponent());
        if (sum_bl - term_bl > work_prec) break;
    }
    
    BigFloat ln_m = sum * BigFloat(2, work_prec);
    
    // 计算 ln(2)
    BigFloat ln2(0, work_prec);
    if (e_val != 0) {
        // ln(2) = 2 * atanh(1/3) = 2 * (1/3 + 1/3^3/3 + 1/3^5/5 + ...)
        BigFloat third = one / BigFloat(3, work_prec);
        BigFloat third2 = third * third;
        BigFloat ln2_sum = third;
        BigFloat ln2_power = third;
        for (int k = 1; ; k++) {
            ln2_power *= third2;
            BigFloat term = ln2_power / BigFloat(2 * k + 1, work_prec);
            ln2_sum += term;
            if (term.is_zero()) break;
            int term_bl = term.mantissa().bit_length() + static_cast<int>(term.exponent());
            int sum_bl = ln2_sum.mantissa().bit_length() + static_cast<int>(ln2_sum.exponent());
            if (sum_bl - term_bl > work_prec) break;
        }
        ln2 = ln2_sum * BigFloat(2, work_prec);
    }
    
    BigFloat result = ln_m + BigFloat(static_cast<int64_t>(e_val), work_prec) * ln2;
    result.set_precision(prec);
    return result;
}

// --- pow(base, exp) = exp(exp * ln(base)) ---
// 对整数指数使用二进制快速幂以避免 ln 引入的精度损失
BigFloat BigFloat::pow(const BigFloat& base, const BigFloat& exponent) {
    if (exponent.is_zero()) return BigFloat(1, base.precision());
    if (base.is_zero()) return BigFloat(0, base.precision());
    
    int prec = std::max(base.precision(), exponent.precision());
    int work_prec = prec + 128;
    
    // 检查指数是否为非负整数
    // 条件: exponent >= 0 且 exponent 为 0 或 exponent_ >= 0 (没有小数位)
    if (!exponent.is_negative() && exponent.exponent() >= 0) {
        // 指数是非负整数，提取整数值
        BigInt mant = exponent.mantissa();
        int64_t shift = exponent.exponent();
        if (shift > 0 && shift <= 64) {
            mant <<= static_cast<int>(shift);
        }
        int64_t total_bits = static_cast<int64_t>(mant.bit_length());
        if (total_bits <= 31) {
            // 安全地转为整数
            std::string n_str = mant.to_decimal_string();
            int64_t n_val = std::stoll(n_str);
            if (n_val >= 0 && n_val <= 1000) {
                // 使用二进制快速幂
                BigFloat result(1, work_prec);
                BigFloat b = base;
                b.set_precision(work_prec);
                int64_t p = n_val;
                while (p > 0) {
                    if (p & 1) result *= b;
                    b *= b;
                    p >>= 1;
                }
                result.set_precision(prec);
                return result;
            }
        }
    }
    
    // 通用路径: exp(exp * ln(base))
    BigFloat ln_base = ln(base.abs());
    ln_base.set_precision(work_prec);
    BigFloat exp_val = exponent;
    exp_val.set_precision(work_prec);
    
    return exp(exp_val * ln_base);
}

// --- factorial: n! ---
BigFloat BigFloat::factorial(int n, int prec_bits) {
    BigInt result(1);
    for (int i = 2; i <= n; i++) {
        result = result.mul_u32(static_cast<uint32_t>(i));
    }
    return from_bigint(result, prec_bits);
}

// --- half_factorial: (k - 0.5)! = Gamma(k + 0.5) ---
// (k - 0.5)! = (2k)! * sqrt(pi) / (4^k * k!)
// 即 Gamma(k + 0.5) = (2k-1)!! / 2^k * sqrt(pi)
// (2k-1)!! = 1 * 3 * 5 * ... * (2k-1)
BigFloat BigFloat::half_factorial(int k, int prec_bits) {
    int work_prec = prec_bits + 64;
    
    if (k == 0) {
        // (-0.5)! = Gamma(0.5) = sqrt(pi)
        return sqrt(pi(work_prec));
    }
    
    // (k - 0.5)! = (2k-1)!! / 2^k * sqrt(pi)
    BigInt double_fact(1);
    for (int i = 1; i <= 2 * k - 1; i += 2) {
        double_fact = double_fact.mul_u32(static_cast<uint32_t>(i));
    }
    
    BigFloat result = from_bigint(double_fact, work_prec);
    result = result.mul_pow2(-k); // / 2^k
    result *= sqrt(pi(work_prec));
    result.set_precision(prec_bits);
    return result;
}

// ============================================
// 十进制输出
// ============================================

std::string BigFloat::to_decimal_string(int decimal_digits, bool scientific_notation) const {
    if (is_zero()) {
        std::string s = "0.";
        s += std::string(decimal_digits, '0');
        if (scientific_notation) s += "e+00";
        return s;
    }
    
    std::string result;
    if (sign_) result = "-";
    
    // 将值转化为: value = mantissa * 2^exponent
    // 策略: 计算需要乘以的 10 的幂 target_exp，以提取所需的有效数字
    // 我们想要提取大约 N 个十进制有效数字，使得结果是一个大整数
    
    // 粗略估计 10 的指数: e10 = floor( exponent * log10(2) + log10(mantissa) )
    double log10_2 = 0.30102999566398119521;
    double log10_val = exponent_ * log10_2 + mantissa_.bit_length() * log10_2;
    int64_t e10 = static_cast<int64_t>(std::floor(log10_val));
    
    // 如果数字非常大或非常小，强制启用科学记数法
    if (e10 >= 40 || e10 <= -10) {
        scientific_notation = true;
    }
    
    int shift_10;
    if (scientific_notation) {
        // 科学记数法: 我们需要小数前面只有1位
        // 提取的有效数字总数为 decimal_digits + 1
        shift_10 = decimal_digits - static_cast<int>(e10);
    } else {
        // 定点数表示
        shift_10 = decimal_digits;
    }
    
    // 计算 scaled = value * 10^shift_10
    BigInt scaled = mantissa_;
    int64_t net_exp = exponent_;
    
    if (shift_10 > 0) {
        // * 10^shift_10 = * 5^shift_10 * 2^shift_10
        BigInt five_power(1);
        for (int i = 0; i < shift_10; i++) {
            five_power = five_power.mul_u32(5);
        }
        scaled = scaled * five_power;
        net_exp += shift_10;
    } else if (shift_10 < 0) {
        // / 10^(-shift_10) = / 5^(-shift_10) / 2^(-shift_10)
        // 为了保持精度，我们最好先扩大然后再除以 5^k
        // 由于这仅仅是输出格式化，我们可以多乘以一个足够大的 2 的幂，做整数除法
        BigInt five_power(1);
        int abs_shift = -shift_10;
        for (int i = 0; i < abs_shift; i++) {
            five_power = five_power.mul_u32(5);
        }
        
        // 扩大 2^k 以保持精度，最后再右移回来
        int extra_bits = abs_shift * 3 + 64; 
        scaled = (scaled << extra_bits) / five_power;
        net_exp = net_exp - abs_shift - extra_bits;
    }
    
    if (net_exp > 0) {
        scaled = scaled << static_cast<int>(net_exp);
    } else if (net_exp < 0) {
        scaled = scaled >> static_cast<int>(-net_exp);
    }
    
    std::string digits = scaled.to_decimal_string();
    
    if (scientific_notation) {
        // 科学记数法: 期望位数为 decimal_digits + 1
        // 注意：由于前面的指数估计可能差1，digits 的长度可能比 decimal_digits + 1 大1或小1
        while (static_cast<int>(digits.size()) < decimal_digits + 1) {
            digits += "0";
            e10--; // 实际数字比估计的更小
        }
        if (static_cast<int>(digits.size()) > decimal_digits + 1) {
            int extra = static_cast<int>(digits.size()) - (decimal_digits + 1);
            digits = digits.substr(0, decimal_digits + 1);
            e10 += extra; // 实际数字比估计的更大
        }
        
        result += digits.substr(0, 1);
        result += ".";
        result += digits.substr(1);
        result += (e10 >= 0) ? "e+" : "e";
        result += std::to_string(e10);
    } else {
        // 定点数
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

std::ostream& operator<<(std::ostream& os, const BigFloat& f) {
    os << f.to_decimal_string(20);
    return os;
}
