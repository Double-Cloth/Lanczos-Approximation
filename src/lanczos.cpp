/**
 * @file lanczos.cpp
 * @brief Lanczos 近似系数计算与 Gamma 函数求值的实现
 *
 * @details
 * 本文件实现了 Godfrey 矩阵方法来计算 Lanczos 近似系数。
 *
 * Lanczos 近似公式（直接计算 Γ(z)）:
 *   base = z + g - 0.5
 *   S(z) = p₀ + Σ_{k=1}^{n-1} p_k / (z + k - 1)
 *   Γ(z) = base^(z-0.5) × e^(−base) × S(z)
 *
 * Godfrey 矩阵分解:
 *   P = D × B × C × F
 *
 *   - B[i][j]: 二项式系数矩阵
 *     B[0][j] = 1，B[i][j] = (−1)^(j−i) × C(i+j−1, j−i) (j ≥ i, i > 0)
 *
 *   - D[i][i]: 对角递推矩阵
 *     D[0]=1, D[1]=−1, D[i] = D[i−1] × 2(2i−1) / (i−1)
 *
 *   - C[i][j]: Chebyshev 系数矩阵
 *     C[0][0]=0.5, C[i][j] = (−1)^(i−j) × Σ C(2i,2k) × C(k,k+j−i)
 *
 *   - F[i]: 高精度浮点向量
 *     F[i] = (2i)!/i! × exp(i+0.5) / 2^(2i−1) / (g+i+0.5)^i / √(g+i+0.5)
 *
 * 参考: Boost C++ Libraries, lanczos_generator.cpp
 */

#include "lanczos.h"
#include <cassert>
#include <cmath>
#include <iostream>

// ============================================
// 辅助函数: 组合数
// ============================================

/**
 * @brief 计算组合数 C(n, k) = n! / (k! × (n−k)!)
 * @details 使用逐步乘除法避免中间结果过大:
 *          C(n,k) = n/1 × (n−1)/2 × ... × (n−k+1)/k
 *          利用 C(n,k) = C(n, n−k) 选择较小的 k
 */
static BigInt comb(int n, int k) {
  if (k < 0 || k > n)
    return BigInt(0);
  if (k == 0 || k == n)
    return BigInt(1);
  if (k > n - k)
    k = n - k; // 利用对称性选择较小的 k
  BigInt result(1);
  for (int i = 0; i < k; i++) {
    result = result.mul_u32(static_cast<uint32_t>(n - i));
    // 除以 (i+1), 由于 C(n,k) 始终为整数，此除法整除无余数
    auto [q, r] = result.divmod_u32(static_cast<uint32_t>(i + 1));
    result = q;
  }
  return result;
}

// ============================================
// 带符号大整数: 用于矩阵 B 的计算
// ============================================

/**
 * @brief 简单的有符号大整数包装
 * @details BigInt 本身只支持无符号运算，
 *          在构造矩阵 B 时需要处理 (−1)^n 的符号，
 *          所以用 val + neg 组合表示带符号值
 */
struct SignedBigInt {
  BigInt val; ///< 绝对值
  bool neg;   ///< true 表示负数

  SignedBigInt() : val(0), neg(false) {}

  SignedBigInt(const BigInt &v, bool n) : val(v), neg(n) {}

  /** @brief 带符号加法: 同号相加，异号做减法 */
  SignedBigInt operator+(const SignedBigInt &other) const {
    if (neg == other.neg) {
      return {val + other.val, neg}; // 同号: 绝对值相加
    }
    int cmp = val.compare(other.val);
    if (cmp == 0)
      return {BigInt(0), false}; // 完全抵消
    if (cmp > 0)
      return {val - other.val, neg};     // |a| > |b|, 符号跟 a
    return {other.val - val, other.neg}; // |b| > |a|, 符号跟 b
  }

  /** @brief 取反 */
  SignedBigInt operator-() const {
    if (val.is_zero())
      return *this;
    return {val, !neg};
  }
};

// ============================================
// 系数计算: P = D × B × C × F
// ============================================

/**
 * @brief 计算 Lanczos 近似系数
 *
 * @details 完整的计算流程:
 *
 *   1. 构造 B 矩阵 (n×n): 二项式系数矩阵
 *   2. 构造 C 矩阵 (n×n): Chebyshev 系数矩阵
 *   3. 构造 D 对角矩阵 (n×n): 递推生成
 *   4. 计算 M = D × B × C
 *   5. 构造 F 向量: 包含 exp、sqrt、阶乘的复合表达式
 *   6. 计算 P = M × F
 *
 *   每一步都使用扩展精度（work_bits = bits + 256）以减少累积误差
 */
std::vector<BigFloat> compute_lanczos_coefficients(int n,
                                                   const std::string &g_str,
                                                   int decimal_digits) {
  // 将十进制精度转换为二进制位数: bits ≈ digits × log2(10) + 64
  int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) + 64;
  int work_bits = bits + 256; // 增加 256 位保护位

  // 从字符串构造 g（保持全精度，不经过 double）
  BigFloat bf_g = BigFloat::from_string(g_str, work_bits);
  std::cout << "[lanczos] n=" << n << " g=" << g_str
            << " decimal_digits=" << decimal_digits << " binary_bits=" << bits
            << std::endl;

  // ========== 步骤 1: 构造矩阵 B (n × n) ==========
  // B[0][j] = 1                                 (第一行全为 1)
  // B[i][j] = 0                   if j < i      (下三角为 0，i > 0)
  // B[i][j] = (−1)^(j−i) × C(i+j−1, j−i)  if j ≥ i, i > 0
  std::cout << "[lanczos] Constructing matrix B..." << std::endl;
  std::vector<std::vector<SignedBigInt>> B(n, std::vector<SignedBigInt>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == 0) {
        B[i][j] = {BigInt(1), false}; // 第一行全 1
      } else if (j >= i) {
        BigInt r = comb(i + j - 1, j - i);
        bool is_neg = ((j - i) % 2) != 0; // (−1)^(j−i)
        B[i][j] = {r, is_neg};
      } else {
        B[i][j] = {BigInt(0), false}; // 下三角为 0
      }
    }
  }

  // ========== 步骤 2: 构造 Chebyshev 矩阵 C (n × n) ==========
  // C[0][0] = 0.5       (唯一的非整数元素)
  // C[i][j] = 0         if j > i (上三角为 0)
  // C[i][j] = (−1)^(i−j) × Σ_{k=0}^{i} C(2i,2k) × C(k, k+j−i)
  //
  // 由于 C[0][0] = 0.5 是分数，使用 BigFloat 矩阵存储
  std::cout << "[lanczos] Constructing Chebyshev matrix C..." << std::endl;
  std::vector<std::vector<BigFloat>> Cmat(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == 0 && j == 0) {
        // 特殊元素: C[0][0] = 0.5
        Cmat[0][0] = BigFloat(1, work_bits) / BigFloat(2, work_bits);
      } else if (j > i) {
        Cmat[i][j] = BigFloat(0, work_bits); // 上三角为 0
      } else {
        // 一般公式: C[i][j] = (−1)^(i−j) × Σ C(2i,2k) × C(k, k+j−i)
        BigInt sum(0);
        for (int k = 0; k <= i; k++) {
          BigInt term = comb(2 * i, 2 * k) * comb(k, k + j - i);
          sum = sum + term;
        }
        BigFloat bf_sum = BigFloat::from_bigint(sum, work_bits);
        if ((i - j) % 2 != 0) {
          bf_sum = -bf_sum; // 应用 (−1)^(i−j) 符号
        }
        Cmat[i][j] = bf_sum;
      }
    }
  }

  // ========== 步骤 3: 构造对角矩阵 D (n × n) ==========
  // D[0] = 1
  // D[1] = −1
  // D[i] = D[i−1] × 2(2i−1) / (i−1)   (i ≥ 2)
  //
  // D 序列的前几项: 1, −1, −6, −30, −140, ...
  // 绝对值快速增长
  std::cout << "[lanczos] Constructing matrix D..." << std::endl;
  std::vector<BigFloat> D(n, BigFloat(0, work_bits));
  D[0] = BigFloat(1, work_bits);
  if (n > 1) {
    D[1] = BigFloat(-1, work_bits);
  }
  for (int i = 2; i < n; i++) {
    // D[i] = D[i−1] × 2(2i−1) / (i−1)
    D[i] = D[i - 1] * BigFloat(2 * (2 * i - 1), work_bits);
    D[i] = D[i] / BigFloat(i - 1, work_bits);
  }

  // ========== 步骤 4: 计算 M = D × B × C ==========
  // 由于 D 是对角矩阵: M[i][j] = D[i] × (B × C)[i][j]
  // 先计算 BC = B × C（标准矩阵乘法）
  std::cout << "[lanczos] Computing M = D * B * C..." << std::endl;

  // BC[i][j] = Σ_k B[i][k] × C[k][j]
  std::vector<std::vector<BigFloat>> BC(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      BigFloat sum(0, work_bits);
      for (int k = 0; k < n; k++) {
        if (B[i][k].val.is_zero())
          continue; // 跳过零元素（稀疏优化）
        // 将 SignedBigInt 转为 BigFloat
        BigFloat b_ik = BigFloat::from_bigint(B[i][k].val, work_bits);
        if (B[i][k].neg)
          b_ik = -b_ik;
        sum += b_ik * Cmat[k][j];
      }
      BC[i][j] = sum;
    }
  }

  // M[i][j] = D[i] × BC[i][j]
  std::vector<std::vector<BigFloat>> M(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[i][j] = D[i] * BC[i][j];
    }
  }

  // ========== 步骤 5: 构造向量 F (高精度浮点) ==========
  // F[i] = (2i)!/i! × exp(i+0.5) / 2^(2i−1) / (g+i+0.5)^i / √(g+i+0.5)
  //
  // 各因子含义:
  //   (2i)!/i!        : 阶乘比（中心二项式系数的一半）
  //   exp(i+0.5)      : 新结构分离出的指数项，不含 g（吸收 e^-g 到系数中）
  //   1/2^(2i−1)      : 归一化因子
  //   1/(g+i+0.5)^i   : 基底的幂次
  //   1/√(g+i+0.5)    : 半整数幂的修正
  std::cout << "[lanczos] Constructing vector F..." << std::endl;

  BigFloat bf_half = BigFloat(1, work_bits) / BigFloat(2, work_bits);

  // 预计算 exp(1) 用于 F[i] 的 exp(i+0.5) 复用
  // power_i = i + 0.5, power_{i+1} = power_i + 1
  // 因此 exp(base_{i+1}) = exp(base_i) * exp(1)
  BigFloat exp_one = BigFloat::exp(BigFloat(1, work_bits));

  std::vector<BigFloat> F(n);
  BigFloat accumulated_exp; // 跨迭代累积的 exp 值
  for (int i = 0; i < n; i++) {
    BigFloat bi(i, work_bits);
    BigFloat base = bf_g + bi + bf_half; // base = g + i + 0.5

    // --- (2i)! / i! ---
    BigFloat fact_ratio(1, work_bits);
    {
      BigInt f2i(1); // 计算 (2i)!
      for (int k = 2; k <= 2 * i; k++) {
        f2i = f2i.mul_u32(static_cast<uint32_t>(k));
      }
      BigInt fi(1); // 计算 i!
      for (int k = 2; k <= i; k++) {
        fi = fi.mul_u32(static_cast<uint32_t>(k));
      }
      if (i == 0) {
        fact_ratio = BigFloat(1, work_bits); // 0!/0! = 1
      } else {
        BigFloat bf_f2i = BigFloat::from_bigint(f2i, work_bits);
        BigFloat bf_fi = BigFloat::from_bigint(fi, work_bits);
        fact_ratio = bf_f2i / bf_fi;
      }
    }

    // --- exp(i + 0.5) ---
    // 优化: 仅 i=0 时完整计算 exp(0.5)，后续乘以 exp(1)
    // 因为 exp(i+0.5) = exp(i-1+0.5) * exp(1)
    if (i == 0) {
      accumulated_exp = BigFloat::exp(bf_half); // exp(0.5)
    } else {
      accumulated_exp = accumulated_exp * exp_one;
    }
    BigFloat exp_term = accumulated_exp;

    // --- 1 / 2^(2i−1): 使用 mul_pow2 高效实现（O(1) 指数调整） ---
    // i=0: 2i-1 = -1, 即乘以 2 (mul_pow2(1))
    // i>0: 2i-1 > 0, 即除以 2^(2i-1) (mul_pow2(-(2i-1)))
    int pow2_shift = -(2 * i - 1);

    // --- (g+i+0.5)^i ---
    BigFloat base_pow_i(1, work_bits);
    if (i > 0) {
      base_pow_i = BigFloat::pow(base, BigFloat(i, work_bits));
    }

    // --- √(g+i+0.5) ---
    BigFloat base_sqrt = BigFloat::sqrt(base);

    // --- 组合: F[i] = fact_ratio × exp × 2^pow2_shift / base^i / √base ---
    F[i] = fact_ratio * exp_term;
    F[i] = F[i].mul_pow2(pow2_shift); // 高效乘除 2 的幂
    F[i] = F[i] / base_pow_i;
    F[i] = F[i] / base_sqrt;

    std::cout << "  F[" << i << "] computed" << std::endl;
  }

  // ========== 步骤 6: 计算 P = M × F (矩阵-向量乘法) ==========
  // P[i] = Σ_j M[i][j] × F[j]
  std::cout << "[lanczos] Computing P = M * F..." << std::endl;
  std::vector<BigFloat> P(n);
  for (int i = 0; i < n; i++) {
    BigFloat sum(0, work_bits);
    for (int j = 0; j < n; j++) {
      sum += M[i][j] * F[j];
    }
    P[i] = sum;
    P[i].set_precision(bits); // 截断到目标精度
  }

  std::cout << "[lanczos] Coefficients computed successfully." << std::endl;
  return P;
}

// ============================================
// Gamma 函数计算
// ============================================

/**
 * @brief 使用 Lanczos 系数直接计算 Γ(z)
 *
 * @details 计算公式:
 *
 *   base = z + g - 0.5
 *   S(z) = p₀ + Σ_{k=1}^{n-1} p_k / (z + k - 1)
 *       = p₀ + p₁/z + p₂/(z+1) + p₃/(z+2) + ...
 *   Γ(z) = (base/e)^(z-0.5) × S(z)
 *
 *   该形不仅通过吸收 e^-g 到系数中分离了部分指数项，
 *   在求值期更利用 base/e 作为底数将其降至一次求大数幂运算，消除了 exp 开销。
 */
BigFloat lanczos_gamma(const BigFloat &z, const std::vector<BigFloat> &coeffs,
                       const std::string &g_str, int decimal_digits) {
  // 将十进制精度转换为二进制位数
  int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) + 64;
  int work_bits = bits + 256; // 256 位保护位
  int n = static_cast<int>(coeffs.size());

  // 构造常用常量
  BigFloat one(1, work_bits);
  BigFloat half = one / BigFloat(2, work_bits); // 0.5
  BigFloat bf_g = BigFloat::from_string(g_str, work_bits);

  // 将 z 升精度到工作精度
  BigFloat zz = z;
  zz.set_precision(work_bits);

  // --- 反射公式支持 ---
  // 对于 z < 0.5，使用 Gamma(z) * Gamma(1-z) = pi / sin(pi * z)
  if (zz < half) {
    BigFloat one_minus_z = one - zz;
    // 递归调用计算大正实数的 Gamma
    BigFloat gamma_1_minus_z =
        lanczos_gamma(one_minus_z, coeffs, g_str, decimal_digits);

    // 计算 sin(pi * z)
    BigFloat pi_val = BigFloat::pi(work_bits);
    BigFloat sin_pi_z = BigFloat::sin(pi_val * zz);

    // 极点拦截：0 和负整数处的 Gamma 函数是无穷大
    if (zz.floor() == zz) {
      throw std::runtime_error("Gamma function pole at 0 or negative integer");
    }

    BigFloat result = pi_val / (sin_pi_z * gamma_1_minus_z);
    result.set_precision(bits);
    return result;
  }

  // 计算 base = z + g - 0.5
  BigFloat base = zz + bf_g - half;

  // 计算 S(z) = p[0] + p[1]/z + p[2]/(z+1) + p[3]/(z+2) + ...
  //            = p[0] + Σ_{k=1}^{n-1} p[k] / (z + k - 1)
  BigFloat S = coeffs[0];
  S.set_precision(work_bits);
  for (int k = 1; k < n; k++) {
    BigFloat ck = coeffs[k];
    ck.set_precision(work_bits);
    BigFloat denom = zz + BigFloat(k - 1, work_bits); // z + k - 1
    S += ck / denom;
  }

  // 计算 Γ(z) = (base/e)^(z-0.5) × S(z)
  BigFloat exp_one = BigFloat::exp(BigFloat(1, work_bits));
  BigFloat base_over_e = base / exp_one;
  BigFloat prefix = BigFloat::pow(base_over_e, zz - half); // (base/e)^(z-0.5)

  BigFloat result = prefix * S;
  result.set_precision(bits);
  return result;
}
