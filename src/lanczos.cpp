#include "lanczos.h"
#include <iostream>
#include <cmath>
#include <cassert>

/**
 * Godfrey 方法计算 Lanczos 系数: P = D * B * C * F
 * 
 * 矩阵定义严格遵循 Boost 实现 (lanczos_generator.cpp):
 *
 * B[i][j]:
 *   - B[0][j] = 1 for all j
 *   - B[i][j] = 0 if j < i
 *   - B[i][j] = (-1)^(j-i) * C(i+j-1, j-i) if j >= i (i > 0)
 *
 * D[i][i] (对角矩阵):
 *   - D[0] = 1
 *   - D[1] = -1
 *   - D[i] = D[i-1] * 2*(2i-1) / (i-1) for i >= 2
 *
 * C[i][j] (Chebyshev 系数矩阵):
 *   - C[0][0] = 0.5
 *   - C[i][j] = 0 if j > i
 *   - C[i][j] = (-1)^(i-j) * sum(k=0..i, C(2i,2k)*C(k,k+j-i))
 *
 * F[i] (高精度浮点向量):
 *   - F[i] = (2i)! / i! * exp(i+g+0.5) / 2^(2i-1) / (g+i+0.5)^i / sqrt(g+i+0.5)
 */

// 计算组合数 C(n, k) 用 BigInt
static BigInt comb(int n, int k) {
    if (k < 0 || k > n) return BigInt(0);
    if (k == 0 || k == n) return BigInt(1);
    if (k > n - k) k = n - k;
    BigInt result(1);
    for (int i = 0; i < k; i++) {
        result = result.mul_u32(static_cast<uint32_t>(n - i));
        auto [q, r] = result.divmod_u32(static_cast<uint32_t>(i + 1));
        result = q;
    }
    return result;
}

// 带符号大整数
struct SignedBigInt {
    BigInt val;
    bool neg;
    SignedBigInt() : val(0), neg(false) {}
    SignedBigInt(const BigInt& v, bool n) : val(v), neg(n) {}
    
    SignedBigInt operator+(const SignedBigInt& other) const {
        if (neg == other.neg) {
            return {val + other.val, neg};
        }
        int cmp = val.compare(other.val);
        if (cmp == 0) return {BigInt(0), false};
        if (cmp > 0) return {val - other.val, neg};
        return {other.val - val, other.neg};
    }
    SignedBigInt operator-() const {
        if (val.is_zero()) return *this;
        return {val, !neg};
    }
};

std::vector<BigFloat> compute_lanczos_coefficients(int n, const std::string& g_str, int decimal_digits) {
    int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) + 64;
    int work_bits = bits + 256; // 增加保护位以提高精度

    BigFloat bf_g = BigFloat::from_string(g_str, work_bits);
    std::cout << "[lanczos] n=" << n << " g=" << g_str 
              << " decimal_digits=" << decimal_digits 
              << " binary_bits=" << bits << std::endl;

    // ========== 构造矩阵 B (n x n) ==========
    // B[0][j] = 1 for all j
    // B[i][j] = 0 if j < i (i > 0)
    // B[i][j] = (-1)^(j-i) * C(i+j-1, j-i) if j >= i (i > 0)
    std::cout << "[lanczos] Constructing matrix B..." << std::endl;
    std::vector<std::vector<SignedBigInt>> B(n, std::vector<SignedBigInt>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                B[i][j] = {BigInt(1), false};
            } else if (j >= i) {
                BigInt r = comb(i + j - 1, j - i);
                bool is_neg = ((j - i) % 2) != 0;
                B[i][j] = {r, is_neg};
            } else {
                B[i][j] = {BigInt(0), false};
            }
        }
    }

    // ========== 构造矩阵 C (n x n) ==========
    // C[0][0] = 0.5 → 用整数表示需要特殊处理
    // 但在 Boost 中它用浮点。对于整数矩阵我们可以乘以2来避免分数
    // 不过 Boost 实际使用浮点矩阵。我们也用浮点矩阵来简化
    // 
    // C[i][j] = 0 if j > i
    // C[0][0] = 0.5
    // C[i][j] = (-1)^(i-j) * sum(k=0..i, C(2i,2k)*C(k,k+j-i))
    std::cout << "[lanczos] Constructing Chebyshev matrix C..." << std::endl;
    
    // 注意 C[0][0] = 0.5, 其余都是整数
    // 因为最终要和浮点 F 相乘, 我们直接用 BigFloat 来存 C 矩阵
    // 但 Boost 实际上 C 和 B 和 D 都是在同一类型下运算的
    // 为了精确, 我们将 C 矩阵乘以 2 使得全部为整数, 最后在 F[0] 中补偿
    // 或者, 我们直接用 SignedBigInt 并记住 C[0][0] 有个额外因子 0.5
    
    // 更简单的方法: 用 BigFloat 矩阵来存 C
    std::vector<std::vector<BigFloat>> Cmat(n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0 && j == 0) {
                // C[0][0] = 0.5
                Cmat[0][0] = BigFloat(1, work_bits) / BigFloat(2, work_bits);
            } else if (j > i) {
                Cmat[i][j] = BigFloat(0, work_bits);
            } else {
                // C[i][j] = (-1)^(i-j) * sum(k=0..i, C(2i,2k)*C(k,k+j-i))
                BigInt sum(0);
                for (int k = 0; k <= i; k++) {
                    BigInt term = comb(2 * i, 2 * k) * comb(k, k + j - i);
                    sum = sum + term;
                }
                BigFloat bf_sum = BigFloat::from_bigint(sum, work_bits);
                if ((i - j) % 2 != 0) {
                    bf_sum = -bf_sum;
                }
                Cmat[i][j] = bf_sum;
            }
        }
    }

    // ========== 构造对角矩阵 D (n x n) ==========
    // D[0] = 1
    // D[1] = -1
    // D[i] = D[i-1] * 2*(2i-1) / (i-1) for i >= 2
    std::cout << "[lanczos] Constructing matrix D..." << std::endl;
    
    // D 的元素可以是负数 (D[1]=-1), 之后的元素通过递推可能也是负数或正数
    // D[1] = -1
    // D[2] = -1 * 2*3 / 1 = -6
    // D[3] = -6 * 2*5 / 2 = -30
    // D[4] = -30 * 2*7 / 3 = -140
    // 等等, 全部都是负数且绝对值递增
    // 用 BigFloat 来存储 D
    std::vector<BigFloat> D(n, BigFloat(0, work_bits));
    D[0] = BigFloat(1, work_bits);
    if (n > 1) {
        D[1] = BigFloat(-1, work_bits);
    }
    for (int i = 2; i < n; i++) {
        // D[i] = D[i-1] * 2*(2i-1) / (i-1)
        D[i] = D[i - 1] * BigFloat(2 * (2 * i - 1), work_bits);
        D[i] = D[i] / BigFloat(i - 1, work_bits);
    }

    // ========== 计算 M = D * B * C ==========
    // 由于 D 是对角的: M[i][j] = D[i] * (B * C)[i][j]
    // 先算 BC = B * C
    std::cout << "[lanczos] Computing M = D * B * C..." << std::endl;
    
    // BC[i][j] = sum_k B[i][k] * C[k][j]
    // B 的非零元素: B[0][k]=1 for all k; B[i][k] for k>=i
    std::vector<std::vector<BigFloat>> BC(n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            BigFloat sum(0, work_bits);
            for (int k = 0; k < n; k++) {
                if (B[i][k].val.is_zero()) continue;
                // B[i][k] * C[k][j]
                BigFloat b_ik = BigFloat::from_bigint(B[i][k].val, work_bits);
                if (B[i][k].neg) b_ik = -b_ik;
                sum += b_ik * Cmat[k][j];
            }
            BC[i][j] = sum;
        }
    }
    
    // M[i][j] = D[i] * BC[i][j]
    std::vector<std::vector<BigFloat>> M(n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = D[i] * BC[i][j];
        }
    }

    // ========== 构造向量 F (高精度浮点) ==========
    // F[i] = (2i)! / i! * exp(i+g+0.5) / 2^(2i-1) / (g+i+0.5)^i / sqrt(g+i+0.5)
    std::cout << "[lanczos] Constructing vector F..." << std::endl;
    
    BigFloat bf_half = BigFloat(1, work_bits) / BigFloat(2, work_bits);
    
    std::vector<BigFloat> F(n);
    for (int i = 0; i < n; i++) {
        BigFloat bi(i, work_bits);
        BigFloat base = bf_g + bi + bf_half;  // g + i + 0.5
        
        // (2i)! / i!
        BigFloat fact_ratio(1, work_bits);
        {
            BigInt f2i(1);
            for (int k = 2; k <= 2 * i; k++) {
                f2i = f2i.mul_u32(static_cast<uint32_t>(k));
            }
            BigInt fi(1);
            for (int k = 2; k <= i; k++) {
                fi = fi.mul_u32(static_cast<uint32_t>(k));
            }
            if (i == 0) {
                fact_ratio = BigFloat(1, work_bits);
            } else {
                BigFloat bf_f2i = BigFloat::from_bigint(f2i, work_bits);
                BigFloat bf_fi = BigFloat::from_bigint(fi, work_bits);
                fact_ratio = bf_f2i / bf_fi;
            }
        }
        
        // exp(i + g + 0.5)
        BigFloat exp_term = BigFloat::exp(base);
        
        // 1 / 2^(2i-1)  (注意 i=0 时 2^(-1) = 0.5)
        BigFloat pow2_term(1, work_bits);
        {
            int shift = 2 * i - 1;
            if (shift >= 0) {
                // 2^shift
                BigInt two_pow(1);
                for (int k = 0; k < shift; k++) {
                    two_pow = two_pow.mul_u32(2);
                }
                pow2_term = BigFloat::from_bigint(two_pow, work_bits);
            } else {
                // 2^(-1) = 0.5
                pow2_term = bf_half;  // shift = -1 for i=0
            }
        }
        
        // (g+i+0.5)^i
        BigFloat base_pow_i(1, work_bits);
        if (i > 0) {
            base_pow_i = BigFloat::pow(base, BigFloat(i, work_bits));
        }
        
        // sqrt(g+i+0.5)
        BigFloat base_sqrt = BigFloat::sqrt(base);
        
        // F[i] = fact_ratio * exp_term / pow2_term / base_pow_i / base_sqrt
        F[i] = fact_ratio * exp_term;
        if (i == 0) {
            // pow2_term = 0.5, so dividing by 0.5 = multiplying by 2
            F[i] = F[i] / pow2_term;
        } else {
            F[i] = F[i] / pow2_term;
        }
        F[i] = F[i] / base_pow_i;
        F[i] = F[i] / base_sqrt;
        
        std::cout << "  F[" << i << "] computed" << std::endl;
    }

    // ========== 计算 P = M * F ==========
    std::cout << "[lanczos] Computing P = M * F..." << std::endl;
    std::vector<BigFloat> P(n);
    for (int i = 0; i < n; i++) {
        BigFloat sum(0, work_bits);
        for (int j = 0; j < n; j++) {
            sum += M[i][j] * F[j];
        }
        P[i] = sum;
        P[i].set_precision(bits);
    }

    std::cout << "[lanczos] Coefficients computed successfully." << std::endl;
    return P;
}

// ============================================
// Gamma 函数计算
// ============================================

BigFloat lanczos_gamma(const BigFloat& z, const std::vector<BigFloat>& coeffs, const std::string& g_str, int decimal_digits) {
    int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) + 64;
    int work_bits = bits + 256; // 增加保护位以提高精度
    int n = static_cast<int>(coeffs.size());
    
    BigFloat one(1, work_bits);
    BigFloat half = one / BigFloat(2, work_bits);
    BigFloat bf_g = BigFloat::from_string(g_str, work_bits);
    
    BigFloat zz = z;
    zz.set_precision(work_bits);
    
    // 根据用户要求修改为新公式:
    // Gamma(z) = base^(z+0.5) * exp(-base) * S(z)
    // 其中 base = z + g + 0.5
    //      S(z) = p[0] + sum_{k=1}^{n-1} p[k] / (z + k)
    
    BigFloat base = zz + bf_g + half;  // z + g + 0.5
    
    // S = p[0] + sum_{k=1}^{n-1} p[k] / (z + k)
    BigFloat S = coeffs[0];
    S.set_precision(work_bits);
    for (int k = 1; k < n; k++) {
        BigFloat ck = coeffs[k];
        ck.set_precision(work_bits);
        BigFloat denom = zz + BigFloat(k, work_bits);
        S += ck / denom;
    }
    
    // Gamma(z+1) = base^(z+0.5) * e^(-base) * S
    BigFloat pow_term = BigFloat::pow(base, zz + half);
    BigFloat exp_term = BigFloat::exp(-base);
    
    // 用户提供的公式根据当前的 Lanczos 系数 (基于 Boost 方案) 
    // 实际上计算的是 Gamma(z+1)。所以计算 Gamma(z) 需要除以 z。
    // 即 Gamma(z) = Gamma(z+1) / z
    BigFloat result = (pow_term * exp_term * S) / zz;
    result.set_precision(bits);
    return result;
}
