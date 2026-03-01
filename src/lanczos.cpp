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
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================
// 辅助函数: 组合数
// ============================================

/**
 * @brief 计算组合数 C(n, k) = n! / (k! × (n−k)!)
 * @details 使用逐步乘除法避免中间结果过大:
 *          C(n,k) = n/1 × (n−1)/2 × ... × (n−k+1)/k
 *          利用 C(n,k) = C(n, n−k) 选择较小的 k
 *          带有全局记忆化缓存提升计算速度
 */
/**
 * @brief 全局组合数内存池（支持单线程填充和跨线程只读安全读取）
 */
static std::vector<std::vector<BigInt>> g_comb_cache;

/**
 * @brief 计算组合数 C(n, k) = n! / (k! × (n−k)!)
 * @details 带有全局记忆化缓存提升计算速度
 */
static BigInt comb(int n, int k) {
  if (k < 0 || k > n)
    return BigInt(0);
  if (k == 0 || k == n)
    return BigInt(1);
  if (k > n - k)
    k = n - k;

  if (n < (int)g_comb_cache.size() && k < (int)g_comb_cache[n].size() &&
      !g_comb_cache[n][k].is_zero()) {
    return g_comb_cache[n][k];
  }

  BigInt result(1);
  for (int i = 0; i < k; i++) {
    result = result.mul_u32(static_cast<uint32_t>(n - i));
    auto [q, r] = result.divmod_u32(static_cast<uint32_t>(i + 1));
    result = q;
  }

  if (n >= (int)g_comb_cache.size()) {
    g_comb_cache.resize(n + 1);
  }
  if (k >= (int)g_comb_cache[n].size()) {
    g_comb_cache[n].resize(k + 1, BigInt(0));
  }
  g_comb_cache[n][k] = result;

  return result;
}

/**
 * @brief 预填充组合数缓存，使得后续多线程读取可以无锁安全执行
 *        包含了二进制文件的序列化落盘与解吸附。
 * @param max_val 最大的组合数 n 值
 */
static void precompute_comb_cache(int max_val) {
  const std::string cache_file = "comb_cache.dat";

  // 尝试从磁盘反序列化
  std::ifstream is(cache_file, std::ios::binary);
  if (is.is_open()) {
    int stored_max;
    if (is.read(reinterpret_cast<char *>(&stored_max), sizeof(stored_max))) {
      if (stored_max >= max_val) {
        std::cout << "[lanczos] Loading combination cache from disk: "
                  << cache_file << " ... ";
        g_comb_cache.resize(stored_max + 1);
        for (int i = 0; i <= stored_max; i++) {
          int size_i;
          is.read(reinterpret_cast<char *>(&size_i), sizeof(size_i));
          g_comb_cache[i].resize(size_i);
          for (int j = 0; j < size_i; j++) {
            uint32_t vec_size;
            is.read(reinterpret_cast<char *>(&vec_size), sizeof(vec_size));
            if (vec_size > 0) {
              // Hack: Const Cast to reconstruct BigInt internal
              auto &digits = const_cast<std::vector<uint32_t> &>(
                  g_comb_cache[i][j].digits());
              digits.resize(vec_size);
              is.read(reinterpret_cast<char *>(digits.data()),
                      vec_size * sizeof(uint32_t));
            } else {
              g_comb_cache[i][j] = BigInt(0);
            }
          }
        }
        std::cout << "Done." << std::endl;
        return;
      } else {
        std::cout << "[lanczos] Stale cache found (stored max " << stored_max
                  << " < req " << max_val << "). Rebuilding..." << std::endl;
      }
    }
  }

  auto start_time = std::chrono::high_resolution_clock::now();
  auto last_time = start_time;
  int last_count = 0;
  // comb generation has internal nested loop up to k, so it's O(n^2), integrate
  // gives power = 1 + 1 = 2
  const int complexity = 1;

  for (int i = 0; i <= max_val; i++) {
    for (int j = 0; j <= i; j++) {
      comb(i, j);
    }

    int current = i + 1;
    int n_total = max_val + 1;
    int progress = (current * 100) / n_total;
    int bar_pos = (50 * current) / n_total;

    static double cached_remaining_ms = 0.0;

    // Check if enough time has passed to avoid micro-jitter
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff_from_last = now - last_time;

    if (diff_from_last.count() >= 500.0 || current == n_total) {
      std::chrono::duration<double, std::milli> elapsed = now - start_time;
      if (elapsed.count() >= 100.0 || current == n_total) {
        // Growth-based ETA Integration: Real work done proportion (O(N^{p+1}))
        double work_ratio =
            std::pow(static_cast<double>(current) / n_total, complexity + 1.0);
        double projected_total_ms = elapsed.count() / work_ratio;
        cached_remaining_ms = projected_total_ms - elapsed.count();
        // Ensure monotonic behavior explicitly on screen
        if (cached_remaining_ms < 0)
          cached_remaining_ms = 0.0;
      }
      last_time = now;
      last_count = current;
    }

    std::cout << "\r[lanczos] cache C(n,k) [";
    for (int p = 0; p < 50; ++p)
      std::cout << (p < bar_pos ? "=" : (p == bar_pos ? ">" : " "));
    std::cout << "] " << progress << "% (" << current << "/" << n_total << ") ";

    // Only display ETA if we've gathered at least a small chunk of elapsed time
    // to be meaningful
    auto total_elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time)
            .count();
    if (current < n_total && total_elapsed >= 500) {
      std::cout << "ETA: " << std::fixed << std::setprecision(1)
                << (cached_remaining_ms / 1000.0) << "s   ";
    }
    std::cout << std::flush;
  }
  std::cout << "\n[lanczos] Serializing cache to disk... ";

  // 保存到本地磁盘
  std::ofstream os(cache_file, std::ios::binary);
  if (os.is_open()) {
    os.write(reinterpret_cast<const char *>(&max_val), sizeof(max_val));
    for (int i = 0; i <= max_val; i++) {
      int size_i = g_comb_cache[i].size();
      os.write(reinterpret_cast<const char *>(&size_i), sizeof(size_i));
      for (int j = 0; j < size_i; j++) {
        const auto &digits = g_comb_cache[i][j].digits();
        uint32_t vec_size = digits.size();
        if (g_comb_cache[i][j].is_zero())
          vec_size = 0;
        os.write(reinterpret_cast<const char *>(&vec_size), sizeof(vec_size));
        if (vec_size > 0) {
          os.write(reinterpret_cast<const char *>(digits.data()),
                   vec_size * sizeof(uint32_t));
        }
      }
    }
  }
  std::cout << "Done." << std::endl;
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
  // 将十进制精度转换为二进制位数
  int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) +
             std::max(64, decimal_digits / 10);
  // 由于项数增加会导致矩阵运算中的阶乘级振荡（灾难性相消）
  // 保护位需随 n 动态极大增加: 至少需要 2n*log2(2n) 位的额外绝对精度
  int guard_bits = 256 + static_cast<int>(2.0 * 2 * n * std::log2(2 * n));
  int work_bits = bits + guard_bits;

  // 从字符串构造 g（保持全精度，不经过 double）
  BigFloat bf_g = BigFloat::from_string(g_str, work_bits);
  std::cout << "[lanczos] n=" << n << " g=" << g_str
            << " decimal_digits=" << decimal_digits << " binary_bits=" << bits
            << std::endl;

  auto t_start = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
  // 获取系统总逻辑核心数
  int total_cores = omp_get_num_procs();
  // 防御性扣除: >= 8 核心截留 2 个，否则截留 1 个核心给操作系统
  int reserve_cores = (total_cores >= 8) ? 2 : 1;
  int target_threads = std::max(1, total_cores - reserve_cores);
  omp_set_num_threads(target_threads);
  std::cout << "[lanczos] CPU Protection: Using " << target_threads << " / "
            << total_cores << " OpenMP threads." << std::endl;
#endif

  // 为多线程任务准备 EMA 滑动窗口状态
  auto omp_last_time = t_start;
  int omp_last_count = 0;
  double omp_ema_ms = 0.0;

  // 进度条辅助函数: 引入多线程安全的计数器和进度控制
  // 增加 complexity 参数 (p): 代表单步耗时 t(i) ~ i^p 增长
  auto print_prog_omp = [&](const char *name, int &completed_count,
                            int complexity = 1) {
    int current;
#pragma omp atomic capture
    current = ++completed_count;

    int progress = (current * 100) / n;
    int bar_pos = (50 * current) / n;

// 控制台输出上锁，防止多线程重影乱码
#pragma omp critical
    {
      auto now = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> diff_from_last =
          now - omp_last_time;

      if (diff_from_last.count() >= 500.0 || current == n) {
        std::chrono::duration<double, std::milli> elapsed = now - t_start;
        if (elapsed.count() >= 100.0 || current == n) {
          // Mathematical extrapolation relying on step complexity integral
          double work_ratio =
              std::pow(static_cast<double>(current) / n, complexity + 1.0);
          double projected_total_ms = elapsed.count() / work_ratio;
          omp_ema_ms = projected_total_ms -
                       elapsed.count(); // store permanently to prevent jitter
          if (omp_ema_ms < 0)
            omp_ema_ms = 0.0;
        }
        omp_last_time = now;
        omp_last_count = current;
      }

      std::cout << "\r[lanczos] " << name << " [";
      for (int p = 0; p < 50; ++p) {
        if (p < bar_pos)
          std::cout << "=";
        else if (p == bar_pos)
          std::cout << ">";
        else
          std::cout << " ";
      }
      std::cout << "] " << progress << "% (" << current << "/" << n << ") ";

      auto total_elapsed =
          std::chrono::duration_cast<std::chrono::milliseconds>(now - t_start)
              .count();
      if (current < n && total_elapsed >= 500) {
        std::cout << "ETA: " << std::fixed << std::setprecision(1)
                  << (omp_ema_ms / 1000.0) << "s   ";
      } else if (current == n) {
        std::cout << "ETA: 0.0s   ";
      }
      std::cout << std::flush;

      if (current == n) {
        omp_ema_ms = 0.0;
        omp_last_count = 0;
        t_start = std::chrono::high_resolution_clock::
            now(); // 重置绝对起点给下一个矩阵任务
        omp_last_time = t_start;
        std::cout << std::endl;
      }
    }
  };

  // ========== 步骤 0: 预计算组合数缓存以免多线程锁竞争 ==========
  precompute_comb_cache(4 * n); // 组合数最大会用到 C(2i, 2k)，即最大 4n

  // ========== 步骤 1: 构造矩阵 B (n × n) ==========
  std::vector<std::vector<SignedBigInt>> B(n, std::vector<SignedBigInt>(n));
  int b_count = 0;
#pragma omp parallel for schedule(dynamic)
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
    print_prog_omp("Matrix B", b_count,
                   0); // 单步为 N 次 O(1) 查表，故步耗时 O(1), p = 0
  }

  // ========== 步骤 2: 构造 Chebyshev 矩阵 C (n × n) ==========
  std::vector<std::vector<BigFloat>> Cmat(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  int c_count = 0;
#pragma omp parallel for schedule(dynamic)
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
    print_prog_omp("Matrix C", c_count,
                   1); // K循环步数随 i 增加，属于 O(i) 步耗时，p = 1
  }

  // ========== 步骤 3: 构造对角矩阵 D (n × n) ==========
  std::vector<BigFloat> D(n, BigFloat(0, work_bits));
  int d_count = 0;
  D[0] = BigFloat(1, work_bits);
  if (n > 1) {
    D[1] = BigFloat(-1, work_bits);
  }
  print_prog_omp("Matrix D", d_count, 0);
  if (n > 1)
    print_prog_omp("Matrix D", d_count, 0);
  for (int i = 2; i < n; i++) {
    // D[i] = D[i−1] × 2(2i−1) / (i−1)
    D[i] = D[i - 1] * BigFloat(2 * (2 * i - 1), work_bits);
    D[i] = D[i] / BigFloat(i - 1, work_bits);
    print_prog_omp("Matrix D", d_count);
  }

  // ========== 步骤 4: 计算 M = D × B × C ==========
  // 由于 D 是对角矩阵: M[i][j] = D[i] × (B × C)[i][j]
  // 先计算 BC = B × C（标准矩阵乘法）
  // BC[i][j] = Σ_k B[i][k] × C[k][j]
  std::vector<std::vector<BigFloat>> BC(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  int bc_count = 0;
#pragma omp parallel for schedule(dynamic)
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
    print_prog_omp("Matrix BC", bc_count,
                   0); // 矩阵乘法单行计算为独立 O(N^2)，不随 i 增大，p = 0
  }

  // M[i][j] = D[i] × BC[i][j]
  std::vector<std::vector<BigFloat>> M(
      n, std::vector<BigFloat>(n, BigFloat(0, work_bits)));
  int m_count = 0;
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[i][j] = D[i] * BC[i][j];
    }
    print_prog_omp("Matrix M", m_count, 0); // O(N) 单行固定耗时，p = 0
  }
  print_prog_omp("Matrix M", m_count, 0); // 手动强制闭环

  // ========== 步骤 5: 构造向量 F (高精度浮点) ==========
  BigFloat bf_half = BigFloat(1, work_bits) / BigFloat(2, work_bits);

  // 在公式 F[i] = (2i)!/(i!) × exp(i+0.5) / 2^(2i-1) / (g+i+0.5)^i /
  // sqrt(g+i+0.5) 中 我们可以通过代数运算彻底避开 exp(i+0.5)
  // 在每一级的真实计算！ 已知 exp(i+0.5) = e^(i+0.5) 且分母有 (g+i+0.5)^i *
  // sqrt(g+i+0.5) = (g+i+0.5)^(i+0.5) 我们如果把它收纳进来： [ e / (g+i+0.5)
  // ]^(i+0.5) 因为这是实数运算且有 e，为了稳定性我们不再使用原版
  // exp，而是让各线程仅算底数与常数。
  // 但为了保留原定精度并利用现成的底数分离法：我们只计算
  // exp(i+0.5)，由于不希望线程串行污染：
  BigFloat exp_one = BigFloat::exp(BigFloat(1, work_bits));
  BigFloat exp_half = BigFloat::exp(bf_half);

  std::vector<BigFloat> F(n);
  int f_count = 0;
#pragma omp parallel for schedule(dynamic)
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

    // --- exp(i + 0.5) 多线程安全生成 ---
    // e^(i+0.5) = (e^1)^i * e^0.5
    BigFloat exp_term = exp_half;
    if (i > 0) {
      exp_term = exp_term * BigFloat::pow(exp_one, bi);
    }

    // --- 1 / 2^(2i−1): 使用 mul_pow2 高效实现（O(1) 指数调整） ---
    // i=0: 2i-1 = -1, 即乘以 2 (mul_pow2(1))
    // i>0: 2i-1 > 0, 即除以 2^(2i-1) (mul_pow2(-(2i-1)))
    int pow2_shift = -(2 * i - 1);

    // --- (g+i+0.5)^i ---
    BigFloat base_pow_i(1, work_bits);
    if (i > 0) {
      base_pow_i = BigFloat::pow(base, bi);
    }

    // --- √(g+i+0.5) ---
    BigFloat base_sqrt = BigFloat::sqrt(base);

    // --- 组合: F[i] = fact_ratio × exp × 2^pow2_shift / base^i / √base ---
    BigFloat Fi = fact_ratio * exp_term;
    Fi = Fi.mul_pow2(pow2_shift); // 高效乘除 2 的幂
    Fi = Fi / base_pow_i;
    Fi = Fi / base_sqrt;
    F[i] = Fi;
    print_prog_omp("Vector F", f_count,
                   1); // f的超越函数极慢，相当于 p=1 级增长
  }

  // ========== 步骤 6: 计算 P = M × F (矩阵-向量乘法) ==========
  std::vector<BigFloat> P(n);
  int p_count = 0;
#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; i++) {
    BigFloat sum(0, work_bits);
    for (int j = 0; j < n; j++) {
      sum += M[i][j] * F[j];
    }
    P[i] = sum;
    P[i].set_precision(bits);               // 截断到目标精度
    print_prog_omp("Vector P", p_count, 0); // p=0 固定耗时
  }

  auto t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms = t_end - t_start;
  std::cout << "[lanczos] Coefficients computed successfully in " << ms.count()
            << " ms." << std::endl;
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
  int bits = static_cast<int>(std::ceil(decimal_digits * 3.3219281)) +
             std::max(64, decimal_digits / 10);
  int n = static_cast<int>(coeffs.size());
  // 求值时也会根据系数大小有精度衰减，补足对应的保护位
  int guard_bits = 256 + static_cast<int>(2.0 * n * std::log2(n + 1));
  int work_bits = bits + guard_bits;

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

  // 计算 S(z) = p[0] + Σ_{k=1}^{n-1} p[k] / (z + k - 1)
  // 使用 Kahan Summation 阻断灾难性大数相消导致的低位精度丢失
  BigFloat S = coeffs[0];
  S.set_precision(work_bits);
  BigFloat c(0, work_bits); // 误差补偿变量

  for (int k = 1; k < n; k++) {
    BigFloat ck = coeffs[k];
    ck.set_precision(work_bits);
    BigFloat denom = zz + BigFloat(k - 1, work_bits); // z + k - 1
    BigFloat term = ck / denom;

    // Kahan 求和逻辑
    BigFloat y = term - c; // 减去上次累加丢失的低位
    BigFloat t = S + y;    // 尝试累加
    c = (t - S) - y;       // 捕获本次累加中被截出的有效低位
    S = t;                 // 更新真实总和
  }

  // 计算 Γ(z) = (base/e)^(z-0.5) × S(z)
  BigFloat exp_one = BigFloat::exp(BigFloat(1, work_bits));
  BigFloat base_over_e = base / exp_one;
  BigFloat prefix = BigFloat::pow(base_over_e, zz - half); // (base/e)^(z-0.5)

  BigFloat result = prefix * S;
  result.set_precision(bits);
  return result;
}
