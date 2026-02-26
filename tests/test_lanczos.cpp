/**
 * @file test_lanczos.cpp
 * @brief Lanczos 近似精度验证测试
 *
 * @details
 * 本测试程序从 CSV 文件（real_gamma.csv）读取已知的 Γ(z) 高精度值，
 * 与 Lanczos 近似计算的结果进行比较，使用 **相对误差百分比** 作为判定标准。
 *
 * CSV 格式: "z","Gamma(z)"
 *   其中 Gamma(z) 是高精度十进制字符串
 *
 * 验证方法: 相对误差百分比
 *   relative_error = |computed - expected| / |expected| × 100%
 *   当相对误差 ≤ 阈值（默认 1e-6%，即约 8 位有效数字）时判定为 PASS。
 *
 * 命令行参数:
 *   test_lanczos [n] [g] [precision] [csv_path] [max_tests] [threshold%]
 *   默认: n=7, g=5.0, precision=16, csv_path=real_gamma.csv, max_tests=50,
 *         threshold=1e-6%
 *
 * 通过标准: 相对误差 ≤ threshold% 为 PASS
 */

#include "lanczos.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <string>
#include <vector>

/**
 * @brief 解析 CSV 行中的一个字段（支持引号包裹）
 * @param line CSV 行内容
 * @param pos  当前解析位置（输入/输出参数）
 * @return 解析出的字段内容
 *
 * @details 支持两种格式:
 *   1. 引号包裹: "value" — 内部的 "" 被转义为单个 "
 *   2. 无引号: 读取直到逗号或行尾
 *   解析完成后 pos 指向下一个字段的起始位置
 */
std::string parse_csv_field(const std::string &line, size_t &pos) {
  std::string result;
  if (pos >= line.size())
    return result;

  if (line[pos] == '"') {
    // 引号包裹的字段
    pos++; // 跳过开头引号
    while (pos < line.size()) {
      if (line[pos] == '"') {
        if (pos + 1 < line.size() && line[pos + 1] == '"') {
          // 转义的引号 "" → "
          result += '"';
          pos += 2;
        } else {
          // 结尾引号
          pos++;
          break;
        }
      } else {
        result += line[pos++];
      }
    }
  } else {
    // 无引号字段: 读取直到逗号
    while (pos < line.size() && line[pos] != ',') {
      result += line[pos++];
    }
  }

  // 跳过字段分隔逗号
  if (pos < line.size() && line[pos] == ',')
    pos++;
  return result;
}

/**
 * @brief 测试主函数
 *
 * @details 流程:
 *   1. 解析命令行参数
 *   2. 计算 Lanczos 系数
 *   3. 逐行读取 CSV 文件
 *   4. 对每个测试点:
 *      a. 跳过 z ≤ 0.5（反射公式未实现）
 *      b. 计算 Γ(z)
 *      c. 使用纯字符串比对法统计匹配位数
 *   5. 输出统计摘要
 */
int main(int argc, char *argv[]) {
  // === 参数解析 ===
  int n = 7;                                         // Lanczos 级数项数
  std::string g_str = "5.0";                         // Lanczos 参数 g
  int precision = 16;                                // 目标精度（十进制位数）
  std::string csv_path = "../assets/real_gamma.csv"; // CSV 验证文件路径
  int max_tests = 50;                                // 最多测试点数
  double threshold_percent =
      1e-6; // 相对误差阈值百分比（默认 1e-6% ≈ 8位有效数字）

  if (argc >= 4) {
    n = std::atoi(argv[1]);
    g_str = argv[2];
    precision = std::atoi(argv[3]);
  }
  if (argc >= 5) {
    csv_path = argv[4];
  }
  if (argc >= 6) {
    max_tests = std::atoi(argv[5]);
  }
  if (argc >= 7) {
    threshold_percent = std::atof(argv[6]);
  }

  // === 显示参数信息 ===
  std::cout << "=== Lanczos Approximation Test ===" << std::endl;
  std::cout << "Parameters: n=" << n << ", g=" << g_str
            << ", precision=" << precision << " decimal digits" << std::endl;
  std::cout << "CSV file: " << csv_path << std::endl;
  std::cout << "Pass threshold: relative error <= " << std::scientific
            << threshold_percent << "%" << std::defaultfloat << std::endl;
  std::cout << std::endl;

  // === 计算 Lanczos 系数 ===
  std::cout << "--- Computing Lanczos Coefficients ---" << std::endl;
  auto coeffs = compute_lanczos_coefficients(n, g_str, precision);

  // 输出系数（多给 5 位以观察截断误差）
  std::cout << std::endl << "Coefficients:" << std::endl;
  for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
    std::cout << "  p[" << i
              << "] = " << coeffs[i].to_decimal_string(precision + 5)
              << std::endl;
  }

  // === 读取 CSV 验证文件 ===
  std::ifstream file(csv_path);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open " << csv_path << std::endl;
    return 1;
  }

  // === 输出表头 ===
  std::cout << std::endl << "--- Gamma Function Verification ---" << std::endl;
  std::cout << std::setw(12) << "z" << std::setw(22) << "relative_error%"
            << std::setw(10) << "status" << std::endl;
  std::cout << std::string(44, '-') << std::endl;

  // === 逐点验证 ===
  int total_tests = 0;     // 总测试点数
  int passed = 0;          // 通过数
  int failed = 0;          // 失败数
  double max_error = 0.0;  // 最大相对误差
  double min_error = 1e30; // 最小相对误差
  double sum_error = 0.0;  // 相对误差之和（用于计算平均值）

  std::string line;
  while (std::getline(file, line) && total_tests < max_tests) {
    if (line.empty())
      continue;

    // --- 解析 CSV 行: "z","gamma_value" ---
    size_t pos = 0;
    std::string z_str = parse_csv_field(line, pos);
    std::string expected_str = parse_csv_field(line, pos);

    if (z_str.empty() || expected_str.empty())
      continue;

    // --- 解析 z 值为 BigFloat ---
    int prec_bits = static_cast<int>(std::ceil(precision * 3.3219281)) + 64;
    BigFloat z;
    try {
      z = BigFloat::from_string(z_str, prec_bits);
    } catch (...) {
      continue; // 无法解析的 z 值直接跳过
    }

    // 已支持反射公式，不再跳过 z <= 0.5 的验证数据

    // --- 计算 Γ(z) ---
    BigFloat computed = lanczos_gamma(z, coeffs, g_str, precision);

    // --- 准备比较: 限制 expected_str 长度以避免解析 3
    // 万位超长数字造成的极大开销 ---
    int max_needed_len = precision + 30; // 包含符号和科学记数法的充足余量
    std::string truncated_expected = expected_str;
    size_t e_pos = expected_str.find_first_of("eE");
    if (e_pos != std::string::npos) {
      // 科学计数法：保留有效数字部分 + 指数部分
      std::string exp_part = expected_str.substr(e_pos);
      if (e_pos > (size_t)max_needed_len) {
        truncated_expected = expected_str.substr(0, max_needed_len) + exp_part;
      }
    } else {
      if (expected_str.length() > (size_t)max_needed_len) {
        truncated_expected = expected_str.substr(0, max_needed_len);
      }
    }

    int expected_prec_bits =
        static_cast<int>(truncated_expected.length() * 3.3219281) +
        std::max(64, static_cast<int>(truncated_expected.length()) / 10);
    int compare_prec = std::min(prec_bits + std::max(64, prec_bits / 10),
                                std::max(prec_bits, expected_prec_bits));

    BigFloat computed_ext = computed;
    computed_ext.set_precision(compare_prec);

    // 修正 CSV 中可能缺失的 "0." 前缀
    std::string fixed_expected = truncated_expected;
    bool is_scientific = (truncated_expected.find('E') != std::string::npos ||
                          truncated_expected.find('e') != std::string::npos);
    if (computed < BigFloat(1, prec_bits) && truncated_expected.length() > 5 &&
        !is_scientific) {
      if (truncated_expected.substr(0, 2) != "0.") {
        if (truncated_expected[0] == '.') {
          fixed_expected = "0" + truncated_expected; // ".886" → "0.886"
        } else {
          fixed_expected = "0." + truncated_expected; // "886..." → "0.886..."
        }
      }
    }

    BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);
    BigFloat abs_err = (computed_ext - expected_val).abs();

    // --- 计算相对误差百分比 ---
    // relative_error = |computed - expected| / |expected| × 100%
    double rel_err_percent = 0.0;
    if (expected_val.is_zero()) {
      rel_err_percent = abs_err.is_zero() ? 0.0 : 100.0;
    } else if (abs_err.is_zero()) {
      rel_err_percent = 0.0; // 完全匹配
    } else {
      BigFloat abs_expected = expected_val.abs();
      BigFloat ratio = abs_err / abs_expected;
      std::string ratio_str = ratio.to_decimal_string(20, true);
      try {
        double ratio_val = std::stod(ratio_str);
        rel_err_percent = ratio_val * 100.0;
      } catch (const std::out_of_range &) {
        rel_err_percent = 0.0;
      }
    }

    // 更新统计量
    sum_error += rel_err_percent;
    if (rel_err_percent > max_error)
      max_error = rel_err_percent;
    if (rel_err_percent < min_error)
      min_error = rel_err_percent;

    // --- 判定通过/失败 ---
    std::string status;
    if (rel_err_percent <= threshold_percent) {
      status = "PASS";
      passed++;
    } else {
      status = "FAIL";
      failed++;
    }

    // 输出当前测试点结果
    std::cout << std::setw(12) << z_str << std::setw(22) << std::scientific
              << rel_err_percent << std::defaultfloat << std::setw(10) << status
              << std::endl;

    total_tests++;
  }

  // === 输出统计摘要 ===
  std::cout << std::endl << "=== Summary ===" << std::endl;
  std::cout << "Total tests: " << total_tests << std::endl;
  std::cout << "Passed: " << passed << " (relative error <= " << std::scientific
            << threshold_percent << "%)" << std::defaultfloat << std::endl;
  std::cout << "Failed: " << failed << std::endl;
  if (total_tests > 0) {
    std::cout << "Max relative error:  " << std::scientific << max_error << "%"
              << std::defaultfloat << std::endl;
    std::cout << "Min relative error:  " << std::scientific << min_error << "%"
              << std::defaultfloat << std::endl;
    std::cout << "Avg relative error:  " << std::scientific
              << sum_error / total_tests << "%" << std::defaultfloat
              << std::endl;
  }

  return failed > 0 ? 1 : 0;
}
