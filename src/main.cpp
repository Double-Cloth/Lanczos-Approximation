/**
 * @file main.cpp
 * @brief Lanczos 近似系数计算器的命令行入口程序
 *
 * @details
 * 支持两种运行模式:
 *
 *   模式 1 (Generate): 计算 Lanczos 系数并用 CSV 数据验证精度
 *     用法: lanczos_app <n> <g> <digits> [output_dir] [csv_path]
 *     示例: lanczos_app 7 5.0 16
 *
 *   模式 2 (Evaluate): 从已保存的系数文件中读取系数，计算指定点的 Γ(z)
 *     用法: lanczos_app eval <output_dir> <z_value> [display_digits]
 *     示例: lanczos_app eval output_n7_g5.0_d16 50.5 50
 *
 * 输出文件:
 *   - coefficients.txt : 系数列表及公式说明
 *   - parameters.txt   : 计算参数
 *   - verification.txt : 验证结果
 */

#include "lanczos.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <string>
#include <vector>

namespace fs = std::filesystem;

// ============================================
// CSV 解析辅助函数
// ============================================

/**
 * @brief 解析 CSV 行中的引号字段
 * @param line       CSV 行内容
 * @param z_str      输出: 第一个字段（z 值）
 * @param gamma_str  输出: 第二个字段（Γ(z) 值）
 * @return true 若成功解析两个非空字段
 *
 * @details 支持格式: "value1","value2"
 *          引号内的内容被提取，支持嵌套引号（""）
 */
static bool parse_csv_line(const std::string &line, std::string &z_str,
                           std::string &gamma_str) {
  if (line.empty() || line[0] != '"')
    return false;

  // 提取第一个引号对内的内容
  size_t first_end = line.find('"', 1);
  if (first_end == std::string::npos)
    return false;
  z_str = line.substr(1, first_end - 1);

  // 提取第二个引号对内的内容
  size_t second_start = line.find('"', first_end + 1);
  if (second_start == std::string::npos)
    return false;
  size_t second_end = line.find('"', second_start + 1);
  if (second_end == std::string::npos) {
    // 值可能很长, 取到行尾
    gamma_str = line.substr(second_start + 1);
  } else {
    gamma_str = line.substr(second_start + 1, second_end - second_start - 1);
  }
  return !z_str.empty() && !gamma_str.empty();
}

/**
 * @brief 修正 CSV 中缺少前导 "0." 的期望值
 * @param computed  计算得到的值（用于判断量级）
 * @param expected  CSV 中的原始期望值
 * @return 修正后的期望值字符串
 *
 * @details 某些 CSV 文件中，小于 1 的值可能省略了 "0." 前缀，
 *          例如 "886..." 实际上应该是 "0.886..."
 *          通过检查计算值是否 < 1 来决定是否补上前缀
 */
static std::string fix_expected_value(const std::string &computed,
                                      const std::string &expected) {
  // 若 expected 已包含小数点, 不需要修正
  if (expected.find('.') != std::string::npos)
    return expected;

  // 科学计数法不需要补 "0."
  if (expected.find('E') != std::string::npos ||
      expected.find('e') != std::string::npos)
    return expected;

  // 若 expected 很短 (如 "1", "2", "24"), 它是精确整数, 不修正
  if (expected.size() <= 10)
    return expected;

  // 检查 computed 值是否 < 1
  bool computed_lt_one = false;
  if (computed.size() >= 2) {
    if (computed[0] == '0' && computed[1] == '.')
      computed_lt_one = true;
    if (computed[0] == '-' && computed.size() >= 3 && computed[1] == '0' &&
        computed[2] == '.')
      computed_lt_one = true;
  }

  if (computed_lt_one) {
    // CSV 中的值省略了 "0.", 补上
    return "0." + expected;
  }
  return expected;
}

// ============================================
// 主程序
// ============================================

int main(int argc, char *argv[]) {
  // 参数不足时显示帮助
  if (argc < 3) {
    std::cerr << "Usage:" << std::endl;
    std::cerr << "  Mode 1 (Generate): " << argv[0]
              << " <n> <g> <digits> [output_dir] [csv_path]" << std::endl;
    std::cerr << "  Mode 2 (Evaluate): " << argv[0]
              << " eval <output_dir> <z_value> [display_digits]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Example:" << std::endl;
    std::cerr << "  " << argv[0] << " 7 5.0 16" << std::endl;
    std::cerr << "  " << argv[0] << " eval output_n7_g5.0_d16 50.5"
              << std::endl;
    std::cerr << "  " << argv[0] << " eval output_n7_g5.0_d16 50.5 50"
              << std::endl;
    return 1;
  }

  // ========== 模式 2: 从已有目录读取参数和系数，计算 Γ(z) ==========
  std::string mode_flag = argv[1];
  if (mode_flag == "eval") {
    if (argc < 4) {
      std::cerr << "Error: 'eval' mode requires <output_dir> and <z_value>."
                << std::endl;
      return 1;
    }
    std::string output_dir = argv[2];
    std::string z_str = argv[3];

    // 可选: 自定义显示精度
    int custom_display_digits = -1;
    if (argc >= 5) {
      custom_display_digits = std::atoi(argv[4]);
    }

    // --- 读取参数文件 ---
    std::string params_path = output_dir + "/parameters.txt";
    std::string coeff_path = output_dir + "/coefficients.txt";

    std::ifstream fparams(params_path);
    if (!fparams.is_open()) {
      std::cerr << "Error: cannot open " << params_path << std::endl;
      return 1;
    }

    int n = 0, digits = 0;
    std::string g_str;
    std::string line;
    while (std::getline(fparams, line)) {
      // 解析 "key = value" 格式
      if (line.find("n = ") == 0)
        n = std::stoi(line.substr(4));
      else if (line.find("g = ") == 0)
        g_str = line.substr(4);
      else if (line.find("precision_decimal_digits = ") == 0)
        digits = std::stoi(line.substr(27));
    }
    fparams.close();

    if (n <= 0 || digits <= 0) {
      std::cerr << "Error: invalid parameters found in " << params_path
                << std::endl;
      return 1;
    }

    // 确定显示精度
    int display_digits = (custom_display_digits > 0) ? custom_display_digits
                                                     : std::min(digits + 5, 40);

    // 确定工作精度（取 digits 和 display_digits 的较大值）
    int work_digits = std::max(digits, display_digits);
    int prec_bits = static_cast<int>(std::ceil(work_digits * 3.3219281)) + 64;

    // --- 读取系数文件 ---
    std::vector<BigFloat> coeffs;
    std::ifstream fcoeffs(coeff_path);
    if (!fcoeffs.is_open()) {
      std::cerr << "Error: cannot open " << coeff_path << std::endl;
      return 1;
    }

    while (std::getline(fcoeffs, line)) {
      if (line.empty() || line[0] == '#')
        continue; // 跳过注释行
      size_t comma_pos = line.find(',');
      if (comma_pos != std::string::npos) {
        // 提取逗号后的系数值
        std::string val_str = line.substr(comma_pos + 1);
        size_t first_non_space = val_str.find_first_not_of(" \t");
        if (first_non_space != std::string::npos) {
          val_str = val_str.substr(first_non_space);
        }
        coeffs.push_back(BigFloat::from_string(val_str, prec_bits));
      }
    }
    fcoeffs.close();

    // 验证系数数量
    if (coeffs.size() != static_cast<size_t>(n)) {
      std::cerr << "Error: expected " << n << " coefficients, but found "
                << coeffs.size() << std::endl;
      return 1;
    }

    // --- 计算并输出 Γ(z) ---
    BigFloat z = BigFloat::from_string(z_str, prec_bits);
    BigFloat gamma_val = lanczos_gamma(z, coeffs, g_str, work_digits);

    std::cout << "Gamma(" << z_str
              << ") = " << gamma_val.to_decimal_string(display_digits)
              << std::endl;

    return 0;
  }

  // ========== 模式 1: 计算系数并验证 ==========
  if (argc < 4) {
    std::cerr << "Error: Generate mode requires <n> <g> <digits>." << std::endl;
    return 1;
  }

  int n = std::atoi(argv[1]);
  std::string g_str = argv[2]; // g 始终保持字符串形式，避免 double 精度丢失
  int digits = std::atoi(argv[3]);

  if (n <= 0 || digits <= 0) {
    std::cerr << "Error: n and digits must be positive integers." << std::endl;
    return 1;
  }

  // 确定输出目录（默认: output_n{n}_g{g}_d{digits}）
  std::string output_dir;
  if (argc >= 5) {
    output_dir = argv[4];
  } else {
    // 使用 g_str 而非 double 避免精度丢失
    output_dir = "output_n" + std::to_string(n) + "_g" + g_str + "_d" +
                 std::to_string(digits);
  }

  // CSV 验证文件路径
  std::string csv_path = "../assets/real_gamma.csv";
  if (argc >= 6) {
    csv_path = argv[5];
  }

  // 创建输出目录
  fs::create_directories(output_dir);

  std::cout << "=== Lanczos Coefficient Calculator ===" << std::endl;
  std::cout << "Parameters: n=" << n << ", g=" << g_str
            << ", precision=" << digits << " decimal digits" << std::endl;
  std::cout << "Output directory: " << output_dir << std::endl;
  std::cout << std::endl;

  // --- 计算系数 ---
  auto coeffs = compute_lanczos_coefficients(n, g_str, digits);

  // 输出系数到控制台
  std::cout << std::endl;
  std::cout << "=== Computed Coefficients ===" << std::endl;
  for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
    std::cout << "  p[" << i << "] = " << coeffs[i].to_decimal_string(digits)
              << std::endl;
  }

  // --- 保存系数文件 ---
  {
    std::string coeff_path = output_dir + "/coefficients.txt";
    std::ofstream fout(coeff_path);
    if (!fout.is_open()) {
      std::cerr << "Error: cannot open " << coeff_path << " for writing."
                << std::endl;
      return 1;
    }
    fout << "# Lanczos Approximation Coefficients" << std::endl;
    fout << "# n = " << n << std::endl;
    fout << "# g = " << g_str << std::endl;
    fout << "# precision = " << digits << " decimal digits" << std::endl;
    fout << "#" << std::endl;
    fout << "# Formula:" << std::endl;
    fout << "#   Gamma(z) = (base/e)^(z-0.5) * S(z)" << std::endl;
    fout << "#   where base = z + g - 0.5" << std::endl;
    fout << "#         S(z) = p[0] + sum_{k=1}^{n-1} p[k] / (z + k - 1)"
         << std::endl;
    fout << "#   Note: sqrt(2*pi) factor is absorbed into p[0]" << std::endl;
    fout << "#" << std::endl;
    fout << "# index, coefficient" << std::endl;
    for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
      fout << i << ", " << coeffs[i].to_decimal_string(digits) << std::endl;
    }
    fout.close();
    std::cout << "\nCoefficients written to: " << coeff_path << std::endl;
  }

  // --- 保存参数文件 ---
  {
    std::string params_path = output_dir + "/parameters.txt";
    std::ofstream fout(params_path);
    fout << "n = " << n << std::endl;
    fout << "g = " << g_str << std::endl;
    fout << "precision_decimal_digits = " << digits << std::endl;
    int bits = static_cast<int>(std::ceil(digits * 3.3219281)) + 64;
    fout << "precision_binary_bits = " << bits << std::endl;
    fout.close();
  }

  // ========== CSV 验证 ==========
  std::cout << "\n=== Loading test data from " << csv_path
            << " ===" << std::endl;

  std::ifstream csv_file(csv_path);
  if (!csv_file.is_open()) {
    std::cerr << "Warning: cannot open " << csv_path
              << ", skipping CSV verification." << std::endl;
    std::cerr << "You can specify the CSV path as the 5th argument."
              << std::endl;
    return 0;
  }

  // 加载所有测试数据
  struct TestEntry {
    std::string z_str;
    std::string expected_str;
  };
  std::vector<TestEntry> test_data;

  std::string line;
  while (std::getline(csv_file, line)) {
    // 去除行尾的 \r（处理 Windows 换行符）
    if (!line.empty() && line.back() == '\r')
      line.pop_back();

    std::string z_str, gamma_str;
    if (parse_csv_line(line, z_str, gamma_str)) {
      test_data.push_back({z_str, gamma_str});
    }
  }
  csv_file.close();

  std::cout << "Loaded " << test_data.size() << " test entries." << std::endl;

  // 筛选测试点: z ∈ (0, 50]
  std::vector<TestEntry> selected;
  for (auto &entry : test_data) {
    double z_val = std::atof(entry.z_str.c_str());
    if (z_val > 0 && z_val <= 50) {
      selected.push_back(entry);
    }
  }
  std::cout << "Selected " << selected.size() << " test points (z <= 50)."
            << std::endl;

  // --- 逐点验证 ---
  // 动态阈值: 基于 Lanczos 近似精度和工作精度自动计算
  // Lanczos 近似的有效数字数约为 n+3（取决于 n 和 g 的选择）
  // 实际精度上限为 min(digits, n+3)，阈值为此精度允许末 2 位偏差
  // 例: n=7,digits=16 → effective=10 → 1e-8%
  //     n=13,digits=25 → effective=16 → 1e-14%
  //     n=24,digits=40 → effective=27 → 1e-25%
  int effective_digits = std::min(digits, n + 3);
  double threshold_percent = std::pow(10.0, -(effective_digits - 2));

  std::cout << "\n=== Verification ===" << std::endl;
  std::cout << "  Pass threshold: relative error <= " << std::scientific
            << threshold_percent << "%" << std::defaultfloat << " (based on "
            << digits << "-digit precision)" << std::endl;
  std::cout << std::endl;

  std::string verify_path = output_dir + "/verification.txt";
  std::ofstream vout(verify_path);
  vout << "# Lanczos Gamma Function Verification Results" << std::endl;
  vout << "# Parameters: n=" << n << ", g=" << g_str << ", digits=" << digits
       << std::endl;
  vout << "# Test data source: " << csv_path << std::endl;
  vout << "# Pass threshold: relative error <= " << std::scientific
       << threshold_percent << "%" << std::defaultfloat << std::endl;
  vout << "#" << std::endl;
  vout << "# z, computed, expected (first 30 digits), relative_error%, status"
       << std::endl;

  int total_tests = static_cast<int>(selected.size());
  int pass_count = 0;
  double max_error = 0.0;
  double min_error = 1e30;
  double sum_error = 0.0;

  int prec_bits = static_cast<int>(std::ceil(digits * 3.3219281)) + 64;
  int display_digits = digits;

  for (int idx = 0; idx < total_tests; idx++) {
    auto &entry = selected[idx];

    // 计算 Γ(z)
    BigFloat z = BigFloat::from_string(entry.z_str, prec_bits);
    BigFloat gamma_val = lanczos_gamma(z, coeffs, g_str, digits);
    std::string computed_str = gamma_val.to_decimal_string(display_digits);

    // --- 准备比较: 截断期望值避免解析几万位造成的无谓性能开销 ---
    int max_needed_len = digits + 30;
    std::string truncated_expected = entry.expected_str;
    size_t e_pos = entry.expected_str.find_first_of("eE");
    if (e_pos != std::string::npos) {
      std::string exp_part = entry.expected_str.substr(e_pos);
      if (e_pos > (size_t)max_needed_len) {
        truncated_expected =
            entry.expected_str.substr(0, max_needed_len) + exp_part;
      }
    } else {
      if (entry.expected_str.length() > (size_t)max_needed_len) {
        truncated_expected = entry.expected_str.substr(0, max_needed_len);
      }
    }

    int expected_prec_bits =
        static_cast<int>(truncated_expected.length() * 3.3219281) + 64;
    int compare_prec =
        std::min(prec_bits + 64, std::max(prec_bits, expected_prec_bits));

    BigFloat gamma_val_ext = gamma_val;
    gamma_val_ext.set_precision(compare_prec);

    // 修正期望值（处理缺失的 "0." 前缀）
    std::string fixed_expected = fix_expected_value(
        gamma_val.to_decimal_string(10, false), truncated_expected);
    BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);

    // 计算绝对误差
    BigFloat abs_err = (gamma_val_ext - expected_val).abs();

    // 取期望值的前 30 位用于显示
    std::string expected_short = entry.expected_str.substr(
        0, std::min(static_cast<size_t>(30), entry.expected_str.size()));

    // --- 计算相对误差百分比 ---
    // relative_error = |computed - expected| / |expected| × 100%
    double rel_err_percent = 0.0;
    if (expected_val.is_zero()) {
      rel_err_percent = abs_err.is_zero() ? 0.0 : 100.0;
    } else if (abs_err.is_zero()) {
      rel_err_percent = 0.0;
    } else {
      BigFloat abs_expected = expected_val.abs();
      BigFloat ratio = abs_err / abs_expected;
      std::string ratio_str = ratio.to_decimal_string(20, true);
      double ratio_val = std::stod(ratio_str);
      rel_err_percent = ratio_val * 100.0;
    }

    // 判定通过/失败
    std::string status =
        (rel_err_percent <= threshold_percent) ? "PASS" : "FAIL";
    if (rel_err_percent <= threshold_percent)
      pass_count++;

    // 更新统计量
    sum_error += rel_err_percent;
    if (rel_err_percent > max_error)
      max_error = rel_err_percent;
    if (rel_err_percent < min_error)
      min_error = rel_err_percent;

    // 输出到控制台和文件
    std::cout << "  Gamma(" << entry.z_str << ") = " << computed_str
              << "  [relative error: " << std::scientific << rel_err_percent
              << "%, " << std::defaultfloat << status << "]" << std::endl;

    vout << entry.z_str << ", " << computed_str << ", " << expected_short
         << "..."
         << ", " << std::scientific << rel_err_percent << "%"
         << std::defaultfloat << ", " << status << std::endl;
  }

  // --- 输出统计摘要 ---
  vout << "#" << std::endl;
  vout << "# ==============================" << std::endl;
  vout << "# Summary" << std::endl;
  vout << "# ==============================" << std::endl;
  vout << "# Total tests:        " << total_tests << std::endl;
  vout << "# Passed (relative error <= " << std::scientific << threshold_percent
       << "%): " << std::defaultfloat << pass_count << "/" << total_tests
       << std::endl;
  if (total_tests > 0) {
    vout << "# Max relative error: " << std::scientific << max_error << "%"
         << std::defaultfloat << std::endl;
    vout << "# Min relative error: " << std::scientific << min_error << "%"
         << std::defaultfloat << std::endl;
    vout << "# Avg relative error: " << std::scientific
         << sum_error / total_tests << "%" << std::defaultfloat << std::endl;
  }
  vout.close();

  std::cout << "\n=== Summary ===" << std::endl;
  std::cout << "  Total tests:         " << total_tests << std::endl;
  std::cout << "  Passed (relative error <= " << std::scientific
            << threshold_percent << "%): " << std::defaultfloat << pass_count
            << "/" << total_tests << std::endl;
  if (total_tests > 0) {
    std::cout << "  Max relative error:  " << std::scientific << max_error
              << "%" << std::defaultfloat << std::endl;
    std::cout << "  Min relative error:  " << std::scientific << min_error
              << "%" << std::defaultfloat << std::endl;
    std::cout << "  Avg relative error:  " << std::scientific
              << sum_error / total_tests << "%" << std::defaultfloat
              << std::endl;
  }

  std::cout << "\n=== Output Files ===" << std::endl;
  std::cout << "  " << output_dir
            << "/coefficients.txt   - Lanczos coefficients" << std::endl;
  std::cout << "  " << output_dir
            << "/parameters.txt     - Computation parameters" << std::endl;
  std::cout << "  " << output_dir
            << "/verification.txt   - Gamma function test results" << std::endl;

  return 0;
}
