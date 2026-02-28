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
static std::string parse_csv_field(const std::string &line, size_t &pos) {
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
// 共享验证逻辑
// ============================================

/**
 * @brief 执行统一的 CSV 验证逻辑
 * @param csv_path      CSV 文件路径
 * @param coeffs        计算得到的 Lanczos 系数
 * @param n             Lanczos 级数项数
 * @param g_str         Lanczos 参数 g (字符串形式)
 * @param precision     验证精度（十进制位数）
 * @param max_tests     最大测试点数，若为负则测试所有不过滤数量的点
 * @param threshold_str 相对误差阈值百分比 (如 "1e-6")
 * @param output_dir    输出目录，若指定则将结果同时写入 "verification.txt"
 * @param filter_z_le_50 是否只测试 z <= 50 的点
 */
static void run_csv_verification(const std::string &csv_path,
                                 const std::vector<BigFloat> &coeffs, int n,
                                 const std::string &g_str, int precision,
                                 int max_tests,
                                 const std::string &threshold_str,
                                 const std::string &output_dir,
                                 bool filter_z_le_50) {
  std::cout << "\n=== Loading test data from " << csv_path
            << " ===" << std::endl;

  std::ifstream csv_file(csv_path);
  if (!csv_file.is_open()) {
    std::cerr << "Warning: cannot open " << csv_path
              << ", skipping CSV verification." << std::endl;
    return;
  }

  // 加载并可能筛选测试数据
  struct TestEntry {
    std::string z_str;
    std::string expected_str;
  };
  std::vector<TestEntry> test_data;

  std::string line;
  while (std::getline(csv_file, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();

    std::string z_str, gamma_str;
    if (parse_csv_line(line, z_str, gamma_str)) {
      test_data.push_back({z_str, gamma_str});
    }
  }
  csv_file.close();

  std::cout << "Loaded " << test_data.size() << " test entries." << std::endl;

  std::vector<TestEntry> selected;
  if (filter_z_le_50) {
    for (auto &entry : test_data) {
      BigFloat z_val = BigFloat::from_string(entry.z_str, 64);
      if (!z_val.is_negative() && !z_val.is_zero() &&
          z_val <= BigFloat(50, 64)) {
        selected.push_back(entry);
      }
    }
    std::cout << "Selected " << selected.size() << " test points (z <= 50)."
              << std::endl;
  } else {
    for (size_t i = 0;
         i < test_data.size() && (max_tests < 0 || (int)i < max_tests); i++) {
      selected.push_back(test_data[i]);
    }
  }

  int prec_bits = static_cast<int>(std::ceil(precision * 3.3219281)) + 64;
  BigFloat threshold = BigFloat::from_string(threshold_str, prec_bits);

  std::cout << "\n=== Verification ===" << std::endl;
  std::cout << "  Pass threshold: relative error <= " << threshold_str
            << "% (based on " << precision << "-digit precision)" << std::endl;
  std::cout << std::endl;

  std::ofstream vout;
  if (!output_dir.empty()) {
    std::string verify_path = output_dir + "/verification.txt";
    vout.open(verify_path);
    vout << "# Lanczos Gamma Function Verification Results" << std::endl;
    vout << "# Parameters: n=" << n << ", g=" << g_str
         << ", digits=" << precision << std::endl;
    vout << "# Test data source: " << csv_path << std::endl;
    vout << "# Pass threshold: relative error <= " << threshold_str << "%"
         << std::endl;
    vout << "#" << std::endl;
    vout << "# z, computed, expected (first 30 digits), abs_error, "
            "relative_error%"
         << std::endl;
  }

  std::cout << std::setw(12) << "z" << std::setw(30) << "computed"
            << std::setw(30) << "expected" << std::setw(22) << "abs_error"
            << std::setw(22) << "relative_error%" << std::endl;
  std::cout << std::string(116, '-') << std::endl;

  int total_tests = static_cast<int>(selected.size());
  int pass_count = 0;
  int failed_count = 0;
  BigFloat max_error(0, prec_bits);
  BigFloat min_error = BigFloat::from_string("1e30", prec_bits);
  BigFloat sum_error(0, prec_bits);
  BigFloat max_abs_error(0, prec_bits);
  BigFloat min_abs_error = BigFloat::from_string("1e30", prec_bits);
  BigFloat sum_abs_error(0, prec_bits);

  int display_digits = precision;

  for (int idx = 0; idx < total_tests; idx++) {
    auto &entry = selected[idx];

    BigFloat z = BigFloat::from_string(entry.z_str, prec_bits);
    BigFloat gamma_val = lanczos_gamma(z, coeffs, g_str, precision);
    std::string computed_str = gamma_val.to_decimal_string(display_digits);

    int max_needed_len = precision + 30;
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
        static_cast<int>(truncated_expected.length() * 3.3219281) +
        std::max(64, static_cast<int>(truncated_expected.length()) / 10);
    int compare_prec = std::min(prec_bits + std::max(64, prec_bits / 10),
                                std::max(prec_bits, expected_prec_bits));

    BigFloat gamma_val_ext = gamma_val;
    gamma_val_ext.set_precision(compare_prec);

    std::string fixed_expected = fix_expected_value(
        gamma_val.to_decimal_string(10, false), truncated_expected);
    BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);
    BigFloat abs_err = (gamma_val_ext - expected_val).abs();

    std::string expected_short = entry.expected_str.substr(
        0, std::min(static_cast<size_t>(30), entry.expected_str.size()));

    BigFloat rel_err_percent(0, compare_prec);
    if (expected_val.is_zero()) {
      rel_err_percent = abs_err.is_zero() ? BigFloat(0, compare_prec)
                                          : BigFloat(100, compare_prec);
    } else if (abs_err.is_zero()) {
      rel_err_percent = BigFloat(0, compare_prec);
    } else {
      BigFloat abs_expected = expected_val.abs();
      BigFloat ratio = abs_err / abs_expected;
      rel_err_percent = ratio * BigFloat(100, compare_prec);
    }

    std::string rel_err_str = rel_err_percent.to_decimal_string(6, true);

    BigFloat threshold_prec = threshold;
    threshold_prec.set_precision(compare_prec);
    if (rel_err_percent <= threshold_prec)
      pass_count++;
    else
      failed_count++;

    sum_error += rel_err_percent;
    if (rel_err_percent > max_error)
      max_error = rel_err_percent;
    if (rel_err_percent < min_error)
      min_error = rel_err_percent;

    sum_abs_error += abs_err;
    if (abs_err > max_abs_error)
      max_abs_error = abs_err;
    if (abs_err < min_abs_error)
      min_abs_error = abs_err;

    std::string abs_err_str = abs_err.to_decimal_string(6, true);
    std::string computed_out_str = gamma_val.to_decimal_string(15, true);
    std::string expected_out_str = expected_val.to_decimal_string(15, true);

    std::cout << "\r" << std::string(80, ' ') << "\r";

    if (!output_dir.empty()) {
      std::cout << "  Gamma(" << entry.z_str << "):" << std::endl
                << "    computed = " << computed_str << std::endl
                << "    expected = " << expected_short << "..." << std::endl
                << "    abs err  = " << abs_err_str << std::endl
                << "    rel err  = " << rel_err_str << "%" << std::endl;
      vout << entry.z_str << ", " << computed_str << ", " << expected_short
           << "..."
           << ", " << abs_err_str << ", " << rel_err_str << "%" << std::endl;
    } else {
      std::cout << std::setw(12) << entry.z_str << std::setw(30)
                << computed_out_str << std::setw(30) << expected_out_str
                << std::setw(22) << abs_err_str << std::setw(22) << rel_err_str
                << std::endl;
    }

    int current = idx + 1;
    int progress_percent = (current * 100) / total_tests;
    int bar_pos = (50 * current) / total_tests;
    std::cout << "[";
    for (int i = 0; i < 50; ++i) {
      if (i < bar_pos)
        std::cout << "=";
      else if (i == bar_pos)
        std::cout << ">";
      else
        std::cout << " ";
    }
    std::cout << "] " << progress_percent << " % (" << current << "/"
              << total_tests << ")" << std::flush;
  }

  std::cout << "\r" << std::string(80, ' ') << "\r";

  if (!output_dir.empty() && vout.is_open()) {
    vout << "#" << std::endl;
    vout << "# ==============================" << std::endl;
    vout << "# Summary" << std::endl;
    vout << "# ==============================" << std::endl;
    vout << "# Total tests:        " << total_tests << std::endl;
    vout << "# Passed (relative error <= " << threshold_str
         << "%): " << pass_count << "/" << total_tests << std::endl;
    if (total_tests > 0) {
      BigFloat avg_error = sum_error / BigFloat(total_tests, prec_bits);
      BigFloat avg_abs_error = sum_abs_error / BigFloat(total_tests, prec_bits);
      vout << "# Max absolute error: "
           << max_abs_error.to_decimal_string(6, true) << std::endl;
      vout << "# Min absolute error: "
           << min_abs_error.to_decimal_string(6, true) << std::endl;
      vout << "# Avg absolute error: "
           << avg_abs_error.to_decimal_string(6, true) << std::endl;
      vout << "# Max relative error: " << max_error.to_decimal_string(6, true)
           << "%" << std::endl;
      vout << "# Min relative error: " << min_error.to_decimal_string(6, true)
           << "%" << std::endl;
      vout << "# Avg relative error: " << avg_error.to_decimal_string(6, true)
           << "%" << std::endl;
    }
    vout.close();
  }

  std::cout << "\n=== Summary ===" << std::endl;
  std::cout << "  Total tests:         " << total_tests << std::endl;
  std::cout << "  Passed (relative error <= " << threshold_str
            << "%): " << pass_count << "/" << total_tests << std::endl;
  if (total_tests > 0) {
    BigFloat avg_error = sum_error / BigFloat(total_tests, prec_bits);
    BigFloat avg_abs_error = sum_abs_error / BigFloat(total_tests, prec_bits);
    std::cout << "\n=== absolute error ===" << std::endl;
    std::cout << "  Max absolute error:  "
              << max_abs_error.to_decimal_string(6, true) << std::endl;
    std::cout << "  Min absolute error:  "
              << min_abs_error.to_decimal_string(6, true) << std::endl;
    std::cout << "  Avg absolute error:  "
              << avg_abs_error.to_decimal_string(6, true) << std::endl;
    std::cout << "=== relative error ===" << std::endl;
    std::cout << "  Max relative error:  "
              << max_error.to_decimal_string(6, true) << "%" << std::endl;
    std::cout << "  Min relative error:  "
              << min_error.to_decimal_string(6, true) << "%" << std::endl;
    std::cout << "  Avg relative error:  "
              << avg_error.to_decimal_string(6, true) << "%" << std::endl;
  }
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
    std::cerr << "  Mode 3 (Test):     " << argv[0]
              << " test <n> <g> <digits> [csv_path] [max_tests] [threshold%]"
              << std::endl;
    std::cerr << std::endl;
    std::cerr << "Example:" << std::endl;
    std::cerr << "  " << argv[0] << " 7 5.0 16" << std::endl;
    std::cerr << "  " << argv[0] << " eval output_n7_g5.0_d16 50.5"
              << std::endl;
    std::cerr << "  " << argv[0] << " eval output_n7_g5.0_d16 50.5 50"
              << std::endl;
    std::cerr << "  " << argv[0]
              << " test 7 5.0 16 ../assets/real_gamma.csv 50 1e-6" << std::endl;
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

  // ========== 模式 3: 测试并输出精度验证表格 ==========
  if (mode_flag == "test") {
    if (argc < 5) {
      std::cerr << "Error: 'test' mode requires <n> <g> <digits>." << std::endl;
      return 1;
    }

    int n = std::atoi(argv[2]);
    std::string g_str = argv[3];
    int precision = std::atoi(argv[4]);
    std::string csv_path = "../assets/real_gamma.csv";
    int max_tests = 50;
    std::string threshold_str = "1e-6";

    if (argc >= 6)
      csv_path = argv[5];
    if (argc >= 7)
      max_tests = std::atoi(argv[6]);
    if (argc >= 8)
      threshold_str = argv[7];

    if (n <= 0 || precision <= 0) {
      std::cerr << "Error: n and digits must be positive integers."
                << std::endl;
      return 1;
    }

    // === 显示参数信息 ===
    std::cout << "=== Lanczos Approximation Test ===" << std::endl;
    std::cout << "Parameters: n=" << n << ", g=" << g_str
              << ", precision=" << precision << " decimal digits" << std::endl;
    std::cout << "CSV file: " << csv_path << std::endl;
    int prec_bits = static_cast<int>(std::ceil(precision * 3.3219281)) + 64;
    BigFloat threshold = BigFloat::from_string(threshold_str, prec_bits);

    std::cout << "Pass threshold: relative error <= " << threshold_str << "%"
              << std::endl;
    std::cout << std::endl;

    // === 计算 Lanczos 系数 ===
    std::cout << "--- Computing Lanczos Coefficients ---" << std::endl;
    auto coeffs = compute_lanczos_coefficients(n, g_str, precision);

    // 输出系数
    std::cout << std::endl << "Coefficients:" << std::endl;
    for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
      std::cout << "  p[" << i
                << "] = " << coeffs[i].to_decimal_string(precision + 5)
                << std::endl;
    }

    // 调用通用验证函数
    run_csv_verification(csv_path, coeffs, n, g_str, precision, max_tests,
                         threshold_str, "", false);

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

  // 动态阈值: 基于 Lanczos 近似精度和工作精度自动计算
  int prec_bits = static_cast<int>(std::ceil(digits * 3.3219281)) + 64;
  int effective_digits = std::min(digits, n + 3);
  BigFloat threshold = BigFloat::from_string(
      "1e-" + std::to_string(effective_digits - 2), prec_bits);
  std::string threshold_str = threshold.to_decimal_string(2, true);

  // 调用通用验证函数
  run_csv_verification(csv_path, coeffs, n, g_str, digits, -1, threshold_str,
                       output_dir, true);

  std::cout << "\n=== Output Files ===" << std::endl;
  std::cout << "  " << output_dir
            << "/coefficients.txt   - Lanczos coefficients" << std::endl;
  std::cout << "  " << output_dir
            << "/parameters.txt     - Computation parameters" << std::endl;
  std::cout << "  " << output_dir
            << "/verification.txt   - Gamma function test results" << std::endl;

  return 0;
}
