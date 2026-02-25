#include "lanczos.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

// 解析 CSV 行: "z","gamma_value"
static bool parse_csv_line(const std::string& line, std::string& z_str, std::string& gamma_str) {
    // 格式: "value1","value2"
    if (line.empty() || line[0] != '"') return false;
    
    size_t first_end = line.find('"', 1);
    if (first_end == std::string::npos) return false;
    z_str = line.substr(1, first_end - 1);
    
    // 找到第二个引号对
    size_t second_start = line.find('"', first_end + 1);
    if (second_start == std::string::npos) return false;
    size_t second_end = line.find('"', second_start + 1);
    if (second_end == std::string::npos) {
        // 值可能很长, 取到行尾
        gamma_str = line.substr(second_start + 1);
    } else {
        gamma_str = line.substr(second_start + 1, second_end - second_start - 1);
    }
    return !z_str.empty() && !gamma_str.empty();
}

// 修正 CSV 中缺少前导 "0." 的 expected 值
static std::string fix_expected_value(const std::string& computed, const std::string& expected) {
    // 如果 expected 已包含小数点, 不需要修正
    if (expected.find('.') != std::string::npos) return expected;
    
    // 如果 expected 很短 (如 "1", "2", "24"), 它是精确整数, 不修正
    if (expected.size() <= 10) return expected;
    
    // 检查 computed 值是否 < 1 (以 "0." 开头或 "-0." 开头)
    bool computed_lt_one = false;
    if (computed.size() >= 2) {
        if (computed[0] == '0' && computed[1] == '.') computed_lt_one = true;
        if (computed[0] == '-' && computed.size() >= 3 && computed[1] == '0' && computed[2] == '.') computed_lt_one = true;
    }
    
    if (computed_lt_one) {
        // CSV 中的值省略了 "0.", 补上
        return "0." + expected;
    }
    return expected;
}

int main(int argc, char* argv[]) {
    if (argc < 4 && argc != 3) {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  Mode 1 (Generate): " << argv[0] << " <n> <g> <digits> [output_dir] [csv_path]" << std::endl;
        std::cerr << "  Mode 2 (Evaluate): " << argv[0] << " eval <output_dir> <z_value> [display_digits]" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Example:" << std::endl;
        std::cerr << "  " << argv[0] << " 7 5.0 16" << std::endl;
        std::cerr << "  " << argv[0] << " eval output_n7_g5_d16 50.5" << std::endl;
        std::cerr << "  " << argv[0] << " eval output_n7_g5_d16 50.5 50" << std::endl;
        return 1;
    }

    // ========== 模式 2: 从已有测试目录读取参数和系数计算 Gamma ==========
    std::string mode_flag = argv[1];
    if (mode_flag == "eval") {
        if (argc < 4) {
            std::cerr << "Error: 'eval' mode requires <output_dir> and <z_value>." << std::endl;
            return 1;
        }
        std::string output_dir = argv[2];
        std::string z_str = argv[3];
        
        int custom_display_digits = -1;
        if (argc >= 5) {
            custom_display_digits = std::atoi(argv[4]);
        }
        
        std::string params_path = output_dir + "/parameters.txt";
        std::string coeff_path = output_dir + "/coefficients.txt";
        
        std::ifstream fparams(params_path);
        if (!fparams.is_open()) {
            std::cerr << "Error: cannot open " << params_path << std::endl;
            return 1;
        }
        
        int n = 0, digits = 0;
        double g = 0.0;
        std::string g_str;
        std::string line;
        while (std::getline(fparams, line)) {
            if (line.find("n = ") == 0) n = std::stoi(line.substr(4));
            else if (line.find("g = ") == 0) {
                g_str = line.substr(4);
                g = std::stod(g_str);
            }
            else if (line.find("precision_decimal_digits = ") == 0) digits = std::stoi(line.substr(27));
        }
        fparams.close();
        
        if (n <= 0 || digits <= 0) {
            std::cerr << "Error: invalid parameters found in " << params_path << std::endl;
            return 1;
        }
        
        int display_digits = (custom_display_digits > 0) ? custom_display_digits : std::min(digits + 5, 40);
        
        int work_digits = std::max(digits, display_digits);
        int prec_bits = static_cast<int>(std::ceil(work_digits * 3.3219281)) + 64;
        
        std::vector<BigFloat> coeffs;
        
        std::ifstream fcoeffs(coeff_path);
        if (!fcoeffs.is_open()) {
            std::cerr << "Error: cannot open " << coeff_path << std::endl;
            return 1;
        }
        
        while (std::getline(fcoeffs, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t comma_pos = line.find(',');
            if (comma_pos != std::string::npos) {
                std::string val_str = line.substr(comma_pos + 1);
                size_t first_non_space = val_str.find_first_not_of(" \t");
                if (first_non_space != std::string::npos) {
                    val_str = val_str.substr(first_non_space);
                }
                coeffs.push_back(BigFloat::from_string(val_str, prec_bits));
            }
        }
        fcoeffs.close();
        
        if (coeffs.size() != static_cast<size_t>(n)) {
            std::cerr << "Error: expected " << n << " coefficients, but found " << coeffs.size() << std::endl;
            return 1;
        }
        
        BigFloat z = BigFloat::from_string(z_str, prec_bits);
        BigFloat gamma_val = lanczos_gamma(z, coeffs, g_str, work_digits);
        
        std::cout << "Gamma(" << z_str << ") = " << gamma_val.to_decimal_string(display_digits) << std::endl;
        
        return 0;
    }

    // ========== 模式 1: 生成并覆盖测试 ==========
    int n = std::atoi(argv[1]);
    std::string g_str = argv[2];
    double g = std::stod(g_str);
    int digits = std::atoi(argv[3]);

    if (n <= 0 || digits <= 0) {
        std::cerr << "Error: n and digits must be positive integers." << std::endl;
        return 1;
    }

    std::string output_dir;
    if (argc >= 5) {
        output_dir = argv[4];
    } else {
        std::ostringstream oss;
        oss << "output_n" << n << "_g" << g << "_d" << digits;
        output_dir = oss.str();
    }

    std::string csv_path = "real_gamma.csv";
    if (argc >= 6) {
        csv_path = argv[5];
    }

    fs::create_directories(output_dir);

    std::cout << "=== Lanczos Coefficient Calculator ===" << std::endl;
    std::cout << "Parameters: n=" << n << ", g=" << g << ", precision=" << digits << " decimal digits" << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << std::endl;

    auto coeffs = compute_lanczos_coefficients(n, g_str, digits);

    std::cout << std::endl;
    std::cout << "=== Computed Coefficients ===" << std::endl;
    for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
        std::cout << "  p[" << i << "] = " << coeffs[i].to_decimal_string(digits) << std::endl;
    }

    {
        std::string coeff_path = output_dir + "/coefficients.txt";
        std::ofstream fout(coeff_path);
        if (!fout.is_open()) {
            std::cerr << "Error: cannot open " << coeff_path << " for writing." << std::endl;
            return 1;
        }
        fout << "# Lanczos Approximation Coefficients" << std::endl;
        fout << "# n = " << n << std::endl;
        fout << "# g = " << g << std::endl;
        fout << "# precision = " << digits << " decimal digits" << std::endl;
        fout << "#" << std::endl;
        fout << "# Formula:" << std::endl;
        fout << "#   Gamma(z) = base^(z+0.5) * exp(-base) * S(z)" << std::endl;
        fout << "#   where base = z + g + 0.5" << std::endl;
        fout << "#         S(z) = p[0] + sum_{k=1}^{n-1} p[k] / (z + k)" << std::endl;
        fout << "#   Note: sqrt(2*pi) factor is absorbed into p[0]" << std::endl;
        fout << "#" << std::endl;
        fout << "# index, coefficient" << std::endl;
        for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
            fout << i << ", " << coeffs[i].to_decimal_string(digits) << std::endl;
        }
        fout.close();
        std::cout << "\nCoefficients written to: " << coeff_path << std::endl;
    }

    {
        std::string params_path = output_dir + "/parameters.txt";
        std::ofstream fout(params_path);
        fout << "n = " << n << std::endl;
        fout << "g = " << g << std::endl;
        fout << "precision_decimal_digits = " << digits << std::endl;
        int bits = static_cast<int>(std::ceil(digits * 3.3219281)) + 64;
        fout << "precision_binary_bits = " << bits << std::endl;
        fout.close();
    }

    std::cout << "\n=== Loading test data from " << csv_path << " ===" << std::endl;

    std::ifstream csv_file(csv_path);
    if (!csv_file.is_open()) {
        std::cerr << "Warning: cannot open " << csv_path << ", skipping CSV verification." << std::endl;
        std::cerr << "You can specify the CSV path as the 5th argument." << std::endl;
        return 0;
    }

    struct TestEntry {
        std::string z_str;
        std::string expected_str;
    };
    std::vector<TestEntry> test_data;

    std::string line;
    while (std::getline(csv_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        
        std::string z_str, gamma_str;
        if (parse_csv_line(line, z_str, gamma_str)) {
            test_data.push_back({z_str, gamma_str});
        }
    }
    csv_file.close();

    std::cout << "Loaded " << test_data.size() << " test entries." << std::endl;

    std::vector<TestEntry> selected;
    for (auto& entry : test_data) {
        double z_val = std::atof(entry.z_str.c_str());
        if (z_val > 0 && z_val <= 50) {
            selected.push_back(entry);
        }
    }
    std::cout << "Selected " << selected.size() << " test points (z <= 50)." << std::endl;

    std::cout << "\n=== Verification ===" << std::endl;

    std::string verify_path = output_dir + "/verification.txt";
    std::ofstream vout(verify_path);
    vout << "# Lanczos Gamma Function Verification Results" << std::endl;
    vout << "# Parameters: n=" << n << ", g=" << g << ", digits=" << digits << std::endl;
    vout << "# Test data source: " << csv_path << std::endl;
    vout << "#" << std::endl;
    vout << "# z, computed, expected (first 30 digits), absolute_error, matching_digits, status" << std::endl;

    int total_tests = static_cast<int>(selected.size());
    int pass_count = 0;
    int min_matching = 9999;
    int max_matching = 0;
    double sum_matching = 0;

    int prec_bits = static_cast<int>(std::ceil(digits * 3.3219281)) + 64;
    int display_digits = digits;

    for (int idx = 0; idx < total_tests; idx++) {
        auto& entry = selected[idx];

        BigFloat z = BigFloat::from_string(entry.z_str, prec_bits);
        BigFloat gamma_val = lanczos_gamma(z, coeffs, g_str, digits);
        std::string computed_str = gamma_val.to_decimal_string(display_digits);
        
        int expected_prec_bits = static_cast<int>(entry.expected_str.length() * 3.3219281) + 64;
        int compare_prec = std::max(prec_bits, expected_prec_bits);
        
        BigFloat gamma_val_ext = gamma_val;
        gamma_val_ext.set_precision(compare_prec);
        
        std::string fixed_expected = fix_expected_value(gamma_val.to_decimal_string(10, false), entry.expected_str);
        BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);
        
        BigFloat abs_err = (gamma_val_ext - expected_val).abs();
        std::string abs_err_str = abs_err.to_decimal_string(10, true);

        // 取期望值的前 30 位用于显示
        std::string expected_short = entry.expected_str.substr(0, std::min((size_t)30, entry.expected_str.size()));

        // 【纯字符串比对法】：完全摒弃 float/double，直接对比每一个有效位字符
        int matching = 0;
        if (expected_val.is_zero()) {
            matching = abs_err.is_zero() ? digits : 0;
        } else if (abs_err.is_zero()) {
            matching = digits;
        } else {
            // 强制将计算值与期望值都格式化为同等精度的科学记数法 (e.g. "X.YYYYe+ZZ")
            std::string c_str = gamma_val_ext.to_decimal_string(digits, true);
            std::string e_str = expected_val.to_decimal_string(digits, true);
            
            size_t c_e_pos = c_str.find('e');
            size_t e_e_pos = e_str.find('e');
            
            if (c_e_pos != std::string::npos && e_e_pos != std::string::npos) {
                // 提取指数部分
                std::string c_exp = c_str.substr(c_e_pos);
                std::string e_exp = e_str.substr(e_e_pos);
                
                // 只有在量级（指数）完全相同时，才对比有效数字
                if (c_exp == e_exp) {
                    size_t limit = std::min(c_e_pos, e_e_pos);
                    for (size_t i = 0; i < limit; i++) {
                        if (c_str[i] == e_str[i]) {
                            // 跳过小数点和负号，只统计数字
                            if (std::isdigit(c_str[i])) {
                                matching++;
                            }
                        } else {
                            // 遇到第一个不一样的位，彻底打断
                            break; 
                        }
                    }
                }
            }
        }

        std::string status = (matching >= 8) ? "PASS" : "FAIL";
        if (matching >= 8) pass_count++;

        if (matching < min_matching) min_matching = matching;
        if (matching > max_matching) max_matching = matching;
        sum_matching += matching;

        std::cout << "  Gamma(" << entry.z_str << ") = " << computed_str
                  << "  [Error: " << abs_err_str << ", " << matching << " digits, " << status << "]" << std::endl;

        vout << entry.z_str << ", " << computed_str
             << ", " << expected_short << "..."
             << ", " << abs_err_str
             << ", " << matching
             << ", " << status << std::endl;
    }

    double avg_matching = (total_tests > 0) ? (sum_matching / total_tests) : 0;

    vout << "#" << std::endl;
    vout << "# ==============================" << std::endl;
    vout << "# Summary" << std::endl;
    vout << "# ==============================" << std::endl;
    vout << "# Total tests:        " << total_tests << std::endl;
    vout << "# Passed (>= 8 digits): " << pass_count << "/" << total_tests << std::endl;
    vout << "# Min matching digits: " << min_matching << std::endl;
    vout << "# Max matching digits: " << max_matching << std::endl;
    vout << "# Avg matching digits: " << avg_matching << std::endl;
    vout.close();

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "  Total tests:         " << total_tests << std::endl;
    std::cout << "  Passed (>= 8 digits): " << pass_count << "/" << total_tests << std::endl;
    std::cout << "  Min matching digits: " << min_matching << std::endl;
    std::cout << "  Max matching digits: " << max_matching << std::endl;
    std::cout << "  Avg matching digits: " << avg_matching << std::endl;

    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "  " << output_dir << "/coefficients.txt   - Lanczos coefficients" << std::endl;
    std::cout << "  " << output_dir << "/parameters.txt     - Computation parameters" << std::endl;
    std::cout << "  " << output_dir << "/verification.txt   - Gamma function test results" << std::endl;

    return 0;
}