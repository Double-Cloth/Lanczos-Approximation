#include "lanczos.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

/**
 * 测试程序: 使用 real_gamma.csv 验证 Lanczos 系数的精度
 * * CSV 格式: "z","Gamma(z)"
 * 其中 Gamma(z) 是高精度十进制字符串
 * * 对每个测试点计算误差大小和匹配的十进制数字个数
 */

// 解析 CSV 中的引号字符串
std::string parse_csv_field(const std::string& line, size_t& pos) {
    std::string result;
    if (pos >= line.size()) return result;
    
    if (line[pos] == '"') {
        pos++; // 跳过开头引号
        while (pos < line.size()) {
            if (line[pos] == '"') {
                if (pos + 1 < line.size() && line[pos + 1] == '"') {
                    result += '"';
                    pos += 2;
                } else {
                    pos++; // 跳过结尾引号
                    break;
                }
            } else {
                result += line[pos++];
            }
        }
    } else {
        while (pos < line.size() && line[pos] != ',') {
            result += line[pos++];
        }
    }
    if (pos < line.size() && line[pos] == ',') pos++; // 跳过逗号
    return result;
}

int main(int argc, char* argv[]) {
    // 默认参数
    int n = 7;
    std::string g_str = "5.0";
    int precision = 16;
    std::string csv_path = "real_gamma.csv";
    int max_tests = 50; // 最多测试多少个点

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

    std::cout << "=== Lanczos Approximation Test ===" << std::endl;
    std::cout << "Parameters: n=" << n << ", g=" << g_str 
              << ", precision=" << precision << " decimal digits" << std::endl;
    std::cout << "CSV file: " << csv_path << std::endl;
    std::cout << std::endl;

    // 计算系数
    std::cout << "--- Computing Lanczos Coefficients ---" << std::endl;
    auto coeffs = compute_lanczos_coefficients(n, g_str, precision);
    
    std::cout << std::endl << "Coefficients:" << std::endl;
    for (int i = 0; i < static_cast<int>(coeffs.size()); i++) {
        std::cout << "  p[" << i << "] = " << coeffs[i].to_decimal_string(precision + 5) << std::endl;
    }

    // 读取 CSV
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << csv_path << std::endl;
        return 1;
    }

    std::cout << std::endl << "--- Gamma Function Verification ---" << std::endl;
    std::cout << std::setw(12) << "z" 
              << std::setw(15) << "matching_digits"
              << std::setw(15) << "status"
              << std::endl;
    std::cout << std::string(42, '-') << std::endl;

    int total_tests = 0;
    int passed = 0;
    int failed = 0;
    int total_matching = 0;

    std::string line;
    while (std::getline(file, line) && total_tests < max_tests) {
        if (line.empty()) continue;
        
        // 解析 CSV: "z","gamma_value"
        size_t pos = 0;
        std::string z_str = parse_csv_field(line, pos);
        std::string expected_str = parse_csv_field(line, pos);
        
        if (z_str.empty() || expected_str.empty()) continue;
        
        // 解析 z (直接转换为高精度)
        int prec_bits = static_cast<int>(std::ceil(precision * 3.3219281)) + 64;
        BigFloat z;
        try {
            z = BigFloat::from_string(z_str, prec_bits);
        } catch (...) {
            continue;
        }
        
        // 跳过 z <= 0.5 (反射公式未实现)
        if (z <= BigFloat(0.5, prec_bits)) continue;
        
        // 计算 Gamma(z)
        BigFloat computed = lanczos_gamma(z, coeffs, g_str, precision);
        
        int expected_prec_bits = static_cast<int>(expected_str.length() * 3.3219281) + 64;
        int compare_prec = std::max(prec_bits, expected_prec_bits);
        
        BigFloat computed_ext = computed;
        computed_ext.set_precision(compare_prec);
        
        // 处理CSV缺省的 "0."
        std::string fixed_expected = expected_str;
        if (computed < BigFloat(1, prec_bits) && expected_str.length() > 5) {
            if (expected_str.substr(0, 2) != "0.") {
                if (expected_str[0] == '.') {
                    fixed_expected = "0" + expected_str;
                } else {
                    fixed_expected = "0." + expected_str;
                }
            }
        }
        
        BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);
        BigFloat abs_err = (computed_ext - expected_val).abs();
        
        // 【纯字符串比对法】：彻底摒弃 double 和指数估算
        int matching = 0;
        if (expected_val.is_zero()) {
            matching = abs_err.is_zero() ? precision : 0;
        } else if (abs_err.is_zero()) {
            matching = precision;
        } else {
            // 将两个 BigFloat 格式化为严格同样位长的科学记数法
            std::string c_str = computed_ext.to_decimal_string(precision, true);
            std::string e_str = expected_val.to_decimal_string(precision, true);
            
            size_t c_e_pos = c_str.find('e');
            size_t e_e_pos = e_str.find('e');
            
            if (c_e_pos != std::string::npos && e_e_pos != std::string::npos) {
                // 提取指数部分并转为整数进行比较（避免格式差异）
                int c_exp_val = std::stoi(c_str.substr(c_e_pos + 1));
                int e_exp_val = std::stoi(e_str.substr(e_e_pos + 1));
                
                // 严格对比尾数
                if (c_exp_val == e_exp_val) {
                    size_t limit = std::min(c_e_pos, e_e_pos);
                    for (size_t i = 0; i < limit; i++) {
                        if (c_str[i] == e_str[i]) {
                            if (std::isdigit(c_str[i])) {
                                matching++;
                            }
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        
        total_matching += matching;
        
        std::string status;
        if (matching >= 8) {
            status = "PASS";
            passed++;
        } else {
            status = "FAIL (" + std::to_string(matching) + "/" + std::to_string(precision) + ")";
            failed++;
        }
        
        std::cout << std::setw(12) << z_str
                  << std::setw(15) << matching
                  << std::setw(15) << status
                  << std::endl;
        
        total_tests++;
    }

    std::cout << std::endl << "=== Summary ===" << std::endl;
    std::cout << "Total tests: " << total_tests << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    if (total_tests > 0) {
        std::cout << "Average matching digits: " 
                  << static_cast<double>(total_matching) / total_tests << std::endl;
    }

    return failed > 0 ? 1 : 0;
}