/**
 * @file test_lanczos.cpp
 * @brief Lanczos 近似精度验证测试
 *
 * @details
 * 本测试程序从 CSV 文件（real_gamma.csv）读取已知的 Γ(z) 高精度值，
 * 与 Lanczos 近似计算的结果进行逐位比对，统计匹配的有效数字位数。
 *
 * CSV 格式: "z","Gamma(z)"
 *   其中 Gamma(z) 是高精度十进制字符串
 *
 * 验证方法: 纯字符串比对法
 *   将计算值和期望值都格式化为相同精度的科学记数法，
 *   逐字符比较尾数部分来统计匹配位数。
 *   这种方法完全避免了 double 精度限制带来的误差。
 *
 * 命令行参数:
 *   test_lanczos [n] [g] [precision] [csv_path] [max_tests]
 *   默认: n=7, g=5.0, precision=16, csv_path=real_gamma.csv, max_tests=50
 *
 * 通过标准: 匹配 ≥ 8 位十进制数字为 PASS
 */

#include "lanczos.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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
    int n = 7; // Lanczos 级数项数
    std::string g_str = "5.0"; // Lanczos 参数 g
    int precision = 16; // 目标精度（十进制位数）
    std::string csv_path = "real_gamma.csv"; // CSV 验证文件路径
    int max_tests = 50; // 最多测试点数

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

    // === 显示参数信息 ===
    std::cout << "=== Lanczos Approximation Test ===" << std::endl;
    std::cout << "Parameters: n=" << n << ", g=" << g_str
            << ", precision=" << precision << " decimal digits" << std::endl;
    std::cout << "CSV file: " << csv_path << std::endl;
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
    std::cout << std::setw(12) << "z" << std::setw(15) << "matching_digits"
            << std::setw(15) << "status" << std::endl;
    std::cout << std::string(42, '-') << std::endl;

    // === 逐点验证 ===
    int total_tests = 0; // 总测试点数
    int passed = 0; // 通过数
    int failed = 0; // 失败数
    int total_matching = 0; // 匹配位数之和（用于计算平均值）

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

        // 跳过 z ≤ 0.5（当前实现不支持反射公式）
        if (z <= BigFloat(0.5, prec_bits))
            continue;

        // --- 计算 Γ(z) ---
        BigFloat computed = lanczos_gamma(z, coeffs, g_str, precision);

        // --- 准备比较: 对齐两个值的精度 ---
        int expected_prec_bits =
                static_cast<int>(expected_str.length() * 3.3219281) + 64;
        int compare_prec = std::max(prec_bits, expected_prec_bits);

        BigFloat computed_ext = computed;
        computed_ext.set_precision(compare_prec);

        // 修正 CSV 中可能缺失的 "0." 前缀
        std::string fixed_expected = expected_str;
        if (computed < BigFloat(1, prec_bits) && expected_str.length() > 5) {
            if (expected_str.substr(0, 2) != "0.") {
                if (expected_str[0] == '.') {
                    fixed_expected = "0" + expected_str; // ".886" → "0.886"
                } else {
                    fixed_expected = "0." + expected_str; // "886..." → "0.886..."
                }
            }
        }

        BigFloat expected_val = BigFloat::from_string(fixed_expected, compare_prec);
        BigFloat abs_err = (computed_ext - expected_val).abs();

        // --- 纯字符串比对法: 统计匹配的有效数字 ---
        // 原理: 将两个 BigFloat 格式化为相同精度的科学记数法，
        //       先比较指数部分，指数相同时逐字符比较尾数
        int matching = 0;
        if (expected_val.is_zero()) {
            matching = abs_err.is_zero() ? precision : 0;
        } else if (abs_err.is_zero()) {
            matching = precision; // 完全匹配
        } else {
            // 格式化为科学记数法字符串
            std::string c_str = computed_ext.to_decimal_string(precision, true);
            std::string e_str = expected_val.to_decimal_string(precision, true);

            // 定位 'e' 分隔符
            size_t c_e_pos = c_str.find('e');
            size_t e_e_pos = e_str.find('e');

            if (c_e_pos != std::string::npos && e_e_pos != std::string::npos) {
                // 比较指数部分（整数比较，避免 "e+0" vs "e+00" 问题）
                int c_exp_val = std::stoi(c_str.substr(c_e_pos + 1));
                int e_exp_val = std::stoi(e_str.substr(e_e_pos + 1));

                // 只在量级（指数）相同时比较有效数字
                if (c_exp_val == e_exp_val) {
                    size_t limit = std::min(c_e_pos, e_e_pos);
                    for (size_t i = 0; i < limit; i++) {
                        if (c_str[i] == e_str[i]) {
                            if (std::isdigit(c_str[i])) {
                                matching++; // 仅统计匹配的数字字符
                            }
                        } else {
                            break; // 遇到第一个不匹配的字符即停止
                        }
                    }
                }
            }
        }

        total_matching += matching;

        // --- 判定通过/失败 ---
        std::string status;
        if (matching >= 8) {
            status = "PASS";
            passed++;
        } else {
            status = "FAIL (" + std::to_string(matching) + "/" +
                     std::to_string(precision) + ")";
            failed++;
        }

        // 输出当前测试点结果
        std::cout << std::setw(12) << z_str << std::setw(15) << matching
                << std::setw(15) << status << std::endl;

        total_tests++;
    }

    // === 输出统计摘要 ===
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
