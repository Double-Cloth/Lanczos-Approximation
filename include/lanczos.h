#pragma once
#include "BigFloat.h"
#include <vector>

/**
 * 计算 Lanczos 近似系数
 * 
 * 使用 Godfrey 矩阵方法: P = D * B * C * F
 * 
 * @param n    级数项数 (输出 n 个系数 p_0 ... p_{n-1})
 * @param g    Lanczos 参数
 * @param bits 计算精度 (有效位数, 十进制)
 * @return     系数列表 {p_0, p_1, ..., p_{n-1}}
 */
std::vector<BigFloat> compute_lanczos_coefficients(int n, const std::string& g_str, int bits);

/**
 * 使用计算出的 Lanczos 系数估算 Gamma(z)
 *
 * @param z      自变量 (z > 0)
 * @param coeffs 系数列表
 * @param g_str  对应的 g 参数 (字符串形式)
 * @param decimal_digits   计算精度 (十进制位数)
 * @return       Gamma(z) 的近似值
 */
BigFloat lanczos_gamma(const BigFloat& z, const std::vector<BigFloat>& coeffs, const std::string& g_str, int decimal_digits);
