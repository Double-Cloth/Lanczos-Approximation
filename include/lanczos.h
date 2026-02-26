/**
 * @file lanczos.h
 * @brief Lanczos 近似的系数计算与 Gamma 函数求值
 *
 * @details
 * Lanczos 近似是计算 Gamma 函数的高效方法。其核心公式表示为:
 *
 *   Γ(z) = ( (z + g - 0.5) / e )^(z-0.5) × S(z)
 *
 * 其中:
 *   S(z) = p₀ + Σ_{k=1}^{n-1} p_k / (z + k - 1)
 *
 * 系数 p₀, p₁, ..., p_{n-1} 通过 Godfrey 矩阵方法计算:
 *   P = D × B × C × F
 *
 * 其中:
 *   - D: 对角矩阵（递推生成）
 *   - B: 二项式系数矩阵
 *   - C: Chebyshev 系数矩阵
 *   - F: 包含 exp、sqrt、阶乘等的高精度浮点向量
 *
 * 参数说明:
 *   - n: 级数项数，越大精度越高（典型值 7~20）
 *   - g: Lanczos 参数，影响近似质量（典型值 5.0~10.5）
 *   - 精度: 十进制有效位数，内部转换为二进制位数
 *
 * 参考文献:
 *   - Paul Godfrey, "A note on the computation of the convergent
 *     Lanczos complex Gamma approximation"
 *   - Boost C++ Libraries, lanczos_generator.cpp
 *   - Viktor Toth, "Lanczos' Approximation for the Gamma Function"
 */
#pragma once
#include "BigFloat.h"
#include <vector>

/**
 * @brief 计算 Lanczos 近似系数
 *
 * 使用 Godfrey 矩阵方法: P = D × B × C × F
 *
 * @param n              级数项数（输出 n 个系数 p_0 ... p_{n-1}）
 * @param g_str          Lanczos 参数 g 的字符串形式（如 "5.0"），
 *                       直接转为 BigFloat 以保持全精度
 * @param decimal_digits 计算精度（十进制有效位数）
 * @return               系数列表 {p_0, p_1, ..., p_{n-1}}
 */
std::vector<BigFloat> compute_lanczos_coefficients(int n,
                                                   const std::string &g_str,
                                                   int decimal_digits);

/**
 * @brief 使用 Lanczos 系数估算 Γ(z)
 *
 * 计算公式:
 *   Γ(z) = ( (z+g-0.5) / e )^(z-0.5) × S(z)
 *
 * @param z              自变量（z > 0.5，因为未实现反射公式）
 * @param coeffs         由 compute_lanczos_coefficients 生成的系数列表
 * @param g_str          对应的 g 参数（字符串形式，必须与生成系数时一致）
 * @param decimal_digits 计算精度（十进制位数）
 * @return               Γ(z) 的高精度近似值
 *
 * @note 当前实现要求 z > 0.5。对于 z ≤ 0.5，需要使用反射公式:
 *       Γ(z) × Γ(1−z) = π / sin(πz)
 */
BigFloat lanczos_gamma(const BigFloat &z, const std::vector<BigFloat> &coeffs,
                       const std::string &g_str, int decimal_digits);
