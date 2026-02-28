# ==============================================================================
# Lanczos Approximation - CMake 构建配置
# ==============================================================================
#
# 项目结构:
#   include/          - 头文件 (BigInt.h, BigFloat.h, lanczos.h)
#   src/              - 源文件
#   tests/            - 测试文件
#
# 构建目标:
#   bignum            - 静态库: BigInt + BigFloat (任意精度算术)
#   lanczos           - 静态库: Lanczos 系数计算与 Gamma 函数
#   lanczos_app       - 可执行文件: 命令行主程序
#   test_bigfloat     - 可执行文件: BigInt/BigFloat 单元测试
#   test_lanczos      - 可执行文件: Lanczos 近似精度验证测试
#
# 构建命令:
#   mkdir build && cd build
#   cmake .. -G "MinGW Makefiles"   (Windows + MinGW)
#   cmake --build .
# ==============================================================================

cmake_minimum_required(VERSION 3.15)
project(LanczosApproximation LANGUAGES CXX)

# 使用 C++17 标准（需要 std::filesystem, 结构化绑定等特性）
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- bignum 静态库: 任意精度算术基础设施 ---
# 包含 BigInt (大整数) 和 BigFloat (高精度浮点数) 的实现
add_library(bignum STATIC
    src/BigInt.cpp
    src/BigFloat.cpp
)
target_include_directories(bignum PUBLIC ${CMAKE_SOURCE_DIR}/include)

# --- lanczos 静态库: Lanczos 近似算法 ---
# 依赖 bignum 库进行高精度计算
add_library(lanczos STATIC
    src/lanczos.cpp
)
target_include_directories(lanczos PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lanczos PUBLIC bignum)

# --- 主程序: 系数计算与 Gamma 函数求值 ---
add_executable(lanczos_app src/main.cpp)
target_link_libraries(lanczos_app PRIVATE lanczos)

# --- 测试: BigInt/BigFloat 单元测试 ---
add_executable(test_bigfloat tests/test_bigfloat.cpp)
target_link_libraries(test_bigfloat PRIVATE bignum)
