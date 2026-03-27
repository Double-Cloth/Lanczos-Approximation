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

# ==============================================================================
# 1. 强制静态链接配置 (C/C++ 运行时库及系统库)
# ==============================================================================
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # 针对 GCC/MinGW/Clang: 全局添加静态链接标志
    # -static: 尽可能使用静态库链接所有依赖
    # -static-libgcc -static-libstdc++: 强制静态链接 C/C++ 核心运行库
    add_link_options(-static -static-libgcc -static-libstdc++)
elseif(MSVC)
    # 针对 MSVC: 强制使用静态多线程运行时库 (/MT 或 /MTd)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
# ==============================================================================

find_package(OpenMP)

# --- bignum 静态库: 任意精度算术基础设施 ---
add_library(bignum STATIC
    src/BigInt.cpp
    src/BigFloat.cpp
)
target_include_directories(bignum PUBLIC ${CMAKE_SOURCE_DIR}/include)

# --- lanczos 静态库: Lanczos 近似算法 ---
add_library(lanczos STATIC
    src/lanczos.cpp
)
target_include_directories(lanczos PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lanczos PUBLIC bignum)

# ==============================================================================
# 2. 处理 OpenMP 的静态链接 (避免引入 libgomp-1.dll)
# ==============================================================================
if(OpenMP_CXX_FOUND)
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # 避开 CMake FindOpenMP 可能找到的动态导入库 (.dll.a)
        # 直接传递编译和链接参数，配合全局的 -static 实现纯静态链接
        target_compile_options(bignum PUBLIC -fopenmp)
        target_link_options(bignum PUBLIC -fopenmp)
        
        target_compile_options(lanczos PUBLIC -fopenmp)
        target_link_options(lanczos PUBLIC -fopenmp)
    else()
        # MSVC 等其他编译器正常使用 CMake 的目标链接
        target_link_libraries(bignum PUBLIC OpenMP::OpenMP_CXX)
        target_link_libraries(lanczos PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()
# ==============================================================================

# --- 主程序: 系数计算与 Gamma 函数求值 ---
add_executable(lanczos_app src/main.cpp)
target_link_libraries(lanczos_app PRIVATE lanczos)
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # 再次在最终目标上强调静态链接，确保万无一失
    target_link_options(lanczos_app PRIVATE -static -static-libgcc -static-libstdc++)
endif()

# --- 测试: BigInt/BigFloat 单元测试 ---
add_executable(test_bigfloat tests/test_bigfloat.cpp)
target_link_libraries(test_bigfloat PRIVATE bignum)
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_options(test_bigfloat PRIVATE -static -static-libgcc -static-libstdc++)
endif()

enable_testing()
add_test(NAME test_bigfloat COMMAND test_bigfloat)