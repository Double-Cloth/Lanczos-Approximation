cmake_minimum_required(VERSION 3.15)
project(LanczosApproximation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(bignum STATIC
    src/BigInt.cpp
    src/BigFloat.cpp
)
target_include_directories(bignum PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_library(lanczos STATIC
    src/lanczos.cpp
)
target_include_directories(lanczos PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lanczos PUBLIC bignum)

add_executable(lanczos_app src/main.cpp)
target_link_libraries(lanczos_app PRIVATE lanczos)

add_executable(test_bigfloat tests/test_bigfloat.cpp)
target_link_libraries(test_bigfloat PRIVATE bignum)

add_executable(test_lanczos tests/test_lanczos.cpp)
target_link_libraries(test_lanczos PRIVATE lanczos)
