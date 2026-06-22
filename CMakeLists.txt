cmake_minimum_required(VERSION 3.16)

project(GeneticAlgos LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(sat
    problems/sat/main.cpp
)

add_executable(cos_function
    problems/cos-function/main.cpp
)

add_executable(radios
    problems/radios/main.cpp
)

add_executable(queens
    problems/queens/main.cpp
)

add_executable(weighted_queens
    problems/weighted-queens/main.cpp
)

add_executable(delivery
    problems/delivery/main.cpp
)

target_compile_options(sat PRIVATE
    -Wall
    -O2
    -Wextra
    -Wpedantic
)
