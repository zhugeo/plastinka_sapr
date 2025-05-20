cmake_minimum_required(VERSION 3.18)

project(PlastinkaSapr VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# TODO refactor
file(GLOB sources1 "src/**/*.cpp")
file(GLOB sources2 "src/*.cpp")

add_executable(Plastinka ${sources1} ${sources2})
target_include_directories(Plastinka PUBLIC include 3rdparty/eigen)
