cmake_minimum_required(VERSION 3.18)

project(PlastinkaSapr VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB sources "src/*.cpp")


add_executable(Plastinka ${sources})
target_include_directories(Plastinka PUBLIC include 3rdparty/eigen)
