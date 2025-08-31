cmake_minimum_required(VERSION 3.21)

project(PlastinkaSapr VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB_RECURSE sources "src/*.cpp")
add_executable(Plastinka ${sources})

include(cmake/CPM.cmake)
CPMAddPackage(
    URI "gh:eigen-mirror/eigen#master"
    OPTIONS "EIGEN_BUILD_TESTING OFF"
)
target_link_libraries(Plastinka PUBLIC Eigen3::Eigen)

option(PLASTINKA_SAPR_BUILD_TESTS "Build Unit Tests for plastinka-sapr" OFF)
if(PLASTINKA_SAPR_BUILD_TESTS)
    enable_testing()
    include(GoogleTest)
    CPMAddPackage("gh:google/googletest#v1.17.0")
    add_subdirectory(tests)
endif()

target_include_directories(Plastinka PUBLIC include)
