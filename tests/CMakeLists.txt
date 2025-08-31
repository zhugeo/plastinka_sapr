file(GLOB_RECURSE sources "*.cpp")
add_executable(PlastinkaTests ${sources})

target_link_libraries(PlastinkaTests GTest::gtest_main)
target_include_directories(PlastinkaTests PUBLIC ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ../include)

gtest_discover_tests(PlastinkaTests)
