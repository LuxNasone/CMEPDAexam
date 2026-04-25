cmake_minimum_required(VERSION 3.10)
project(Cross Section for Z production in p-p collision)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)

# GoogleTest
add_subdirectory(googletest)

include_directories(
    ${ROOT_INCLUDE_DIRS}
    src
)

enable_testing()

add_executable(test_vars
    tests/test_vars.cpp
)

target_link_libraries(test_vars
    gtest
    gtest_main
    ${ROOT_LIBRARIES}
)

add_test(NAME VarsTests COMMAND test_vars)