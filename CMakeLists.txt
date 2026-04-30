cmake_minimum_required(VERSION 3.10)
project(ZAnalysis)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)

add_subdirectory(external/googletest)

include_directories(
    ${ROOT_INCLUDE_DIRS}
    external/googletest/googletest/include
    /home/lux_n/CMEPDA/Exam/RooUnfold/build
    src
    include
)

enable_testing()

add_executable(test_utils
    tests/test_utils.cpp
    src/Utils.cpp
)

target_link_libraries(test_utils
    gtest
    gtest_main
    ${ROOT_LIBRARIES}
)

add_test(NAME VarsTests COMMAND test_vars)
