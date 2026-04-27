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
)

enable_testing()

add_executable(test_vars
    tests/test_vars.cpp
    src/Vars.cpp
)

target_link_libraries(test_vars
    gtest
    gtest_main
    ${ROOT_LIBRARIES}
)

add_test(NAME VarsTests COMMAND test_vars)

add_executable(test_analysis
    tests/test_analysis.cpp
    src/Analysis.cpp
    src/
)

target_link_libraries(test_analysis
    gtest
    gtest_main
    ${ROOT_LIBRARIES}
    /home/lux_n/CMEPDA/Exam/RooUnfold/build/libRooUnfold.so
)

add_test(NAME AnalysisTests COMMAND test_analysis)