# CMake generated Testfile for 
# Source directory: /home/lux_n/CMEPDA/Exam/Repo
# Build directory: /home/lux_n/CMEPDA/Exam/Repo/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(VarsTests "/home/lux_n/CMEPDA/Exam/Repo/build/test_vars")
set_tests_properties(VarsTests PROPERTIES  _BACKTRACE_TRIPLES "/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;30;add_test;/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;0;")
add_test(AnalysisTests "/home/lux_n/CMEPDA/Exam/Repo/build/test_analysis")
set_tests_properties(AnalysisTests PROPERTIES  _BACKTRACE_TRIPLES "/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;45;add_test;/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;0;")
subdirs("external/googletest")
