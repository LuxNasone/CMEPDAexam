# CMake generated Testfile for 
# Source directory: /home/lux_n/CMEPDA/Exam/Repo
# Build directory: /home/lux_n/CMEPDA/Exam/Repo/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(VarsTests "/home/lux_n/CMEPDA/Exam/Repo/build/test_vars")
set_tests_properties(VarsTests PROPERTIES  _BACKTRACE_TRIPLES "/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;29;add_test;/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;0;")
add_test(MainTests "/home/lux_n/CMEPDA/Exam/Repo/build/test_main")
set_tests_properties(MainTests PROPERTIES  _BACKTRACE_TRIPLES "/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;42;add_test;/home/lux_n/CMEPDA/Exam/Repo/CMakeLists.txt;0;")
subdirs("external/googletest")
