# CMake generated Testfile for 
# Source directory: /home/roger/Projects/ksl
# Build directory: /home/roger/Projects/ksl/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(check_axis "/home/roger/Projects/ksl/src/test/check_axis")
add_test(check_coscrew "/home/roger/Projects/ksl/src/test/check_coscrew")
add_test(check_inertia "/home/roger/Projects/ksl/src/test/check_inertia")
add_test(check_linalg "/home/roger/Projects/ksl/src/test/check_linalg")
add_test(check_matrix "/home/roger/Projects/ksl/src/test/check_matrix")
add_test(check_quaternion "/home/roger/Projects/ksl/src/test/check_quaternion")
add_test(check_screw "/home/roger/Projects/ksl/src/test/check_screw")
add_test(check_util "/home/roger/Projects/ksl/src/test/check_util")
add_test(check_vector "/home/roger/Projects/ksl/src/test/check_vector")
add_test(check_print "/home/roger/Projects/ksl/src/test/check_print")
subdirs("test")
subdirs("doc")
