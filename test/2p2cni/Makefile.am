# tests where program to build and program to run are equal
NORMALTESTS = new_leakywell new_test_2p2cni test_2p2cni

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_2p2cnidir = $(EXTRA_DIST=CMakeLists.txt
includedir)/dune/dumux/test/2p2cni
test_2p2cni_HEADERS = brineco2problem.hh co2_soilproperties.hh waterairproblem.hh

new_test_2p2cni_SOURCES = new_test_2p2cni.cc

test_2p2cni_SOURCES = test_2p2cni.cc

EXTRA_DIST=CMakeLists.txt
include $(top_srcdir)/am/global-rules
