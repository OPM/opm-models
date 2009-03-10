# tests where program to build and program to run are equal
NORMALTESTS = test_twophase test_twophase_pnsw

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = test_twophase.cc test_twophase_pnsw.cc test_dmt.cc

test_twophase_SOURCES = test_twophase.cc

test_twophase_pnsw_SOURCES = test_twophase_pnsw.cc

EXTRA_DIST=CMakeLists.txt
noinst_HEADERS=*.hh
include $(top_srcdir)/am/global-rules
