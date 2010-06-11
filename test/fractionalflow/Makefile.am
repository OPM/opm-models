# tests where program to build and program to run are equal
NORMALTESTS = test_fractionalflow test_decoupled2p2c

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_fractionalflowdir = $(EXTRA_DIST=CMakeLists.txt
noinst_HEADERS=*.hh
includedir)/dune/dumux/test/fractionalflow
test_fractionalflow_SOURCES = test_fractionalflow.cc

EXTRA_DIST=CMakeLists.txt
noinst_HEADERS=*.hh
include $(top_srcdir)/am/global-rules
