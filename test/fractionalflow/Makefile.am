dist_noinst_DATA = test_fractionalflow.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_fractionalflow test_decoupled2p2c test_heterogeneoussoil

test_fractionalflowdir = $(EXTRA_DIST=CMakeLists.txt
includedir)/dune/dumux/test/fractionalflow
test_fractionalflow_HEADERS = heterogeneoussoil_testproblem.hh soilproperties.hh
test_fractionalflow_SOURCES = test_fractionalflow.cc

test_decoupled2p2c_SOURCES = test_decoupled2p2c.cc

EXTRA_DIST=CMakeLists.txt
include $(top_srcdir)/am/global-rules
