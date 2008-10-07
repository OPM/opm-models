dist_noinst_DATA = test_fractionalflow.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_fractionalflow test_decoupled2p2c test_heterogeneoussoil

test_fractionalflow_SOURCES = test_fractionalflow.cc
test_fractionalflow_DEPENDENCIES = $(LOCAL_LIBS)

test_decoupled2p2c_SOURCES = test_decoupled2p2c.cc
test_decoupled2p2c_DEPENDENCIES = $(LOCAL_LIBS)

EXTRA_DIST=CMakeLists.txt
noinst_HEADERS=*.hh
include $(top_srcdir)/am/global-rules
