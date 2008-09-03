dist_noinst_DATA = test_fractionalflow.cc

AM_LDFLAGS = @AM_LDFLAGS@ $(LOCAL_LIBS)

noinst_PROGRAMS = test_fractionalflow test_2p2c_decoupled test_heterogeneoussoil
test_fractionalflow_SOURCES = test_fractionalflow.cc
test_fractionalflow_DEPENDENCIES = $(LOCAL_LIBS)

noinst_HEADERS=*.hh
include $(top_srcdir)/am/global-rules
