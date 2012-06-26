check_PROGRAMS = test_mpnc obstacle_pvs
noinst_HEADERS = *.hh

EXTRA_DIST=*reference.vtu *.input CMakeLists.txt grids/*.dgf
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/*.dgf

test_mpnc_SOURCES = test_mpnc.cc
obstacle_pvs_SOURCES = obstacle_pvs.cc

include $(top_srcdir)/am/global-rules
