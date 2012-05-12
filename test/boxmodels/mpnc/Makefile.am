check_PROGRAMS = test_mpnc obstacle_2p2c
noinst_HEADERS = *.hh

EXTRA_DIST=*reference.vtu *.input CMakeLists.txt grids/*.dgf
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/*.dgf

test_mpnc_SOURCES = test_mpnc.cc
obstacle_2p2c_SOURCES = obstacle_2p2c.cc

include $(top_srcdir)/am/global-rules
