check_PROGRAMS = obstacle_ncp obstacle_pvs
noinst_HEADERS = *.hh

EXTRA_DIST=*reference.vtu *.input CMakeLists.txt grids/*.dgf
gridsdir=$(datadir)/dumux/grids
grids_DATA=grids/*.dgf

obstacle_ncp_SOURCES = obstacle_ncp.cc
obstacle_pvs_SOURCES = obstacle_pvs.cc

include $(top_srcdir)/am/global-rules
