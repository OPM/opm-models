SUBDIRS = problems grids quadrature

EXTRA_DIST = CMakeLists.txt

check_PROGRAMS = \
	test_1p \
	test_2p \
	test_2pni \
	test_1p2c \
	test_2p2c \
	test_2p2cni \
	test_3p3c \
	test_3p3cni \
	obstacle_ncp \
	obstacle_pvs \
	test_blackoil \
	test_richards

test_1p_SOURCES=test_1p.cc
test_2p_SOURCES=test_2p.cc
test_2pni_SOURCES=test_2pni.cc
test_1p2c_SOURCES=test_1p2c.cc
test_2p2c_SOURCES=test_2p2c.cc
test_2p2cni_SOURCES=test_2p2cni.cc
test_3p3c_SOURCES=test_3p3c.cc
test_3p3cni_SOURCES=test_3p3cni.cc
obstacle_ncp_SOURCES=obstacle_ncp.cc
obstacle_pvs_SOURCES=obstacle_pvs.cc
test_blackoil_SOURCES=test_blackoil.cc
test_richards_SOURCES=test_richards.cc

include $(top_srcdir)/am/global-rules

