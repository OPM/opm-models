# don't follow the full GNU-standard
# we need automake 1.9 or later
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS= CXX="$(CXX)" CC="$(CC)"

# to the tarball, add all files required to make dunecontrol work
# nicely, the files for the CMake build system, some README files and
# a set of default .opts files for dunecontrol.
EXTRA_DIST = dune.module \
	     .gitignore \
             cmake/Modules/*.cmake \
             bin/runtest.sh bin/fuzzycomparevtu.py \
             CMakeLists.txt CTestConfig.cmake config.h.cmake \
             CHANGELOG COPYING FAQ INSTALL.cmake README \
             EwomsConfig.cmake.in \
             *.opts 

SUBDIRS = doc ewoms m4 test tutorial

VCFV_TESTS= \
co2injection_immiscible \
co2injection_flash \
co2injection_ncp \
co2injection_pvs \
co2injection_immiscible_ni \
co2injection_flash_ni \
co2injection_ncp_ni \
co2injection_pvs_ni \
cuvette_pvs \
diffusion_ncp \
diffusion_pvs \
diffusion_flash \
finger_immiscible \
groundwater_immiscible \
infiltration_pvs \
lens_immiscible \
lens_richards \
obstacle_immiscible \
obstacle_ncp \
obstacle_pvs \
outflow_pvs \
powerinjection_darcy \
powerinjection_forchheimer \
reservoir_blackoil \
test_navierstokes \
test_quadrature \
test_stokes \
test_stokes2c \
test_stokesni \
waterair_pvs_ni

IMPET_TESTS= \
test_1p \
test_diffusion \
test_transport \
test_impes \
test_mpfa2p \
test_impesadaptive \
test_dec2p2c \
test_adaptive2p2c \
test_multiphysics2p2c

OTHER_TESTS= \
test_immiscibleflash \
test_pengrobinson \
test_tabulation \
test_fluidsystems \
test_ncpflash \
test_spline \
test_propertysystem

TUTORIALS= \
tutorial_coupled \
tutorial_decoupled

check_PROGRAMS = $(VCFV_TESTS) $(IMPET_TESTS) $(OTHER_TESTS) $(TUTORIALS)

# VCFV tests
co2injection_immiscible_SOURCES=test/boxmodels/co2injection_immiscible.cc
co2injection_flash_SOURCES=test/boxmodels/co2injection_flash.cc
co2injection_ncp_SOURCES=test/boxmodels/co2injection_ncp.cc
co2injection_pvs_SOURCES=test/boxmodels/co2injection_pvs.cc
co2injection_immiscible_ni_SOURCES=test/boxmodels/co2injection_immiscible_ni.cc
co2injection_flash_ni_SOURCES=test/boxmodels/co2injection_flash_ni.cc
co2injection_ncp_ni_SOURCES=test/boxmodels/co2injection_ncp_ni.cc
co2injection_pvs_ni_SOURCES=test/boxmodels/co2injection_pvs_ni.cc
cuvette_pvs_SOURCES=test/boxmodels/cuvette_pvs.cc
diffusion_flash_SOURCES=test/boxmodels/diffusion_flash.cc
diffusion_ncp_SOURCES=test/boxmodels/diffusion_ncp.cc
diffusion_pvs_SOURCES=test/boxmodels/diffusion_pvs.cc
finger_immiscible_SOURCES=test/boxmodels/finger_immiscible.cc
groundwater_immiscible_SOURCES=test/boxmodels/groundwater_immiscible.cc 
infiltration_pvs_SOURCES=test/boxmodels/infiltration_pvs.cc
lens_immiscible_SOURCES=test/boxmodels/lens_immiscible.cc
lens_richards_SOURCES=test/boxmodels/lens_richards.cc
obstacle_immiscible_SOURCES=test/boxmodels/obstacle_immiscible.cc
obstacle_ncp_SOURCES=test/boxmodels/obstacle_ncp.cc
obstacle_pvs_SOURCES=test/boxmodels/obstacle_pvs.cc
outflow_pvs_SOURCES=test/boxmodels/outflow_pvs.cc
powerinjection_darcy_SOURCES=test/boxmodels/powerinjection_darcy.cc
powerinjection_forchheimer_SOURCES=test/boxmodels/powerinjection_forchheimer.cc
reservoir_blackoil_SOURCES=test/boxmodels/reservoir_blackoil.cc
test_navierstokes_SOURCES=test/boxmodels/test_navierstokes.cc
test_quadrature_SOURCES=test/boxmodels/test_quadrature.cc
test_stokes_SOURCES=test/boxmodels/test_stokes.cc
test_stokes2c_SOURCES=test/boxmodels/test_stokes2c.cc
test_stokesni_SOURCES=test/boxmodels/test_stokesni.cc
waterair_pvs_ni_SOURCES=test/boxmodels/waterair_pvs_ni.cc

# IMPET tests
test_1p_SOURCES=test/decoupled/1p/test_1p.cc
test_diffusion_SOURCES=test/decoupled/1p/test_diffusion.cc
test_transport_SOURCES=test/decoupled/2p/test_transport.cc
test_impes_SOURCES=test/decoupled/2p/test_impes.cc
test_mpfa2p_SOURCES=test/decoupled/2p/test_mpfa2p.cc
test_impesadaptive_SOURCES=test/decoupled/2p/test_impesadaptive.cc
test_dec2p2c_SOURCES=test/decoupled/2p2c/test_dec2p2c.cc
test_adaptive2p2c_SOURCES=test/decoupled/2p2c/test_adaptive2p2c.cc
test_multiphysics2p2c_SOURCES=test/decoupled/2p2c/test_multiphysics2p2c.cc

# other tests
test_immiscibleflash_SOURCES=test/material/immiscibleflash/test_immiscibleflash.cc
test_pengrobinson_SOURCES=test/material/pengrobinson/test_pengrobinson.cc
test_tabulation_SOURCES=test/material/tabulation/test_tabulation.cc
test_fluidsystems_SOURCES=test/material/fluidsystems/test_fluidsystems.cc
test_ncpflash_SOURCES=test/material/ncpflash/test_ncpflash.cc
test_spline_SOURCES=test/common/spline/test_spline.cc
test_propertysystem_SOURCES=test/common/propertysystem/test_propertysystem.cc

# tutorials
tutorial_coupled_SOURCES=tutorial/tutorial_coupled.cc
tutorial_decoupled_SOURCES=tutorial/tutorial_decoupled.cc

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

