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
test/boxmodels/co2injection_immiscible \
test/boxmodels/co2injection_flash \
test/boxmodels/co2injection_ncp \
test/boxmodels/co2injection_pvs \
test/boxmodels/co2injection_immiscible_ni \
test/boxmodels/co2injection_flash_ni \
test/boxmodels/co2injection_ncp_ni \
test/boxmodels/co2injection_pvs_ni \
test/boxmodels/cuvette_pvs \
test/boxmodels/diffusion_ncp \
test/boxmodels/diffusion_pvs \
test/boxmodels/diffusion_flash \
test/boxmodels/finger_immiscible \
test/boxmodels/groundwater_immiscible \
test/boxmodels/infiltration_pvs \
test/boxmodels/lens_immiscible \
test/boxmodels/lens_richards \
test/boxmodels/obstacle_immiscible \
test/boxmodels/obstacle_ncp \
test/boxmodels/obstacle_pvs \
test/boxmodels/outflow_pvs \
test/boxmodels/powerinjection_darcy \
test/boxmodels/powerinjection_forchheimer \
test/boxmodels/reservoir_blackoil \
test/boxmodels/test_navierstokes \
test/boxmodels/test_quadrature \
test/boxmodels/test_stokes \
test/boxmodels/test_stokes2c \
test/boxmodels/test_stokesni \
test/boxmodels/waterair_pvs_ni

IMPET_TESTS= \
test/decoupled/1p/test_1p \
test/decoupled/1p/test_diffusion \
test/decoupled/2p/test_transport \
test/decoupled/2p/test_impes \
test/decoupled/2p/test_mpfa2p \
test/decoupled/2p/test_impesadaptive \
test/decoupled/2p2c/test_dec2p2c \
test/decoupled/2p2c/test_adaptive2p2c \
test/decoupled/2p2c/test_multiphysics2p2c

OTHER_TESTS= \
test/material/immiscibleflash/test_immiscibleflash \
test/material/pengrobinson/test_pengrobinson \
test/material/tabulation/test_tabulation \
test/material/fluidsystems/test_fluidsystems \
test/material/ncpflash/test_ncpflash \
test/common/spline/test_spline \
test/common/propertysystem/test_propertysystem

TUTORIALS= \
tutorial/tutorial_coupled \
tutorial/tutorial_decoupled

check_PROGRAMS = $(VCFV_TESTS) $(IMPET_TESTS) $(OTHER_TESTS) $(TUTORIALS)

# VCFV tests
test_boxmodels_co2injection_immiscible_SOURCES=test/boxmodels/co2injection_immiscible.cc
test_boxmodels_co2injection_flash_SOURCES=test/boxmodels/co2injection_flash.cc
test_boxmodels_co2injection_ncp_SOURCES=test/boxmodels/co2injection_ncp.cc
test_boxmodels_co2injection_pvs_SOURCES=test/boxmodels/co2injection_pvs.cc
test_boxmodels_co2injection_immiscible_ni_SOURCES=test/boxmodels/co2injection_immiscible_ni.cc
test_boxmodels_co2injection_flash_ni_SOURCES=test/boxmodels/co2injection_flash_ni.cc
test_boxmodels_co2injection_ncp_ni_SOURCES=test/boxmodels/co2injection_ncp_ni.cc
test_boxmodels_co2injection_pvs_ni_SOURCES=test/boxmodels/co2injection_pvs_ni.cc
test_boxmodels_cuvette_pvs_SOURCES=test/boxmodels/cuvette_pvs.cc
test_boxmodels_diffusion_flash_SOURCES=test/boxmodels/diffusion_flash.cc
test_boxmodels_diffusion_ncp_SOURCES=test/boxmodels/diffusion_ncp.cc
test_boxmodels_diffusion_pvs_SOURCES=test/boxmodels/diffusion_pvs.cc
test_boxmodels_finger_immiscible_SOURCES=test/boxmodels/finger_immiscible.cc
test_boxmodels_groundwater_immiscible_SOURCES=test/boxmodels/groundwater_immiscible.cc 
test_boxmodels_infiltration_pvs_SOURCES=test/boxmodels/infiltration_pvs.cc
test_boxmodels_lens_immiscible_SOURCES=test/boxmodels/lens_immiscible.cc
test_boxmodels_lens_richards_SOURCES=test/boxmodels/lens_richards.cc
test_boxmodels_obstacle_immiscible_SOURCES=test/boxmodels/obstacle_immiscible.cc
test_boxmodels_obstacle_ncp_SOURCES=test/boxmodels/obstacle_ncp.cc
test_boxmodels_obstacle_pvs_SOURCES=test/boxmodels/obstacle_pvs.cc
test_boxmodels_outflow_pvs_SOURCES=test/boxmodels/outflow_pvs.cc
test_boxmodels_powerinjection_darcy_SOURCES=test/boxmodels/powerinjection_darcy.cc
test_boxmodels_powerinjection_forchheimer_SOURCES=test/boxmodels/powerinjection_forchheimer.cc
test_boxmodels_reservoir_blackoil_SOURCES=test/boxmodels/reservoir_blackoil.cc
test_boxmodels_test_navierstokes_SOURCES=test/boxmodels/test_navierstokes.cc
test_boxmodels_test_quadrature_SOURCES=test/boxmodels/test_quadrature.cc
test_boxmodels_test_stokes_SOURCES=test/boxmodels/test_stokes.cc
test_boxmodels_test_stokes2c_SOURCES=test/boxmodels/test_stokes2c.cc
test_boxmodels_test_stokesni_SOURCES=test/boxmodels/test_stokesni.cc
test_boxmodels_waterair_pvs_ni_SOURCES=test/boxmodels/waterair_pvs_ni.cc

# IMPET tests
test_decoupled_1p_test_1p_SOURCES=test/decoupled/1p/test_1p.cc
test_decoupled_1p_test_diffusion_SOURCES=test/decoupled/1p/test_diffusion.cc
test_decoupled_2p_test_transport_SOURCES=test/decoupled/2p/test_transport.cc
test_decoupled_2p_test_impes_SOURCES=test/decoupled/2p/test_impes.cc
test_decoupled_2p_test_mpfa2p_SOURCES=test/decoupled/2p/test_mpfa2p.cc
test_decoupled_2p_test_impesadaptive_SOURCES=test/decoupled/2p/test_impesadaptive.cc
test_decoupled_2p2c_test_dec2p2c_SOURCES=test/decoupled/2p2c/test_dec2p2c.cc
test_decoupled_2p2c_test_adaptive2p2c_SOURCES=test/decoupled/2p2c/test_adaptive2p2c.cc
test_decoupled_2p2c_test_multiphysics2p2c_SOURCES=test/decoupled/2p2c/test_multiphysics2p2c.cc

# other tests
test_material_immiscibleflash_test_immiscibleflash_SOURCES=test/material/immiscibleflash/test_immiscibleflash.cc
test_material_pengrobinson_test_pengrobinson_SOURCES=test/material/pengrobinson/test_pengrobinson.cc
test_material_tabulation_test_tabulation_SOURCES=test/material/tabulation/test_tabulation.cc
test_material_fluidsystems_test_fluidsystems_SOURCES=test/material/fluidsystems/test_fluidsystems.cc
test_material_ncpflash_test_ncpflash_SOURCES=test/material/ncpflash/test_ncpflash.cc
test_common_spline_test_spline_SOURCES=test/common/spline/test_spline.cc
test_common_propertysystem_test_propertysystem_SOURCES=test/common/propertysystem/test_propertysystem.cc

# tutorials
tutorial_tutorial_coupled_SOURCES=tutorial/tutorial_coupled.cc
tutorial_tutorial_decoupled_SOURCES=tutorial/tutorial_decoupled.cc

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

