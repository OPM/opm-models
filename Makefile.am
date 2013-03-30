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
test/implicit/co2injection_immiscible \
test/implicit/co2injection_flash \
test/implicit/co2injection_ncp \
test/implicit/co2injection_pvs \
test/implicit/co2injection_immiscible_ni \
test/implicit/co2injection_flash_ni \
test/implicit/co2injection_ncp_ni \
test/implicit/co2injection_pvs_ni \
test/implicit/cuvette_pvs \
test/implicit/diffusion_ncp \
test/implicit/diffusion_pvs \
test/implicit/diffusion_flash \
test/implicit/finger_immiscible \
test/implicit/fracture_discretefracture \
test/implicit/groundwater_immiscible \
test/implicit/infiltration_pvs \
test/implicit/lens_immiscible \
test/implicit/lens_richards \
test/implicit/obstacle_immiscible \
test/implicit/obstacle_ncp \
test/implicit/obstacle_pvs \
test/implicit/outflow_pvs \
test/implicit/powerinjection_darcy \
test/implicit/powerinjection_forchheimer \
test/implicit/reservoir_blackoil \
test/implicit/test_navierstokes \
test/implicit/test_quadrature \
test/implicit/test_stokes \
test/implicit/test_stokes2c \
test/implicit/test_stokesni \
test/implicit/waterair_pvs_ni

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
test_implicit_co2injection_immiscible_SOURCES=test/implicit/co2injection_immiscible.cc
test_implicit_co2injection_flash_SOURCES=test/implicit/co2injection_flash.cc
test_implicit_co2injection_ncp_SOURCES=test/implicit/co2injection_ncp.cc
test_implicit_co2injection_pvs_SOURCES=test/implicit/co2injection_pvs.cc
test_implicit_co2injection_immiscible_ni_SOURCES=test/implicit/co2injection_immiscible_ni.cc
test_implicit_co2injection_flash_ni_SOURCES=test/implicit/co2injection_flash_ni.cc
test_implicit_co2injection_ncp_ni_SOURCES=test/implicit/co2injection_ncp_ni.cc
test_implicit_co2injection_pvs_ni_SOURCES=test/implicit/co2injection_pvs_ni.cc
test_implicit_cuvette_pvs_SOURCES=test/implicit/cuvette_pvs.cc
test_implicit_diffusion_flash_SOURCES=test/implicit/diffusion_flash.cc
test_implicit_diffusion_ncp_SOURCES=test/implicit/diffusion_ncp.cc
test_implicit_diffusion_pvs_SOURCES=test/implicit/diffusion_pvs.cc
test_implicit_finger_immiscible_SOURCES=test/implicit/finger_immiscible.cc
test_implicit_fracture_discretefracture_SOURCES=test/implicit/fracture_discretefracture.cc
test_implicit_groundwater_immiscible_SOURCES=test/implicit/groundwater_immiscible.cc 
test_implicit_infiltration_pvs_SOURCES=test/implicit/infiltration_pvs.cc
test_implicit_lens_immiscible_SOURCES=test/implicit/lens_immiscible.cc
test_implicit_lens_richards_SOURCES=test/implicit/lens_richards.cc
test_implicit_obstacle_immiscible_SOURCES=test/implicit/obstacle_immiscible.cc
test_implicit_obstacle_ncp_SOURCES=test/implicit/obstacle_ncp.cc
test_implicit_obstacle_pvs_SOURCES=test/implicit/obstacle_pvs.cc
test_implicit_outflow_pvs_SOURCES=test/implicit/outflow_pvs.cc
test_implicit_powerinjection_darcy_SOURCES=test/implicit/powerinjection_darcy.cc
test_implicit_powerinjection_forchheimer_SOURCES=test/implicit/powerinjection_forchheimer.cc
test_implicit_reservoir_blackoil_SOURCES=test/implicit/reservoir_blackoil.cc
test_implicit_test_navierstokes_SOURCES=test/implicit/test_navierstokes.cc
test_implicit_test_quadrature_SOURCES=test/implicit/test_quadrature.cc
test_implicit_test_stokes_SOURCES=test/implicit/test_stokes.cc
test_implicit_test_stokes2c_SOURCES=test/implicit/test_stokes2c.cc
test_implicit_test_stokesni_SOURCES=test/implicit/test_stokesni.cc
test_implicit_waterair_pvs_ni_SOURCES=test/implicit/waterair_pvs_ni.cc

# tutorials
tutorial_tutorial_coupled_SOURCES=tutorial/tutorial_coupled.cc
tutorial_tutorial_decoupled_SOURCES=tutorial/tutorial_decoupled.cc

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

# convenience shortcuts
co2injection_immiscible: test/implicit/co2injection_immiscible

co2injection_flash: test/implicit/co2injection_flash

co2injection_ncp: test/implicit/co2injection_ncp

co2injection_pvs: test/implicit/co2injection_pvs

co2injection_immiscible_ni: test/implicit/co2injection_immiscible_ni

co2injection_flash_ni: test/implicit/co2injection_flash_ni

co2injection_ncp_ni: test/implicit/co2injection_ncp_ni

co2injection_pvs_ni: test/implicit/co2injection_pvs_ni

cuvette_pvs: test/implicit/cuvette_pvs

diffusion_flash: test/implicit/diffusion_flash

diffusion_ncp: test/implicit/diffusion_ncp

diffusion_pvs: test/implicit/diffusion_pvs

finger_immiscible: test/implicit/finger_immiscible

fracture_discretefracture: test/implicit/fracture_discretefracture

groundwater_immiscible: test/implicit/groundwater_immiscible 

infiltration_pvs: test/implicit/infiltration_pvs

lens_immiscible: test/implicit/lens_immiscible

lens_richards: test/implicit/lens_richards

obstacle_immiscible: test/implicit/obstacle_immiscible

obstacle_ncp: test/implicit/obstacle_ncp

obstacle_pvs: test/implicit/obstacle_pvs

outflow_pvs: test/implicit/outflow_pvs

powerinjection_darcy: test/implicit/powerinjection_darcy

powerinjection_forchheimer: test/implicit/powerinjection_forchheimer

reservoir_blackoil: test/implicit/reservoir_blackoil

test_navierstokes: test/implicit/test_navierstokes

test_quadrature: test/implicit/test_quadrature

test_stokes: test/implicit/test_stokes

test_stokes2c: test/implicit/test_stokes2c

test_stokesni: test/implicit/test_stokesni

waterair_pvs_ni: test/implicit/waterair_pvs_ni

test_1p:test/decoupled/1p/test_1p

test_diffusion:test/decoupled/1p/test_diffusion

test_transport:test/decoupled/2p/test_transport

test_impes:test/decoupled/2p/test_impes

test_mpfa2p:test/decoupled/2p/test_mpfa2p

test_impesadaptive:test/decoupled/2p/test_impesadaptive

test_dec2p2c:test/decoupled/2p2c/test_dec2p2c

test_adaptive2p2c:test/decoupled/2p2c/test_adaptive2p2c

test_multiphysics2p2c:test/decoupled/2p2c/test_multiphysics2p2c

test_immiscibleflash:test/material/immiscibleflash/test_immiscibleflash

test_pengrobinson:test/material/pengrobinson/test_pengrobinson

test_tabulation:test/material/tabulation/test_tabulation

test_fluidsystems:test/material/fluidsystems/test_fluidsystems

test_ncpflash:test/material/ncpflash/test_ncpflash

test_spline:test/common/spline/test_spline

test_propertysystem:test/common/propertysystem/test_propertysystem

tutorial_coupled:tutorial/tutorial_coupled

tutorial_decoupled:tutorial/tutorial_decoupled

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/top-rules
