# tests where program to build and program to run are equal
NORMALTESTS = test_diffusion test_diffusion3d test_mpfa \
			  test_equilibrated test_gradedmesh test_finesubdomainmesh \
			  test_cpgrid test_normals test_dgf 

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

test_diffusion_SOURCES = test_diffusion.cc

                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la 

                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la \
			-lblas -llapack

                    ../../../dune-cornerpoint/grid/preprocess/libpreprocess.la 

EXTRA_DIST=CMakeLists.txt
noinst_HEADERS=*.hh
include $(top_srcdir)/am/global-rules
