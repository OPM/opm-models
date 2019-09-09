# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.
#
# ATTIC_FILES           Unmaintained files. This for the projects developers
#                       only. Don't expect these files to build.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
#list(APPEND MAIN_SOURCE_FILES "")

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
#list(APPEND TEST_SOURCE_FILES "")

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
file(GLOB_RECURSE TMP_GRIDS RELATIVE "${PROJECT_SOURCE_DIR}" "tests/*.dgf")
file(GLOB_RECURSE TMP_VTUS RELATIVE "${PROJECT_SOURCE_DIR}" "*tests/*.vtu")
file(GLOB_RECURSE TMP_VTPS RELATIVE "${PROJECT_SOURCE_DIR}" "tests/*.vtp")

list(APPEND TEST_DATA_FILES
	${TMP_GRIDS}
	${TMP_VTPS}
	${TMP_VTUS}
	)

list(APPEND TEST_SOURCE_FILES)

list(APPEND TEST_SOURCE_FILES)

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list(APPEND EXAMPLE_SOURCE_FILES
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	)

list (APPEND PUBLIC_HEADER_FILES
             ewoms/nonlinear/nullconvergencewriter.hh
             ewoms/nonlinear/newtonmethod.hh
             ewoms/io/vtktensorfunction.hh
             ewoms/io/dgfvanguard.hh
             ewoms/io/vtkscalarfunction.hh
             ewoms/io/vtkenergymodule.hh
             ewoms/io/restart.hh
             ewoms/io/cubegridvanguard.hh
             ewoms/io/baseoutputwriter.hh
             ewoms/io/vtkmultiwriter.hh
             ewoms/io/vtkmultiphasemodule.hh
             ewoms/io/vtkdiscretefracturemodule.hh
             ewoms/io/vtkdiffusionmodule.hh
             ewoms/io/vtkphasepresencemodule.hh
             ewoms/io/vtkblackoilmodule.hh
             ewoms/io/vtkblackoilsolventmodule.hh
             ewoms/io/vtkvectorfunction.hh
             ewoms/io/vtkprimaryvarsmodule.hh
             ewoms/io/simplexvanguard.hh
             ewoms/io/basevanguard.hh
             ewoms/io/vtktemperaturemodule.hh
             ewoms/io/vtkcompositionmodule.hh
             ewoms/io/structuredgridvanguard.hh
             ewoms/io/vtkblackoilenergymodule.hh
             ewoms/io/baseoutputmodule.hh
             ewoms/io/vtkblackoilpolymermodule.hh
             ewoms/common/start.hh
             ewoms/common/timerguard.hh
             ewoms/common/propertysystem.hh
             ewoms/common/pffgridvector.hh
             ewoms/common/prefetch.hh
             ewoms/common/parametersystem.hh
             ewoms/common/simulator.hh
             ewoms/common/quadraturegeometries.hh
             ewoms/common/alignedallocator.hh
             ewoms/common/timer.hh
             ewoms/common/signum.hh
             ewoms/common/genericguard.hh
             ewoms/common/basicproperties.hh
             ewoms/parallel/tasklets.hh
             ewoms/parallel/threadmanager.hh
             ewoms/parallel/gridcommhandles.hh
             ewoms/parallel/mpibuffer.hh
             ewoms/parallel/threadedentityiterator.hh
             ewoms/models/flash/flashmodel.hh
             ewoms/models/flash/flashintensivequantities.hh
             ewoms/models/flash/flashindices.hh
             ewoms/models/flash/flashlocalresidual.hh
             ewoms/models/flash/flashratevector.hh
             ewoms/models/flash/flashboundaryratevector.hh
             ewoms/models/flash/flashprimaryvariables.hh
             ewoms/models/flash/flashextensivequantities.hh
             ewoms/models/flash/flashproperties.hh
             ewoms/models/richards/richardsmodel.hh
             ewoms/models/richards/richardsextensivequantities.hh
             ewoms/models/richards/richardsratevector.hh
             ewoms/models/richards/richardsprimaryvariables.hh
             ewoms/models/richards/richardsnewtonmethod.hh
             ewoms/models/richards/richardsindices.hh
             ewoms/models/richards/richardsboundaryratevector.hh
             ewoms/models/richards/richardsproperties.hh
             ewoms/models/richards/richardsintensivequantities.hh
             ewoms/models/richards/richardslocalresidual.hh
             ewoms/models/discretefracture/discretefractureproblem.hh
             ewoms/models/discretefracture/discretefractureprimaryvariables.hh
             ewoms/models/discretefracture/discretefractureproperties.hh
             ewoms/models/discretefracture/fracturemapper.hh
             ewoms/models/discretefracture/discretefractureextensivequantities.hh
             ewoms/models/discretefracture/discretefracturemodel.hh
             ewoms/models/discretefracture/discretefractureintensivequantities.hh
             ewoms/models/discretefracture/discretefracturelocalresidual.hh
             ewoms/models/common/multiphasebaseproperties.hh
             ewoms/models/common/multiphasebasemodel.hh
             ewoms/models/common/quantitycallbacks.hh
             ewoms/models/common/multiphasebaseextensivequantities.hh
             ewoms/models/common/multiphasebaseproblem.hh
             ewoms/models/common/diffusionmodule.hh
             ewoms/models/common/flux.hh
             ewoms/models/common/forchheimerfluxmodule.hh
             ewoms/models/common/darcyfluxmodule.hh
             ewoms/models/common/energymodule.hh
             ewoms/models/pvs/pvsboundaryratevector.hh
             ewoms/models/pvs/pvsratevector.hh
             ewoms/models/pvs/pvsindices.hh
             ewoms/models/pvs/pvsproperties.hh
             ewoms/models/pvs/pvsnewtonmethod.hh
             ewoms/models/pvs/pvsprimaryvariables.hh
             ewoms/models/pvs/pvsextensivequantities.hh
             ewoms/models/pvs/pvsintensivequantities.hh
             ewoms/models/pvs/pvslocalresidual.hh
             ewoms/models/pvs/pvsmodel.hh
             ewoms/models/immiscible/immisciblelocalresidual.hh
             ewoms/models/immiscible/immiscibleproperties.hh
             ewoms/models/immiscible/immisciblemodel.hh
             ewoms/models/immiscible/immiscibleboundaryratevector.hh
             ewoms/models/immiscible/immiscibleratevector.hh
             ewoms/models/immiscible/immiscibleindices.hh
             ewoms/models/immiscible/immiscibleextensivequantities.hh
             ewoms/models/immiscible/immiscibleprimaryvariables.hh
             ewoms/models/immiscible/immiscibleintensivequantities.hh
             ewoms/models/ncp/ncpmodel.hh
             ewoms/models/ncp/ncpindices.hh
             ewoms/models/ncp/ncpextensivequantities.hh
             ewoms/models/ncp/ncpnewtonmethod.hh
             ewoms/models/ncp/ncpratevector.hh
             ewoms/models/ncp/ncpprimaryvariables.hh
             ewoms/models/ncp/ncpintensivequantities.hh
             ewoms/models/ncp/ncpproperties.hh
             ewoms/models/ncp/ncplocalresidual.hh
             ewoms/models/ncp/ncpboundaryratevector.hh
             ewoms/models/blackoil/blackoilmodel.hh
             ewoms/models/blackoil/blackoilextensivequantities.hh
             ewoms/models/blackoil/blackoilintensivequantities.hh
             ewoms/models/blackoil/blackoildarcyfluxmodule.hh
             ewoms/models/blackoil/blackoilratevector.hh
             ewoms/models/blackoil/blackoilfoammodules.hh
             ewoms/models/blackoil/blackoilindices.hh
             ewoms/models/blackoil/blackoillocalresidual.hh
             ewoms/models/blackoil/blackoilnewtonmethod.hh
             ewoms/models/blackoil/blackoilsolventmodules.hh
             ewoms/models/blackoil/blackoilproperties.hh
             ewoms/models/blackoil/blackoilprimaryvariables.hh
             ewoms/models/blackoil/blackoilproblem.hh
             ewoms/models/blackoil/blackoilenergymodules.hh
             ewoms/models/blackoil/blackoiltwophaseindices.hh
             ewoms/models/blackoil/blackoilpolymermodules.hh
             ewoms/models/blackoil/blackoilboundaryratevector.hh
             ewoms/disc/vcfv/vcfvbaseoutputmodule.hh
             ewoms/disc/vcfv/vcfvdiscretization.hh
             ewoms/disc/vcfv/p1fegradientcalculator.hh
             ewoms/disc/vcfv/vcfvgridcommhandlefactory.hh
             ewoms/disc/vcfv/vcfvproperties.hh
             ewoms/disc/vcfv/vcfvstencil.hh
             ewoms/disc/common/fvbasenewtonmethod.hh
             ewoms/disc/common/fvbasenewtonconvergencewriter.hh
             ewoms/disc/common/fvbaseintensivequantities.hh
             ewoms/disc/common/fvbaseconstraintscontext.hh
             ewoms/disc/common/baseauxiliarymodule.hh
             ewoms/disc/common/fvbaseelementcontext.hh
             ewoms/disc/common/fvbaselocalresidual.hh
             ewoms/disc/common/fvbasefdlocallinearizer.hh
             ewoms/disc/common/fvbaseboundarycontext.hh
             ewoms/disc/common/fvbaseadlocallinearizer.hh
             ewoms/disc/common/fvbaseconstraints.hh
             ewoms/disc/common/fvbaseproperties.hh
             ewoms/disc/common/fvbaseextensivequantities.hh
             ewoms/disc/common/fvbaselinearizer.hh
             ewoms/disc/common/restrictprolong.hh
             ewoms/disc/common/fvbasediscretization.hh
             ewoms/disc/common/fvbasegradientcalculator.hh
             ewoms/disc/common/fvbaseproblem.hh
             ewoms/disc/common/fvbaseprimaryvariables.hh
             ewoms/disc/ecfv/ecfvgridcommhandlefactory.hh
             ewoms/disc/ecfv/ecfvstencil.hh
             ewoms/disc/ecfv/ecfvbaseoutputmodule.hh
             ewoms/disc/ecfv/ecfvdiscretization.hh
             ewoms/disc/ecfv/ecfvproperties.hh
             opm/simulators/linalg/parallelistlbackend.hh
             opm/simulators/linalg/weightedresidreductioncriterion.hh
             opm/simulators/linalg/vertexborderlistfromgrid.hh
             opm/simulators/linalg/linearsolverreport.hh
             opm/simulators/linalg/istlsparsematrixadapter.hh
             opm/simulators/linalg/istlpreconditionerwrappers.hh
             opm/simulators/linalg/residreductioncriterion.hh
             opm/simulators/linalg/overlappingbcrsmatrix.hh
             opm/simulators/linalg/blacklist.hh
             opm/simulators/linalg/parallelbasebackend.hh
             opm/simulators/linalg/overlappingblockvector.hh
             opm/simulators/linalg/parallelbicgstabbackend.hh
             opm/simulators/linalg/nullborderlistmanager.hh
             opm/simulators/linalg/overlappingoperator.hh
             opm/simulators/linalg/elementborderlistfromgrid.hh
             opm/simulators/linalg/combinedcriterion.hh
             opm/simulators/linalg/bicgstabsolver.hh
             opm/simulators/linalg/globalindices.hh
             opm/simulators/linalg/superlubackend.hh
             opm/simulators/linalg/matrixblock.hh
             opm/simulators/linalg/istlsolverwrappers.hh
             opm/simulators/linalg/overlaptypes.hh
             opm/simulators/linalg/overlappingpreconditioner.hh
             opm/simulators/linalg/domesticoverlapfrombcrsmatrix.hh
             opm/simulators/linalg/fixpointcriterion.hh
             opm/simulators/linalg/parallelamgbackend.hh
             opm/simulators/linalg/foreignoverlapfrombcrsmatrix.hh
             opm/simulators/linalg/overlappingscalarproduct.hh
             opm/simulators/linalg/convergencecriterion.hh)
