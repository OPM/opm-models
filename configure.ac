# -*- Autoconf -*-
DUNE_AC_INIT
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([ewoms/common/spline.hh])
AC_CONFIG_HEADERS([config.h])

DUNE_CHECK_ALL
EWOMS_CHECKS_PRIVATE

AC_CONFIG_FILES([Makefile
    doc/Makefile 
    doc/doxygen/Makefile 
    doc/doxygen/Doxyfile
    doc/handbook/Makefile
    ewoms/Makefile 
    ewoms/models/Makefile 
    ewoms/models/blackoil/Makefile 
    ewoms/models/discretefracture/Makefile 
    ewoms/models/flash/Makefile 
    ewoms/models/immiscible/Makefile 
    ewoms/models/modules/Makefile 
    ewoms/models/modules/energy/Makefile 
    ewoms/models/modules/diffusion/Makefile 
    ewoms/models/modules/velocity/Makefile 
    ewoms/models/ncp/Makefile 
    ewoms/models/pvs/Makefile 
    ewoms/models/richards/Makefile 
    ewoms/models/stokes/Makefile
    ewoms/common/Makefile 
    ewoms/decoupled/Makefile
    ewoms/decoupled/1p/Makefile
    ewoms/decoupled/1p/diffusion/Makefile
    ewoms/decoupled/1p/diffusion/fv/Makefile 
    ewoms/decoupled/2p/Makefile 
    ewoms/decoupled/2p/diffusion/Makefile 
    ewoms/decoupled/2p/diffusion/fv/Makefile 
    ewoms/decoupled/2p/diffusion/fvmpfa/Makefile
    ewoms/decoupled/2p/diffusion/fvmpfa/lmethod/Makefile
    ewoms/decoupled/2p/diffusion/fvmpfa/omethod/Makefile 
    ewoms/decoupled/2p/diffusion/mimetic/Makefile 
    ewoms/decoupled/2p/impes/Makefile 
    ewoms/decoupled/2p/transport/Makefile 
    ewoms/decoupled/2p/transport/fv/Makefile
    ewoms/decoupled/2p2c/Makefile 
    ewoms/decoupled/common/Makefile
    ewoms/decoupled/common/fv/Makefile
    ewoms/decoupled/common/fv/mpfa/Makefile
    ewoms/decoupled/common/mimetic/Makefile
    ewoms/decoupled/constraintsolvers/Makefile 
    ewoms/decoupled/spatialparams/Makefile 
    ewoms/disc/Makefile 
    ewoms/disc/vcfv/Makefile 
    ewoms/io/Makefile
    ewoms/istl/Makefile
    ewoms/linear/Makefile
    ewoms/material/Makefile 
    ewoms/material/binarycoefficients/Makefile 
    ewoms/material/components/Makefile 
    ewoms/material/components/iapws/Makefile 
    ewoms/material/constraintsolvers/Makefile 
    ewoms/material/eos/Makefile 
    ewoms/material/fluidmatrixinteractions/Makefile 
    ewoms/material/fluidmatrixinteractions/2p/Makefile 
    ewoms/material/fluidmatrixinteractions/3p/Makefile 
    ewoms/material/fluidmatrixinteractions/mp/Makefile 
    ewoms/material/fluidstates/Makefile 
    ewoms/material/fluidsystems/Makefile 
    ewoms/material/heatconduction/Makefile 
    ewoms/nonlinear/Makefile
    ewoms/parallel/Makefile 
    ewoms/vtk/Makefile 
    m4/Makefile
    test/Makefile
    test/implicit/Makefile 
    test/implicit/problems/Makefile 
    test/implicit/grids/Makefile 
    test/common/Makefile
    test/common/propertysystem/Makefile
    test/common/spline/Makefile
    test/decoupled/Makefile
    test/decoupled/1p/Makefile
    test/decoupled/2p/Makefile
    test/decoupled/2p2c/Makefile
    test/material/Makefile
    test/material/fluidsystems/Makefile
    test/material/immiscibleflash/Makefile
    test/material/ncpflash/Makefile
    test/material/pengrobinson/Makefile
    test/material/tabulation/Makefile
    tutorial/Makefile
    ewoms.pc
    EwomsConfig.cmake
])

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, "$ALL_PKG_CPPFLAGS" )
AC_SUBST(AM_LDFLAGS, "$ALL_PKG_LDFLAGS" )

LIBS="$DUNE_LIBS $ALL_PKG_LIBS"

AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
