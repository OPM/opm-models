# -*- Autoconf -*-
DUNE_AC_INIT
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dumux/common/spline.hh])
AM_CONFIG_HEADER([config.h])

DUNE_CHECK_ALL
EWOMS_CHECKS_PRIVATE

AC_CONFIG_FILES([ewoms.pc
    Makefile
    doc/Makefile 
    doc/doxygen/Makefile 
    doc/doxygen/Doxyfile
    doc/handbook/Makefile
    dumux/Makefile 
    dumux/boxmodels/Makefile 
    dumux/boxmodels/blackoil/Makefile 
    dumux/boxmodels/common/Makefile 
    dumux/boxmodels/flash/Makefile 
    dumux/boxmodels/immiscible/Makefile 
    dumux/boxmodels/modules/Makefile 
    dumux/boxmodels/modules/energy/Makefile 
    dumux/boxmodels/ncp/Makefile 
    dumux/boxmodels/ncp/diffusion/Makefile
    dumux/boxmodels/pvs/Makefile 
    dumux/boxmodels/richards/Makefile 
    dumux/boxmodels/vtk/Makefile 
    dumux/common/Makefile 
    dumux/decoupled/Makefile
    dumux/decoupled/1p/Makefile
    dumux/decoupled/1p/diffusion/Makefile
    dumux/decoupled/1p/diffusion/fv/Makefile 
    dumux/decoupled/2p/Makefile 
    dumux/decoupled/2p/diffusion/Makefile 
    dumux/decoupled/2p/diffusion/fv/Makefile 
    dumux/decoupled/2p/diffusion/fvmpfa/Makefile
    dumux/decoupled/2p/diffusion/fvmpfa/lmethod/Makefile
    dumux/decoupled/2p/diffusion/fvmpfa/omethod/Makefile 
    dumux/decoupled/2p/diffusion/mimetic/Makefile 
    dumux/decoupled/2p/impes/Makefile 
    dumux/decoupled/2p/transport/Makefile 
    dumux/decoupled/2p/transport/fv/Makefile
    dumux/decoupled/2p2c/Makefile 
    dumux/decoupled/common/Makefile
    dumux/decoupled/common/fv/Makefile
    dumux/decoupled/common/fv/mpfa/Makefile
    dumux/decoupled/common/mimetic/Makefile
    dumux/decoupled/constraintsolvers/Makefile 
    dumux/decoupled/spatialparams/Makefile 
    dumux/freeflow/Makefile
    dumux/freeflow/stokes/Makefile
    dumux/io/Makefile
    dumux/linear/Makefile
    dumux/material/Makefile 
    dumux/material/binarycoefficients/Makefile 
    dumux/material/components/Makefile 
    dumux/material/components/iapws/Makefile 
    dumux/material/constraintsolvers/Makefile 
    dumux/material/eos/Makefile 
    dumux/material/fluidmatrixinteractions/Makefile 
    dumux/material/fluidmatrixinteractions/2p/Makefile 
    dumux/material/fluidmatrixinteractions/3p/Makefile 
    dumux/material/fluidmatrixinteractions/mp/Makefile 
    dumux/material/fluidstates/Makefile 
    dumux/material/fluidsystems/Makefile 
    dumux/material/heatconduction/Makefile 
    dumux/nonlinear/Makefile
    dumux/parallel/Makefile 
    m4/Makefile
    test/Makefile
    test/boxmodels/Makefile 
    test/boxmodels/problems/Makefile 
    test/boxmodels/grids/Makefile 
    test/boxmodels/quadrature/Makefile
    test/common/Makefile
    test/common/propertysystem/Makefile
    test/common/spline/Makefile
    test/decoupled/Makefile
    test/decoupled/1p/Makefile
    test/decoupled/2p/Makefile
    test/decoupled/2padaptive/Makefile
    test/decoupled/2p2c/Makefile
    test/freeflow/Makefile
    test/freeflow/navierstokes/Makefile
    test/freeflow/stokes/Makefile
    test/freeflow/stokes2c/Makefile
    test/freeflow/stokesni/Makefile
    test/material/Makefile
    test/material/fluidsystems/Makefile
    test/material/immiscibleflash/Makefile
    test/material/ncpflash/Makefile
    test/material/pengrobinson/Makefile
    test/material/tabulation/Makefile
    tutorial/Makefile
])

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, "$ALL_PKG_CPPFLAGS" )
AC_SUBST(AM_LDFLAGS, "$ALL_PKG_LDFLAGS" )

LIBS="$DUNE_LIBS $ALL_PKG_LIBS"

AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
