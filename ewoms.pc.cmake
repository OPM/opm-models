prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_INSTALL_PREFIX@/bin
libdir=@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_LIBDIR@
includedir=@CMAKE_INSTALL_PREFIX@/include
CXX=@CMAKE_CXX_COMPILER@ -std=c++0x
CC=@CMAKE_C_COMPILER@
DEPENDENCIES=@ProjectDepends@

Name: @ProjectName@
Version: @ProjectVersion@
Description: @ProjectDescription@
URL: @ProjectUrl@
Requires: ${DEPENDENCIES}
Libs:
Cflags:  -I${includedir}
