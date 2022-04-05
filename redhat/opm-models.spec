%define tag final
%define rtype release
%define toolset devtoolset-9
%define build_openmpi 1
%define build_openmpi3 1
%define build_mpich 1

Name: opm-models
Summary: OPM - Fully implicit models for flow and transport in porous media
Version:        2018.10
Release: 0
License: GPL-3.0+
Group:   Development/Libraries/C and C++
URL: 	 http://opm-project.org/ewoms
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: environment-modules openblas-devel
BuildRequires: make pkgconfig cmake3
BuildRequires: doxygen zlib-devel
BuildRequires: %{toolset}-toolchain
BuildRequires: boost-devel python3-devel tbb-devel
BuildRequires: dune-common-devel
BuildRequires: dune-uggrid-devel
BuildRequires: dune-grid-devel
BuildRequires: dune-istl-devel
BuildRequires: dune-localfunctions-devel
BuildRequires: opm-material-devel
BuildRequires: opm-grid-devel
BuildRequires: opm-common-devel

%if %{build_openmpi}
BuildRequires: zoltan-openmpi-devel
BuildRequires: openmpi-devel
BuildRequires: dune-common-openmpi-devel
BuildRequires: dune-istl-openmpi-devel
BuildRequires: dune-geometry-openmpi-devel
BuildRequires: dune-uggrid-openmpi-devel
BuildRequires: dune-grid-openmpi-devel
BuildRequires: dune-localfunctions-openmpi-devel
BuildRequires: opm-material-openmpi-devel
BuildRequires: opm-grid-openmpi-devel
BuildRequires: opm-common-openmpi-devel
%endif

%if %{build_openmpi3}
BuildRequires: zoltan-openmpi3-devel
BuildRequires: openmpi3-devel
BuildRequires: dune-common-openmpi3-devel
BuildRequires: dune-istl-openmpi3-devel
BuildRequires: dune-geometry-openmpi3-devel
BuildRequires: dune-uggrid-openmpi3-devel
BuildRequires: dune-grid-openmpi3-devel
BuildRequires: dune-localfunctions-openmpi3-devel
BuildRequires: opm-material-openmpi3-devel
BuildRequires: opm-grid-openmpi3-devel
BuildRequires: opm-common-openmpi3-devel
%endif

%if %{build_mpich}
BuildRequires: zoltan-mpich-devel
BuildRequires: mpich-devel
BuildRequires: dune-common-mpich-devel
BuildRequires: dune-istl-mpich-devel
BuildRequires: dune-geometry-mpich-devel
BuildRequires: dune-uggrid-mpich-devel
BuildRequires: dune-grid-mpich-devel
BuildRequires: dune-localfunctions-mpich-devel
BuildRequires: opm-material-mpich-devel
BuildRequires: opm-grid-mpich-devel
BuildRequires: opm-common-mpich-devel
%endif

%description
opm-models is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.

%package devel
License:      GPL
Requires:     %name = %version
Summary:     	opm-models development files
Group:        Development/Libraries/C and C++

%description devel
opm-models is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.

%if %{build_openmpi}
%package openmpi-devel
License:        GPL
Requires:      %name = %version
Summary:     	opm-models development files
Group:          Development/Libraries/C and C++

%description openmpi-devel
opm-models is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.
%endif

%if %{build_openmpi3}
%package openmpi3-devel
License:        GPL
Requires:      %name = %version
Summary:     	opm-models development files
Group:          Development/Libraries/C and C++

%description openmpi3-devel
opm-models is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.
%endif

%if %{build_mpich}
%package mpich-devel
License:        GPL
Requires:      %name = %version
Summary:     	opm-models development files
Group:          Development/Libraries/C and C++

%description mpich-devel
opm-models is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.
%endif

%prep
%setup -q -n %{name}-%{rtype}-%{version}-%{tag}

%build
mkdir serial
pushd serial
scl enable %{toolset} 'cmake3 -DENABLE_MPI=0 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
popd

%if %{build_openmpi}
mkdir openmpi
pushd openmpi
module load mpi/openmpi-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DZOLTAN_ROOT=/usr/lib64/openmpi -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
module unload mpi/openmpi-x86_64
popd
%endif

%if %{build_openmpi3}
mkdir openmpi3
pushd openmpi3
module load mpi/openmpi3-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi3 -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DZOLTAN_ROOT=/usr/lib64/openmpi3 -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
module unload mpi/openmpi3-x86_64
popd
%endif

%if %{build_mpich}
mkdir mpich
pushd mpich
module load mpi/mpich-x86_64
scl enable %{toolset} 'cmake3 -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DZOLTAN_ROOT=/usr/lib64/mpich -DCMAKE_INSTALL_SYSCONFDIR=/etc ..'
scl enable %{toolset} 'make %{?_smp_mflags}'
module unload mpi/mpich-x86_64
popd
%endif

# No symbols in a template-only library
%global debug_package %{nil}

%install
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C serial'

%if %{build_openmpi}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
%endif

%if %{build_openmpi3}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C openmpi3'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/openmpi3-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi3/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi3-x86_64/
%endif

%if %{build_mpich}
scl enable %{toolset} 'make install DESTDIR=${RPM_BUILD_ROOT} -C mpich'
mkdir -p ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
%endif

%clean
rm -fr %buildroot

%files
%defattr(-,root,root)
%doc README

%files devel
%defattr(-,root,root)
%_includedir/*
/usr/lib/dunecontrol/*
/usr/lib/pkgconfig/*
%{_datadir}/*

%if %{build_openmpi}
%files openmpi-devel
%defattr(-,root,root)
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_libdir}/openmpi/share/*
%endif

%if %{build_openmpi3}
%files openmpi3-devel
%defattr(-,root,root)
%{_includedir}/openmpi3-x86_64/*
%{_libdir}/openmpi3/lib/dunecontrol/*
%{_libdir}/openmpi3/lib/pkgconfig/*
%{_libdir}/openmpi3/share/*
%endif

%if %{build_mpich}
%files mpich-devel
%defattr(-,root,root)
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_libdir}/mpich/share/*
%endif
