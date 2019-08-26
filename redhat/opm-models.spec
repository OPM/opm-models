%define tag final

Name: opm-models
Summary: OPM - Fully implicit models for flow and transport in porous media
Version: 2018.10
Release: 0
License: GPL-3.0+
Group:   Development/Libraries/C and C++
URL: 	 http://opm-project.org/ewoms
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: dune-common-devel openmpi environment-modules valgrind
BuildRequires: make pkgconfig openmpi-devel dune-grid-devel
BuildRequires: opm-material-devel opm-grid-devel opm-common-devel
BuildRequires: opm-material-openmpi-devel opm-grid-openmpi-devel opm-common-openmpi-devel
BuildRequires: opm-material-mpich-devel opm-grid-mpich-devel opm-common-mpich-devel
BuildRequires: dune-istl-devel dune-localfunctions-devel doxygen zlib-devel
BuildRequires: devtoolset-6-toolchain ecl-devel
BuildRequires:  openmpi-devel trilinos-openmpi-devel ptscotch-openmpi-devel scotch-devel
BuildRequires:  mpich-devel trilinos-mpich-devel ptscotch-mpich-devel
Requires:      opm-models-openmpi-devel opm-models-mpich-devel
%{?el6:BuildRequires: cmake3 boost148-devel}
%{!?el6:BuildRequires: cmake boost-devel}

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

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

%build
scl enable devtoolset-6 bash
mkdir serial
cd serial
%{?el6:cmake3} %{!?el6:cmake} -DENABLE_MPI=0 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}
cd ..

mkdir openmpi
cd openmpi
%{?el6:module load openmpi-x86_64}
%{?!el6:module load mpi/openmpi-x86_64}
%{?el6:cmake3} %{!?el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/openmpi -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DZOLTAN_ROOT=/usr/lib64/openmpi -DCMAKE_CXX_FLAGS=-I/usr/include/openmpi-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/openmpi-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/openmpi -DPTSCOTCH_INCLUDE_DIR=/usr/include/openmpi-x86_64 -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}
cd ..

mkdir mpich
cd mpich
%{?el6:module rm openmpi-x86_64}
%{?el6:module load mpich-x86_64}
%{?!el6:module rm mpi/openmpi-x86_64}
%{?!el6:module load mpi/mpich-x86_64}
%{?el6:cmake3} %{!?el6:cmake} -DUSE_MPI=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix}/lib64/mpich -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148} -DZOLTAN_ROOT=/usr/lib64/mpich -DCMAKE_CXX_FLAGS=-I/usr/include/mpich-x86_64/trilinos -DZOLTAN_INCLUDE_DIRS=/usr/include/mpich-x86_64/trilinos -DPTSCOTCH_ROOT=/usr/lib64/mpich -DPTSCOTCH_INCLUDE_DIR=/usr/include/mpich-x86_64 -DCMAKE_INSTALL_SYSCONFDIR=/etc ..
make %{?_smp_mflags}

# No symbols in a template-only library
%global debug_package %{nil}

%install
cd serial
make install DESTDIR=${RPM_BUILD_ROOT}
cd ..

cd openmpi
%{!?el6:make doc}
make install DESTDIR=${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/openmpi/include/* ${RPM_BUILD_ROOT}/usr/include/openmpi-x86_64/
cd ..

cd mpich
make install DESTDIR=${RPM_BUILD_ROOT}
mkdir -p ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/
mv ${RPM_BUILD_ROOT}/usr/lib64/mpich/include/* ${RPM_BUILD_ROOT}/usr/include/mpich-x86_64/

%clean
rm -fr %buildroot

%files
%defattr(-,root,root)
%{!?el6:
%doc README
%doc openmpi/doc/doxygen/html}

%files devel
%defattr(-,root,root)
%_includedir/*
/usr/lib/dunecontrol/*
/usr/lib/pkgconfig/*
%{_datadir}/*

%files openmpi-devel
%defattr(-,root,root)
%{_includedir}/openmpi-x86_64/*
%{_libdir}/openmpi/lib/dunecontrol/*
%{_libdir}/openmpi/lib/pkgconfig/*
%{_libdir}/openmpi/share/*

%files mpich-devel
%defattr(-,root,root)
%{_includedir}/mpich-x86_64/*
%{_libdir}/mpich/lib/dunecontrol/*
%{_libdir}/mpich/lib/pkgconfig/*
%{_libdir}/mpich/share/*
