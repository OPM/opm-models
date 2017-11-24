%define tag rc1

Name: ewoms
Summary: OPM - Fully implicit models for flow and transport in porous media
Version: 2017.10
Release: 0
License: GPL-3.0+
Group:   Development/Libraries/C and C++
URL: 	 http://opm-project.org/ewoms
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: dune-common-devel openmpi environment-modules valgrind
BuildRequires: make pkgconfig openmpi-devel dune-grid-devel opm-common-devel
BuildRequires: opm-core-devel opm-material-devel opm-grid-devel opm-output-devel
BuildRequires: dune-istl-devel dune-localfunctions-devel doxygen zlib-devel
BuildRequires: devtoolset-6-toolchain
Requires:      ewoms-devel
%{?el6:BuildRequires: cmake3 boost148-devel}
%{!?el6:BuildRequires: cmake boost-devel}

%description
eWoms is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.

%package devel
License:        GPL
Requires:      %name = %version
Summary:     	ewoms development files
Group:          Development/Libraries/C and C++

%description devel
eWoms is an simulation framework which is primary focused on fully implicit
models for flow and transport in porous media. Its main objectives
are to provide a easily usable, well maintainable, high performance
framework which is capable of capturing all macro-scale scenarios
relevant for academic research and industrial applications involving
flow and transport processes in porous media.

%package -n ebos
Summary:        ebos is an ECL simulator.
Group:          Scientific

%description -n ebos
ebos is an ECL simulator.

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

%build
module add openmpi-%{_arch}
scl enable devtoolset-6 bash
%{?el6:cmake3} %{!?el6:cmake} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF -DWITH_NATIVE=OFF -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++ -DCMAKE_C_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/gcc %{?el6:-DBOOST_LIBRARYDIR=%{_libdir}/boost148 -DBOOST_INCLUDEDIR=%{_includedir}/boost148}
%__make %{?jobs:-j%{jobs}}

# No symbols in a template-only library
%global debug_package %{nil}

%install
make install DESTDIR=${RPM_BUILD_ROOT}
%{!?el6:make doc}

%clean
rm -fr %buildroot

%files
%defattr(-,root,root)
%{!?el6:
%doc README
%doc doc/doxygen/html}

%files devel
%defattr(-,root,root)
%_includedir/*
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_datadir}/*

%files -n ebos
%{_bindir}/ebos
