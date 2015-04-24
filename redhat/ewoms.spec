%define tag final

Name: ewoms
Summary: OPM - Fully implicit models for flow and transport in porous media
Version: 2015.04
Release: 0
License: GPL-3.0+
Group:   Development/Libraries/C and C++
URL: 	 http://opm-project.org/ewoms
Source0:        https://github.com/OPM/%{name}/archive/release/%{version}/%{tag}.tar.gz#/%{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: cmake28 dune-common-devel openmpi environment-modules valgrind
BuildRequires: make pkgconfig openmpi-devel dune-grid-devel ALUGrid-devel 
BuildRequires: alberta-devel opm-core-devel opm-material-devel
BuildRequires: dune-istl-devel dune-localfunctions-devel doxygen
Requires:      ewoms-devel
%{?el5:BuildRequires: gcc44 gcc44-c++}
%{!?el5:BuildRequires: gcc gcc-c++}
%{?el5:BuildRequires: boost141-devel}
%{!?el5:BuildRequires: boost-devel}

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

%prep
%setup -q -n %{name}-release-%{version}-%{tag}

%build
%{!?el5: module add openmpi-%{_arch}}
cmake28 -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DSTRIP_DEBUGGING_SYMBOLS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DCMAKE_INSTALL_DOCDIR=share/doc/%{name}-%{version} -DUSE_RUNPATH=OFF %{?el5:-DCMAKE_CXX_COMPILER=g++44 -DCMAKE_C_COMPILER=gcc44 -DBOOST_LIBRARYDIR=%{_libdir}/boost141 -DBOOST_INCLUDEDIR=/usr/include/boost141}
%__make %{?jobs:-j%{jobs}}

# No symbols in a template-only library
%global debug_package %{nil}

%install
make install DESTDIR=${RPM_BUILD_ROOT}
make doc

%clean
rm -fr %buildroot

%files
%defattr(-,root,root)
%doc README FAQ
%doc doc/doxygen/html

%files devel
%defattr(-,root,root)
%_includedir/*
%{_libdir}/dunecontrol/*
%{_libdir}/pkgconfig/*
%{_datadir}/*
