Name: ewoms
Summary: OPM - Fully implicit models for flow and transport in porous media
Version: 2.3_git
Release: 0
License: GPL-3.0+
Group:   Development/Libraries/C and C++
URL: 	 http://opm-project.org/ewoms
Source:  %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: cmake dune-common-devel openmpi environment-modules valgrind
BuildRequires: make gcc gcc-c++ pkgconfig openmpi-devel dune-grid-devel ALUGrid-devel alberta-devel
BuildRequires: dune-istl-devel dune-localfunctions-devel doxygen texlive-latex
Requires:      ewoms-devel

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
%setup -q -n %{name}-%{version}

%build
module add openmpi-%{_arch}
cmake . -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%{_prefix}
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
