// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#define EWOMS_VERSION_MAJOR 2016
#define EWOMS_VERSION_MINOR 04
#define EWOMS_VERSION_REVISION -1 // -1 means that this is a version from the development branch...

#define EWOMS_VERSION_SUFFIX "pre"
#define EWOMS_VERSION_CODENAME "Yolanda"

#include <string>
#include <iostream>
#include <iomanip>

namespace Ewoms {

inline std::string versionString()
{
    std::ostringstream oss;
    oss << EWOMS_VERSION_MAJOR << "."
        << std::setfill('0') << std::setw(2)  << EWOMS_VERSION_MINOR;

    if (EWOMS_VERSION_REVISION >= 0)
        oss << "." << EWOMS_VERSION_REVISION;

    // append the version suffix to the version string
#ifdef EWOMS_VERSION_SUFFIX
    oss << "-" << EWOMS_VERSION_SUFFIX;
#endif

    // append the code name to the version string
#ifdef EWOMS_VERSION_CODENAME
    oss << " (\"" << EWOMS_VERSION_CODENAME "\")";
#endif

    return oss.str();
}

} // namespace Ewoms
