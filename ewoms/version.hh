// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#define EWOMS_VERSION_MAJOR 2015
#define EWOMS_VERSION_MINOR 10
#define EWOMS_VERSION_REVISION 0 // -1 means that this is a version from the development branch...

#define EWOMS_VERSION_SUFFIX "git"
#define EWOMS_VERSION_CODENAME "Pumpkin"

#include <string>

namespace Ewoms {

inline std::string versionString()
{
    std::string tmp;
    tmp =
        std::to_string(static_cast<long long>(EWOMS_VERSION_MAJOR)) + "."
        + std::to_string(static_cast<long long>(EWOMS_VERSION_MINOR));

    if (EWOMS_VERSION_REVISION >= 0) {
        tmp += ".";
        tmp +=  std::to_string(static_cast<long long>(EWOMS_VERSION_REVISION));
    }

    // append the version suffix to the version string
#ifdef EWOMS_VERSION_SUFFIX
    tmp += "-";
    tmp += EWOMS_VERSION_SUFFIX;
#endif

    // append the code name to the version string
#ifdef EWOMS_VERSION_CODENAME
    tmp += " (\"";
    tmp += EWOMS_VERSION_CODENAME;
    tmp += "\")";
#endif

    return tmp;
}

} // namespace Ewoms
