#define EWOMS_VERSION_MAJOR 2015
#define EWOMS_VERSION_MINOR 04
#define EWOMS_VERSION_REVISION -1 // -1 means that this is a version from the development branch...

#define EWOMS_VERSION_SUFFIX "git"
#define EWOMS_VERSION_CODENAME "Ringo"

#include <string>

namespace Ewoms {

inline std::string versionString()
{
    std::string tmp;
    tmp =
        std::to_string(EWOMS_VERSION_MAJOR) + "."
        + std::to_string(EWOMS_VERSION_MINOR);

    if (EWOMS_VERSION_REVISION >= 0) {
        tmp += ".";
        tmp +=  std::to_string(EWOMS_VERSION_REVISION);
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
