// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#define EWOMS_VERSION_MAJOR 2016
#define EWOMS_VERSION_MINOR 10
#define EWOMS_VERSION_REVISION 0 // -1 means that this is a version from the development branch...

#define EWOMS_VERSION_SUFFIX "rc1"
#define EWOMS_VERSION_CODENAME "Butch"

#include <string>
#include <sstream>
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
