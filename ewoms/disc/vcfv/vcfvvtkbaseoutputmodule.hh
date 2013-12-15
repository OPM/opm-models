// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \copydoc Ewoms::VcfvVtkBaseOutputModule
 */
#ifndef EWOMS_VCFV_VTK_BASE_OUTPUT_MODULE_HH
#define EWOMS_VCFV_VTK_BASE_OUTPUT_MODULE_HH

#include "vcfvproperties.hh"

#include <ewoms/io/vtkmultiwriter.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup VcfvVtkBaseOutputModule
 *
 * \brief Implements the discretization specific parts of writing VTK files.
 */
template<class TypeTag>
class VcfvVtkBaseOutputModule
{
public:
    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    template <class MultiWriter, class Buffer>
    static void attachDofData_(MultiWriter &writer, Buffer &buffer, const std::string &name, int numComponents)
    { writer.attachVertexData(buffer, name.c_str(), numComponents); }
};

} // namespace Ewoms

#endif
