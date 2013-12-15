// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2013 by Andreas Lauser

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
 * \copydoc Ewoms::EcfvVtkBaseOutputModule
 */
#ifndef EWOMS_ECFV_VTK_BASE_OUTPUT_MODULE_HH
#define EWOMS_ECFV_VTK_BASE_OUTPUT_MODULE_HH

#include "ecfvproperties.hh"

#include <ewoms/io/vtkmultiwriter.hh>

namespace Ewoms {
/*!
 * \ingroup EcfvVtkBaseOutputModule
 *
 * \brief Implements the discretization specific parts of writing VTK files.
 */
template<class TypeTag>
class EcfvVtkBaseOutputModule
{
public:
    /*!
     * \brief Add a buffer where the data is associated with the
     *        degrees of freedom to the current VTK output file.
     */
    template <class MultiWriter, class Buffer>
    static void attachDofData_(MultiWriter &writer, Buffer &buffer, const std::string &name, int numComponents)
    { writer.attachElementData(buffer, name, numComponents); }
};

} // namespace Ewoms

#endif
