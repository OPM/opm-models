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
 * \ingroup EcfvDiscretization
 *
 * \brief Declare the basic properties used by the common infrastructure of
 *        the vertex-centered finite volume discretization.
 */
#ifndef EWOMS_ECFV_PROPERTIES_HH
#define EWOMS_ECFV_PROPERTIES_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <opm/core/utility/PropertySystem.hpp>

namespace Opm {
namespace Properties {
/*!
 * \ingroup EcfvDiscretization
 */
// \{

//! The type tag for models based on the ECFV-scheme
NEW_TYPE_TAG(EcfvDiscretization, INHERITS_FROM(FvBaseDiscretization));

// \}
}} // namespace Properties, Opm


#endif
