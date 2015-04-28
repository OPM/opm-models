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
 * \ingroup VcfvDiscretization
 *
 * \brief Declare the basic properties used by the common infrastructure of
 *        the vertex-centered finite volume discretization.
 */
#ifndef EWOMS_VCFV_PROPERTIES_HH
#define EWOMS_VCFV_PROPERTIES_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <ewoms/common/propertysystem.hh>

namespace Ewoms {
namespace Properties {
/*!
 * \ingroup VcfvDiscretization
 */
// \{

//! The type tag for models based on the VCFV-scheme
NEW_TYPE_TAG(VcfvDiscretization, INHERITS_FROM(FvBaseDiscretization));

//! Use the two-point gradient approximation scheme instead of first
//! order finite element ones
NEW_PROP_TAG(UseTwoPointGradients);

// \}
}} // namespace Properties, Opm


#endif
