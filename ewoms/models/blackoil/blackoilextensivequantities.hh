/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Ewoms::BlackOilExtensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_EXTENSIVE_QUANTITIES_HH
#define EWOMS_BLACK_OIL_EXTENSIVE_QUANTITIES_HH

#include "blackoilproperties.hh"

#include <ewoms/models/common/multiphasebaseextensivequantities.hh>

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 * \ingroup ExtensiveQuantities
 *
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the black-oil model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class BlackOilExtensiveQuantities
    : public MultiPhaseBaseExtensiveQuantities<TypeTag>
{};

} // namespace Ewoms

#endif
