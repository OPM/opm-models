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
 * \ingroup ImmiscibleModel
 *
 * \brief Defines the properties required for the immiscible
 *        multi-phase model.
 */
#ifndef EWOMS_IMMISCIBLE_PROPERTIES_HH
#define EWOMS_IMMISCIBLE_PROPERTIES_HH

#include <ewoms/models/common/multiphasebaseproperties.hh>
#include <ewoms/vtk/vtkenergymodule.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about the phases
NEW_PROP_TAG(EnableEnergy); //!< Specify whether energy should be considered as a conservation quantity or not

// these properties only make sense for the ImmiscibleTwoPhase type tag
NEW_PROP_TAG(WettingPhase); //!< The wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase); //!< The non-wetting phase for two-phase models

// these properties only make sense for the ImmiscibleOnePhase type tag
NEW_PROP_TAG(Fluid); //!< The fluid used by the model

} // namespace Properties
} // namespace Opm

#endif
