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
 * \ingroup NcpModel
 *
 * \brief Declares the properties required for the NCP compositional
 *        multi-phase model.
 */
#ifndef EWOMS_NCP_PROPERTIES_HH
#define EWOMS_NCP_PROPERTIES_HH

#include <ewoms/models/common/multiphasebaseproperties.hh>
#include <ewoms/vtk/vtkcompositionmodule.hh>
#include <ewoms/vtk/vtkenergymodule.hh>
#include <ewoms/vtk/vtkdiffusionmodule.hh>

namespace Opm {
namespace Properties {
//! Enable the energy equation?
NEW_PROP_TAG(EnableEnergy);

//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

NEW_PROP_TAG(NcpPressureBaseWeight); //!< The unmodified weight for the pressure primary variable
NEW_PROP_TAG(NcpSaturationsBaseWeight); //!< The weight for the saturation primary variables
NEW_PROP_TAG(NcpFugacitiesBaseWeight); //!< The unmodified weight for the fugacity primary variables

//! The themodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
NEW_PROP_TAG(NcpCompositionFromFugacitiesSolver);

//! Number of Newton iterations per time step where the update gets chopped?
NEW_PROP_TAG(NcpNewtonNumChoppedIterations);
} // namespace Properties
} // namespace Opm

#endif
