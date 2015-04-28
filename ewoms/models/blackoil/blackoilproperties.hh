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
 * \ingroup BlackOilModel
 *
 * \brief Declares the properties required by the black oil model.
 */
#ifndef EWOMS_BLACK_OIL_PROPERTIES_HH
#define EWOMS_BLACK_OIL_PROPERTIES_HH

#include <ewoms/models/common/multiphasebaseproperties.hh>
#include <ewoms/io/vtkcompositionmodule.hh>
#include <ewoms/io/vtkblackoilmodule.hh>

namespace Ewoms {
namespace Properties {
//! Specifies if the simulation should write output files that are
//! compatible with those produced by the commercial Eclipse simulator
NEW_PROP_TAG(EnableEclipseOutput);
//! The fluid state used by the model
NEW_PROP_TAG(BlackOilFluidState);
//! The material law for heat conduction
NEW_PROP_TAG(HeatConductionLaw);
//! The parameters of the material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams);
//! Number of Newton-Raphson iterations for which the update should be chopped
NEW_PROP_TAG(BlackoilNumChoppedIterations);

}} // namespace Properties, Opm

#endif
