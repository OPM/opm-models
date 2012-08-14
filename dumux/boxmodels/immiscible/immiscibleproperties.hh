// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup ImmiscibleBoxModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the twophase BOX model.
 */

#ifndef DUMUX_IMMISCIBLE_PROPERTIES_HH
#define DUMUX_IMMISCIBLE_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>
#include <dumux/boxmodels/vtk/boxvtkmultiphasemodule.hh>
#include <dumux/boxmodels/vtk/boxvtktemperaturemodule.hh>
#include <dumux/boxmodels/vtk/boxvtkenergymodule.hh>

namespace Dumux
{

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(BoxImmiscible, INHERITS_FROM(BoxModel, VtkMultiPhase, VtkTemperature, VtkEnergy));
//! The type tag for single-phase immiscible problems 
NEW_TYPE_TAG(BoxImmiscibleOnePhase, INHERITS_FROM(BoxImmiscible));
//! The type tag for two-phase immiscible problems
NEW_TYPE_TAG(BoxImmiscibleTwoPhase, INHERITS_FROM(BoxImmiscible));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of chemical species in the system
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(HeatConductionLaw); //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law for heat conduction
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about the phases
NEW_PROP_TAG(FluidState); //!<The phases state
NEW_PROP_TAG(EnableEnergy); //!< Specify whether energy should be considered as a conservation quantity or not

// these properties only make sense for the BoxImmiscibleTwoPhase type tag
NEW_PROP_TAG(WettingPhase); //!< The wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase); //!< The non-wetting phase for two-phase models

// these properties only make sense for the BoxImmiscibleOnePhase type tag
NEW_PROP_TAG(Fluid); //!< The fluid used by the model

} // namespace Properties
} // namespace Dumux

#endif
