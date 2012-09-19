// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
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
 * \file
 *
 * \brief Contains the property declarations for the Richards box
 *        model.
 */
#ifndef DUMUX_RICHARDS_PROPERTIES_HH
#define DUMUX_RICHARDS_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>
#include <dumux/boxmodels/vtk/boxvtkmultiphasemodule.hh>
#include <dumux/boxmodels/vtk/boxvtktemperaturemodule.hh>

namespace Dumux
{
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal richards model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for problems discretized using the Richards model
NEW_TYPE_TAG(BoxRichards, INHERITS_FROM(BoxModel, VtkMultiPhase, VtkTemperature));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of components in the system
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used
NEW_PROP_TAG(MaterialLawParams); //!< The type of the parameter object for the material law
NEW_PROP_TAG(HeatConductionLaw); //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law for heat conduction
NEW_PROP_TAG(FluidSystem); //!< The fluid system to be used for the Richards model
NEW_PROP_TAG(LiquidPhaseIndex); //!< Index of the fluid which represents the wetting phase
NEW_PROP_TAG(WettingPhase); //!< Fluid which represents the wetting phase
NEW_PROP_TAG(NonwettingPhase); //!< Fluid which represents the non-wetting phase
NEW_PROP_TAG(VelocityModule); //!< Specifies the relation used for velocity
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specifies whether the smooth upwinding method should be used
// \}
}

} // end namepace

#endif
