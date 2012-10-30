// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \ingroup NcpModel
 *
 * \brief Declares the properties required for the NCP compositional
 *        multi-phase box model.
 */
#ifndef DUMUX_NCP_PROPERTIES_HH
#define DUMUX_NCP_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>
#include <dumux/boxmodels/vtk/boxvtkmultiphasemodule.hh>
#include <dumux/boxmodels/vtk/boxvtkcompositionmodule.hh>
#include <dumux/boxmodels/vtk/boxvtktemperaturemodule.hh>
#include <dumux/boxmodels/vtk/boxvtkenergymodule.hh>
#include <dumux/boxmodels/vtk/boxvtkdiffusionmodule.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

/*!
 * \brief Define the type tag for the compositional NCP model.
 */
NEW_TYPE_TAG(BoxNcp, INHERITS_FROM(BoxModel, VtkMultiPhase, VtkComposition, VtkTemperature, VtkEnergy, VtkDiffusion));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(NcpEnergyIndices); //!< Enumerations for the 2pNc model

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the soil)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the soil)

NEW_PROP_TAG(HeatConductionLaw); //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law for heat conduction

//! The compositional twophase system of fluids which is considered
NEW_PROP_TAG(FluidSystem);

//! The themodynamic constraint solver which calculates the
//! composition of any phase given all component fugacities.
NEW_PROP_TAG(CompositionFromFugacitiesSolver);

//! Enable the energy equation?
NEW_PROP_TAG(EnableEnergy);

//! Enable diffusive fluxes?
NEW_PROP_TAG(EnableDiffusion);

//! Enable gravity?
NEW_PROP_TAG(EnableGravity);

//! Use the smooth upwinding method?
NEW_PROP_TAG(EnableSmoothUpwinding);

//! Specifies the relation used for velocity
NEW_PROP_TAG(VelocityModule);

//! Number of Newton iterations per time step where the update gets chopped?
NEW_PROP_TAG(NewtonChoppedIterations);
}
}

#endif
