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
 *        multi-phase VCVF discretization.
 */
#ifndef EWOMS_NCP_PROPERTIES_HH
#define EWOMS_NCP_PROPERTIES_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <ewoms/vtk/vcfvvtkmultiphasemodule.hh>
#include <ewoms/vtk/vcfvvtkcompositionmodule.hh>
#include <ewoms/vtk/vcfvvtktemperaturemodule.hh>
#include <ewoms/vtk/vcfvvtkenergymodule.hh>
#include <ewoms/vtk/vcfvvtkdiffusionmodule.hh>

namespace Opm {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

/*!
 * \brief Define the type tag for the compositional NCP model.
 */
NEW_TYPE_TAG(VcfvNcp, INHERITS_FROM(VcfvModel, VtkMultiPhase, VtkComposition, VtkTemperature, VtkEnergy, VtkDiffusion));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations used by the model

NEW_PROP_TAG(NcpPressureBaseWeight); //!< The unmodified weight for the pressure primary variable
NEW_PROP_TAG(NcpSaturationsBaseWeight); //!< The weight for the saturation primary variables
NEW_PROP_TAG(NcpFugacitiesBaseWeight); //!< The unmodified weight for the fugacity primary variables

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
