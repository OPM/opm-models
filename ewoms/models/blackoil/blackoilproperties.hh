// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \ingroup BlackOilVcfvModel
 *
 * \brief Declares the properties required by the black oil model.
 */
#ifndef EWOMS_BLACK_OIL_PROPERTIES_HH
#define EWOMS_BLACK_OIL_PROPERTIES_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <ewoms/vtk/vcfvvtkmultiphasemodule.hh>
#include <ewoms/vtk/vcfvvtkcompositionmodule.hh>
#include <ewoms/vtk/vcfvvtktemperaturemodule.hh>
#include <ewoms/vtk/vcfvvtkblackoilmodule.hh>

namespace Opm {
namespace Properties {
//! The type tag for the black-oil problems
NEW_TYPE_TAG(VcfvBlackOil, INHERITS_FROM(VcfvModel, VtkBlackOil, VtkMultiPhase, VtkComposition, VtkTemperature));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of chemical species in the system
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specifies whether the smooth upwinding method should be used
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(BlackOilFluidState); //!< The fluid state used by the model
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(HeatConductionLaw); //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law for heat conduction
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about the phases
NEW_PROP_TAG(FluidState); //!<The phases state
NEW_PROP_TAG(VelocityModule); //!< Specifies the relation used for velocity
NEW_PROP_TAG(FluidState); //!<The type of the fluid state
}
}

#endif
