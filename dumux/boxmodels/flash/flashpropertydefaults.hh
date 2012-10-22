// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \brief Defines default values for the properties required by the
 *        compositional multi-phase box model based on flash calculations.
 */
#ifndef DUMUX_FLASH_PROPERTY_DEFAULTS_HH
#define DUMUX_FLASH_PROPERTY_DEFAULTS_HH

#include "flashmodel.hh"
#include "flashprimaryvariables.hh"
#include "flashratevector.hh"
#include "flashboundaryratevector.hh"
#include "flashvolumevariables.hh"
#include "flashfluxvariables.hh"
#include "flashindices.hh"
#include "flashproperties.hh"

#include <dumux/boxmodels/modules/velocity/boxvelocitymodules.hh>
#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include <dumux/material/fluidmatrixinteractions/mp/nullmateriallaw.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 */
SET_INT_PROP(BoxFlash, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_INT_PROP(BoxFlash, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the number of PDEs to the number of compontents
 */
SET_INT_PROP(BoxFlash, 
             NumEq,
             GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material law to the dummy law.
 */
SET_TYPE_PROP(BoxFlash,
              MaterialLaw, 
              Dumux::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases), typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxFlash,
              MaterialLawParams, 
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BoxFlash,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxFlash,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the FlashLocalResidual function for the flash model
SET_TYPE_PROP(BoxFlash,
              LocalResidual,
              FlashLocalResidual<TypeTag>);

//! Use the NCP flash solver by default
SET_TYPE_PROP(BoxFlash,
              FlashSolver,
              Dumux::NcpFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
                              typename GET_PROP_TYPE(TypeTag, FluidSystem)>);

//! Let the flash solver choose its tolerance by default
SET_SCALAR_PROP(BoxFlash, FlashTolerance, 0.0);

//! the Model property
SET_TYPE_PROP(BoxFlash, Model, FlashModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxFlash, BaseProblem, BoxMultiPhaseProblem<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(BoxFlash, VelocityModule, Dumux::BoxDarcyVelocityModule<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxFlash, PrimaryVariables, FlashPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BoxFlash, RateVector, FlashRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxFlash, BoundaryRateVector, FlashBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxFlash, VolumeVariables, FlashVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxFlash, FluxVariables, FlashFluxVariables<TypeTag>);

//! The indices required by the flash-baseed isothermal compositional model
SET_TYPE_PROP(BoxFlash, Indices, FlashIndices<TypeTag, /*PVIdx=*/0>);

// use an isothermal model by default
SET_BOOL_PROP(BoxFlash, EnableEnergy, false);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxFlash, EnableSmoothUpwinding, false);

} // end namespace Properties
} // end namespace Dumux

#endif
