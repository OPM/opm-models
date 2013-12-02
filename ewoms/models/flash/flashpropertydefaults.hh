// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2012 by Andreas Lauser

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
 * \ingroup FlashModel
 *
 * \brief Defines default values for the properties required by the
 *        compositional multi-phase VCVF discretization based on flash
 *calculations.
 */
#ifndef EWOMS_FLASH_PROPERTY_DEFAULTS_HH
#define EWOMS_FLASH_PROPERTY_DEFAULTS_HH

#include "flashmodel.hh"
#include "flashprimaryvariables.hh"
#include "flashlocalresidual.hh"
#include "flashratevector.hh"
#include "flashboundaryratevector.hh"
#include "flashvolumevariables.hh"
#include "flashfluxvariables.hh"
#include "flashindices.hh"
#include "flashproperties.hh"

#include <ewoms/models/modules/velocity/vcfvvelocitymodules.hh>
#include <ewoms/disc/vcfv/vcfvmultiphaseproblem.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

namespace Opm {
namespace Properties {
/*!
 * \brief Set the property for the number of components.
 */
SET_INT_PROP(VcfvFlash, NumComponents,
             GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_INT_PROP(VcfvFlash, NumPhases,
             GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the number of PDEs to the number of compontents
 */
SET_INT_PROP(VcfvFlash, NumEq, GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material law to the dummy law.
 */
SET_PROP(VcfvFlash, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef NullMaterialTraits<Scalar, FluidSystem::numPhases> Traits;

public:
    typedef Opm::NullMaterial<Traits> type;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(VcfvFlash, MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(VcfvFlash, HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(VcfvFlash, HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the FlashLocalResidual function for the flash model
SET_TYPE_PROP(VcfvFlash, LocalResidual, Ewoms::FlashLocalResidual<TypeTag>);

//! Use the NCP flash solver by default
SET_TYPE_PROP(VcfvFlash, FlashSolver,
              Opm::NcpFlash<typename GET_PROP_TYPE(TypeTag, Scalar),
                            typename GET_PROP_TYPE(TypeTag, FluidSystem)>);

//! Let the flash solver choose its tolerance by default
SET_SCALAR_PROP(VcfvFlash, FlashTolerance, 0.0);

//! the Model property
SET_TYPE_PROP(VcfvFlash, Model, Ewoms::FlashModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(VcfvFlash, BaseProblem, Ewoms::VcfvMultiPhaseProblem<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvFlash, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvFlash, PrimaryVariables, Ewoms::FlashPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(VcfvFlash, RateVector, Ewoms::FlashRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(VcfvFlash, BoundaryRateVector,
              Ewoms::FlashBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(VcfvFlash, VolumeVariables, Ewoms::FlashVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(VcfvFlash, FluxVariables, Ewoms::FlashFluxVariables<TypeTag>);

//! The indices required by the flash-baseed isothermal compositional model
SET_TYPE_PROP(VcfvFlash, Indices, Ewoms::FlashIndices<TypeTag, /*PVIdx=*/0>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(VcfvFlash, EnableSmoothUpwinding, false);

// use an isothermal model by default
SET_BOOL_PROP(VcfvFlash, EnableEnergy, false);

// disable molecular diffusion by default
SET_BOOL_PROP(VcfvFlash, EnableDiffusion, false);

} // namespace Properties
} // namespace Opm

#endif
