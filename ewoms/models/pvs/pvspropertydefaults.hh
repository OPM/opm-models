// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \ingroup PvsModel
 *
 * \brief Defines default values for the properties required for the
 *        primary variable switching model.
 */
#ifndef EWOMS_PVS_PROPERTY_DEFAULTS_HH
#define EWOMS_PVS_PROPERTY_DEFAULTS_HH

#include "pvslocalresidual.hh"
#include "pvsmodel.hh"
#include "pvsnewtonmethod.hh"
#include "pvsprimaryvariables.hh"
#include "pvsratevector.hh"
#include "pvsboundaryratevector.hh"
#include "pvsvolumevariables.hh"
#include "pvsfluxvariables.hh"
#include "pvsindices.hh"
#include "pvsproperties.hh"

#include <ewoms/models/common/multiphasebaseproblem.hh>
#include <ewoms/models/common/velocity.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Opm {
namespace Properties {
/*!
 * \brief Set the property for the number of components.
 */
SET_INT_PROP(PvsModel, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_INT_PROP(PvsModel, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the number of PDEs to the number of compontents
 */
SET_INT_PROP(PvsModel,
             NumEq,
             GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material law to the dummy law.
 */
SET_PROP(PvsModel, MaterialLaw)
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
SET_TYPE_PROP(PvsModel,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(PvsModel,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(PvsModel,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the PVS local jacobian operator for the PVS model
SET_TYPE_PROP(PvsModel,
              LocalResidual,
              Ewoms::PvsLocalResidual<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(PvsModel, VelocityModule, Ewoms::DarcyVelocityModule<TypeTag>);

//! Use the PVS specific newton method for the PVS model
SET_TYPE_PROP(PvsModel, NewtonMethod, Ewoms::PvsNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(PvsModel, Model, Ewoms::PvsModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(PvsModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(PvsModel, PrimaryVariables, Ewoms::PvsPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(PvsModel, RateVector, Ewoms::PvsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(PvsModel, BoundaryRateVector, Ewoms::PvsBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(PvsModel, VolumeVariables, Ewoms::PvsVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(PvsModel, FluxVariables, Ewoms::PvsFluxVariables<TypeTag>);

//! The indices required by the isothermal PVS model
SET_TYPE_PROP(PvsModel, Indices, Ewoms::PvsIndices<TypeTag, /*PVIdx=*/0>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(PvsModel, EnableSmoothUpwinding, false);

// set the model to a medium verbosity
SET_INT_PROP(PvsModel, PvsVerbosity, 1);

// disable the energy equation by default
SET_BOOL_PROP(PvsModel, EnableEnergy, false);

// disable molecular diffusion by default
SET_BOOL_PROP(PvsModel, EnableDiffusion, false);

//! The basis value for the weight of the pressure primary variable
SET_SCALAR_PROP(PvsModel, PvsPressureBaseWeight, 1.0);

//! The basis value for the weight of the saturation primary variables
SET_SCALAR_PROP(PvsModel, PvsSaturationsBaseWeight, 1.0);

//! The basis value for the weight of the mole fraction primary variables
SET_SCALAR_PROP(PvsModel, PvsMoleFractionsBaseWeight, 1.0);
}} // namespace Properties, Opm

#endif
