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

#include <ewoms/disc/vcfv/vcfvmultiphaseproblem.hh>
#include <ewoms/models/modules/velocity/vcfvvelocitymodules.hh>

#include <opm/material/fluidmatrixinteractions/mp/NullMaterialLaw.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Ewoms {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 */
SET_INT_PROP(VcfvPvs, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_INT_PROP(VcfvPvs, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the number of PDEs to the number of compontents
 */
SET_INT_PROP(VcfvPvs,
             NumEq,
             GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material law to the dummy law.
 */
SET_TYPE_PROP(VcfvPvs,
              MaterialLaw,
              Opm::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases), typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(VcfvPvs,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(VcfvPvs,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(VcfvPvs,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the PVS local jacobian operator for the PVS model
SET_TYPE_PROP(VcfvPvs,
              LocalResidual,
              PvsLocalResidual<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvPvs, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

//! Use the PVS specific newton method for the PVS model
SET_TYPE_PROP(VcfvPvs, NewtonMethod, PvsNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(VcfvPvs, Model, PvsModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(VcfvPvs, BaseProblem, VcfvMultiPhaseProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvPvs, PrimaryVariables, PvsPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(VcfvPvs, RateVector, PvsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(VcfvPvs, BoundaryRateVector, PvsBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(VcfvPvs, VolumeVariables, PvsVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(VcfvPvs, FluxVariables, PvsFluxVariables<TypeTag>);

//! The indices required by the isothermal PVS model
SET_TYPE_PROP(VcfvPvs, Indices, PvsIndices<TypeTag, /*PVIdx=*/0>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(VcfvPvs, EnableSmoothUpwinding, false);

// set the model to a medium verbosity
SET_INT_PROP(VcfvPvs, PvsVerbosity, 1);

// disable the energy equation by default
SET_BOOL_PROP(VcfvPvs, EnableEnergy, false);

// disable molecular diffusion by default
SET_BOOL_PROP(VcfvPvs, EnableDiffusion, false);

//! The basis value for the weight of the pressure primary variable
SET_SCALAR_PROP(VcfvPvs, PvsPressureBaseWeight, 1.0);

//! The basis value for the weight of the saturation primary variables
SET_SCALAR_PROP(VcfvPvs, PvsSaturationsBaseWeight, 1.0);

//! The basis value for the weight of the mole fraction primary variables
SET_SCALAR_PROP(VcfvPvs, PvsMoleFractionsBaseWeight, 1.0);
}

}

#endif
