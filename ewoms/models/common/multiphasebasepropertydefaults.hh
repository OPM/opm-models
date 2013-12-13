// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2013 by Andreas Lauser

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
 * \ingroup MultiPhaseBaseModel
 *
 * \brief Defines default values for the common properties of the
 *        porous-media multi-phase models.
 */
#ifndef EWOMS_MULTI_PHASE_BASE_PROPERTY_DEFAULTS_HH
#define EWOMS_MULTI_PHASE_BASE_PROPERTY_DEFAULTS_HH

#include "multiphasebaseproblem.hh"
#include "multiphasebasefluxvariables.hh"

#include <ewoms/models/common/velocity.hh>
#include <ewoms/disc/common/fvbasepropertydefaults.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Opm {
namespace Properties {
SET_INT_PROP(MultiPhaseBaseModel, NumEq, GET_PROP_TYPE(TypeTag, Indices)::numEq); //!< set the number of equations to the number of phases
SET_INT_PROP(MultiPhaseBaseModel, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases); //!< The number of phases is determined by the fluid system
SET_INT_PROP(MultiPhaseBaseModel, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents); //!< Number of chemical species in the system

//! The type of the base base class for actual problems
SET_TYPE_PROP(MultiPhaseBaseModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(MultiPhaseBaseModel, VelocityModule, Ewoms::DarcyVelocityModule<TypeTag>);

/*!
 * \brief Set the material law to the null law by default.
 */
SET_PROP(MultiPhaseBaseModel, MaterialLaw)
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
SET_TYPE_PROP(MultiPhaseBaseModel,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(MultiPhaseBaseModel,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(MultiPhaseBaseModel,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! disable the smooth upwinding method by default
SET_BOOL_PROP(MultiPhaseBaseModel, EnableSmoothUpwinding, false);

//! disable gravity by default
SET_BOOL_PROP(MultiPhaseBaseModel, EnableGravity, false);

}} // namespace Properties, Opm

#endif
