// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \ingroup BlackOilModel
 *
 * \brief Defines default values for the properties used by the
 *        black oil model.
 */
#ifndef EWOMS_BLACK_OIL_PROPERTY_DEFAULTS_HH
#define EWOMS_BLACK_OIL_PROPERTY_DEFAULTS_HH

#include "blackoilmodel.hh"
#include "blackoilindices.hh"
#include "blackoilfluxvariables.hh"
#include "blackoilprimaryvariables.hh"
#include "blackoilvolumevariables.hh"
#include "blackoilratevector.hh"
#include "blackoilboundaryratevector.hh"
#include "blackoilproperties.hh"

#include <ewoms/models/common/multiphasebaseproblem.hh>
#include <ewoms/models/modules/velocity.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Opm {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BlackOilModel, NumEq, 3); //!< set the number of equations to 3
SET_INT_PROP(BlackOilModel, NumPhases, 3); //!< The number of phases considered by the model
SET_INT_PROP(BlackOilModel, NumComponents, 3); //!< Number of chemical species in the system

//! Set the local residual function
SET_TYPE_PROP(BlackOilModel,
              LocalResidual,
              Ewoms::BlackOilLocalResidual<TypeTag>);

//! The Model property
SET_TYPE_PROP(BlackOilModel, Model, Ewoms::BlackOilModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BlackOilModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! The BlackOilFluidState property
SET_PROP(BlackOilModel, BlackOilFluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Opm::CompositionalFluidState<Scalar,
                                           FluidSystem,
                                           /*enableEnthalpy=*/false> type;
};

//! Use the Darcy relation by default
SET_TYPE_PROP(BlackOilModel, VelocityModule, Ewoms::DarcyVelocityModule<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BlackOilModel, RateVector, Ewoms::BlackOilRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BlackOilModel, BoundaryRateVector, Ewoms::BlackOilBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BlackOilModel, PrimaryVariables, Ewoms::BlackOilPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BlackOilModel, VolumeVariables, Ewoms::BlackOilVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BlackOilModel, FluxVariables, Ewoms::BlackOilFluxVariables<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(BlackOilModel, Indices, Ewoms::BlackOilIndices</*PVOffset=*/0>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BlackOilModel,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BlackOilModel,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BlackOilModel,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Set the fluid system to the black-oil fluid system by default
SET_TYPE_PROP(BlackOilModel,
              FluidSystem,
              Opm::FluidSystems::BlackOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BlackOilModel, EnableSmoothUpwinding, false);
}} // namespace Properties, Opm

#endif
