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
 * \brief Defines default values for the properties used by the
 *        black oil VCVF discretization.
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

#include <ewoms/models/modules/velocity/vcfvvelocitymodules.hh>
#include <ewoms/disc/vcfv/vcfvmultiphaseproblem.hh>
#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>
#include <ewoms/material/heatconduction/dummyheatconductionlaw.hh>

namespace Ewoms {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(VcfvBlackOil, NumEq, 3); //!< set the number of equations to 3
SET_INT_PROP(VcfvBlackOil, NumPhases, 3); //!< The number of phases considered by the model
SET_INT_PROP(VcfvBlackOil, NumComponents, 3);   //!< Number of chemical species in the system

//! Set the local residual function
SET_TYPE_PROP(VcfvBlackOil,
              LocalResidual,
              BlackOilLocalResidual<TypeTag>);

//! The Model property
SET_TYPE_PROP(VcfvBlackOil, Model, BlackOilModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(VcfvBlackOil, BaseProblem, VcfvMultiPhaseProblem<TypeTag>);

//! The BlackOilFluidState property
SET_PROP(VcfvBlackOil, BlackOilFluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Ewoms::CompositionalFluidState<Scalar,
                                           FluidSystem,
                                           /*enableEnthalpy=*/false> type;
};

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvBlackOil, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(VcfvBlackOil, RateVector, BlackOilRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(VcfvBlackOil, BoundaryRateVector, BlackOilBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvBlackOil, PrimaryVariables, BlackOilPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(VcfvBlackOil, VolumeVariables, BlackOilVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(VcfvBlackOil, FluxVariables, BlackOilFluxVariables<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(VcfvBlackOil, Indices, BlackOilIndices</*PVOffset=*/0>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(VcfvBlackOil,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(VcfvBlackOil,
              HeatConductionLaw,
              Ewoms::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(VcfvBlackOil,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Set the fluid system to the black-oil fluid system by default
SET_TYPE_PROP(VcfvBlackOil,
              FluidSystem,
              Ewoms::FluidSystems::BlackOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(VcfvBlackOil, EnableSmoothUpwinding, false);

}

}

#endif
