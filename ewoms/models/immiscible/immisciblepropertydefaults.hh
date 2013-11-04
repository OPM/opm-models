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
 * \ingroup ImmiscibleVcfvModel
 *
 * \brief Defines default values for the properties required by the
 *        immiscible multi-phase model.
 */
#ifndef EWOMS_IMMISCIBLE_PROPERTY_DEFAULTS_HH
#define EWOMS_IMMISCIBLE_PROPERTY_DEFAULTS_HH

#include "immisciblemodel.hh"
#include "immiscibleindices.hh"
#include "immisciblefluxvariables.hh"
#include "immiscibleprimaryvariables.hh"
#include "immisciblevolumevariables.hh"
#include "immiscibleratevector.hh"
#include "immiscibleboundaryratevector.hh"
#include "immiscibleproperties.hh"

#include <ewoms/models/modules/velocity/vcfvvelocitymodules.hh>
#include <ewoms/disc/vcfv/vcfvmultiphaseproblem.hh>

#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/1pFluidSystem.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterialLaw.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Opm {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(VcfvImmiscible, NumEq, GET_PROP_TYPE(TypeTag, Indices)::numEq); //!< set the number of equations to the number of phases
SET_INT_PROP(VcfvImmiscible, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases); //!< The number of phases is determined by the fluid system
SET_INT_PROP(VcfvImmiscible, NumComponents, GET_PROP_VALUE(TypeTag, NumPhases));   //!< Number of chemical species in the system

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(VcfvImmiscible,
              LocalResidual,
              Ewoms::ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(VcfvImmiscible, Model, Ewoms::ImmiscibleModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(VcfvImmiscible, BaseProblem, Ewoms::VcfvMultiPhaseProblem<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvImmiscible, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(VcfvImmiscible, RateVector, Ewoms::ImmiscibleRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(VcfvImmiscible, BoundaryRateVector, Ewoms::ImmiscibleBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvImmiscible, PrimaryVariables, Ewoms::ImmisciblePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(VcfvImmiscible, VolumeVariables, Ewoms::ImmiscibleVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(VcfvImmiscible, FluxVariables, Ewoms::ImmiscibleFluxVariables<TypeTag>);

//! The indices required by the isothermal immiscible multi-phase model
SET_TYPE_PROP(VcfvImmiscible, Indices, Ewoms::ImmiscibleIndices<TypeTag, /*PVOffset=*/0>);

/*!
 * \brief Set the material law to the null law by default.
 */
SET_TYPE_PROP(VcfvImmiscible,
              MaterialLaw,
              Opm::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases), typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(VcfvImmiscible,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(VcfvImmiscible,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(VcfvImmiscible,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

// disable the smooth upwinding method by default
SET_BOOL_PROP(VcfvImmiscible, EnableSmoothUpwinding, false);

// disable the energy equation by default
SET_BOOL_PROP(VcfvImmiscible, EnableEnergy, false);

/////////////////////
// set slightly different properties for the single-phase case
/////////////////////

//! The fluid system to use by default
SET_TYPE_PROP(VcfvImmiscibleOnePhase, FluidSystem, Opm::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(VcfvImmiscibleOnePhase, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
SET_BOOL_PROP(VcfvImmiscibleOnePhase, VtkWriteSaturations, false);
SET_BOOL_PROP(VcfvImmiscibleOnePhase, VtkWriteMobilities, false);
SET_BOOL_PROP(VcfvImmiscibleOnePhase, VtkWriteRelativePermeabilities, false);

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
SET_PROP(VcfvImmiscibleTwoPhase, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(VcfvImmiscibleTwoPhase, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(VcfvImmiscibleTwoPhase, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Opm::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

}
}

#endif
