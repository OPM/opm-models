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
 * \ingroup ImmiscibleModel
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

#include <ewoms/models/common/multiphasebasepropertydefaults.hh>

#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/1pFluidSystem.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>

namespace Opm {
namespace Properties {
//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, LocalResidual,
              Ewoms::ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(ImmiscibleModel, Model, Ewoms::ImmiscibleModel<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(ImmiscibleModel, RateVector, Ewoms::ImmiscibleRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(ImmiscibleModel, BoundaryRateVector, Ewoms::ImmiscibleBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(ImmiscibleModel, PrimaryVariables, Ewoms::ImmisciblePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(ImmiscibleModel, VolumeVariables, Ewoms::ImmiscibleVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(ImmiscibleModel, FluxVariables, Ewoms::ImmiscibleFluxVariables<TypeTag>);

//! The indices required by the isothermal immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, Indices, Ewoms::ImmiscibleIndices<TypeTag, /*PVOffset=*/0>);

//! disable the energy equation by default
SET_BOOL_PROP(ImmiscibleModel, EnableEnergy, false);

// set slightly different properties for the single-phase case

//! The fluid system to use by default
SET_TYPE_PROP(ImmiscibleOnePhaseModel, FluidSystem, Opm::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(ImmiscibleOnePhaseModel, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteSaturations, false);
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteMobilities, false);
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteRelativePermeabilities, false);

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
SET_PROP(ImmiscibleTwoPhaseModel, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(ImmiscibleTwoPhaseModel, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(ImmiscibleTwoPhaseModel, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Opm::FluidSystems::TwoPImmiscible<Scalar, WettingPhase,
                                              NonwettingPhase> type;
};

}} // namespace Properties, Opm

#endif
