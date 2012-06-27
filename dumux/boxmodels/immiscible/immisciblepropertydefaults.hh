// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \ingroup ImmiscibleBoxModel
 * \ingroup BoxProperties
 * \ingroup Properties
 * \file
 *
 * \brief Defines default values for the properties required by the
 *        twophase box model.
 */
#ifndef DUMUX_IMMISCIBLE_PROPERTY_DEFAULTS_HH
#define DUMUX_IMMISCIBLE_PROPERTY_DEFAULTS_HH

#include "immisciblemodel.hh"
#include "immiscibleindices.hh"
#include "immisciblefluxvariables.hh"
#include "immiscibleprimaryvariables.hh"
#include "immisciblevolumevariables.hh"
#include "immiscibleratevector.hh"
#include "immiscibleboundaryratevector.hh"
#include "immiscibleproperties.hh"

#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>

#include <dumux/material/fluidmatrixinteractions/mp/nullmateriallaw.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>

namespace Dumux {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxImmiscible, NumEq, GET_PROP_VALUE(TypeTag, NumPhases)); //!< set the number of equations to the number of phases
SET_INT_PROP(BoxImmiscible, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases); //!< The number of phases is determined by the fluid system
SET_INT_PROP(BoxImmiscible, NumComponents, GET_PROP_VALUE(TypeTag, NumPhases));   //!< Number of chemical species in the system

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(BoxImmiscible,
              LocalResidual,
              ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxImmiscible, Model, ImmiscibleModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxImmiscible, BaseProblem, BoxMultiPhaseProblem<TypeTag>);

//! the FluidState property
SET_PROP(BoxImmiscible, FluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::ImmiscibleFluidState<Scalar,
                                        FluidSystem,
                                        /*enableEnthalpy=*/false> type;
};

//! the RateVector property
SET_TYPE_PROP(BoxImmiscible, RateVector, ImmiscibleRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxImmiscible, BoundaryRateVector, ImmiscibleBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxImmiscible, PrimaryVariables, ImmisciblePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxImmiscible, VolumeVariables, ImmiscibleVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxImmiscible, FluxVariables, ImmiscibleFluxVariables<TypeTag>);

//! The indices required by the isothermal immiscible multi-phase model
SET_TYPE_PROP(BoxImmiscible, Indices, ImmiscibleIndices</*PVOffset=*/0>);

/*!
 * \brief Set the material law to the null law by default.
 */
SET_TYPE_PROP(BoxImmiscible,
              MaterialLaw, 
              Dumux::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases), typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxImmiscible,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BoxImmiscible,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxImmiscible,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxImmiscible, EnableSmoothUpwinding, false);

/////////////////////
// set slightly different properties for the single-phase case
/////////////////////

//! The fluid system to use by default
SET_TYPE_PROP(BoxImmiscibleOnePhase, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(BoxImmiscibleOnePhase, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
SET_BOOL_PROP(BoxImmiscibleOnePhase, VtkWriteSaturations, false);
SET_BOOL_PROP(BoxImmiscibleOnePhase, VtkWriteMobilities, false);
SET_BOOL_PROP(BoxImmiscibleOnePhase, VtkWriteRelativePermeabilities, false);

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
SET_PROP(BoxImmiscibleTwoPhase, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxImmiscibleTwoPhase, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxImmiscibleTwoPhase, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

}

}

#endif
