// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup OnePBoxModel
 * \file
 *
 * \brief Defines the properties required for the one-phase BOX model.
 */
#ifndef DUMUX_1P_PROPERTY_DEFAULTS_HH
#define DUMUX_1P_PROPERTY_DEFAULTS_HH

#include <dumux/boxmodels/common/boxproperties.hh>

#include "1pmodel.hh"
#include "1plocalresidual.hh"
#include "1pindices.hh"
#include "1pprimaryvariables.hh"
#include "1pratevector.hh"
#include "1pvolumevariables.hh"
#include "1pfluxvariables.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>


namespace Dumux
{
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(BoxOneP, NumEq, 1); //!< set the number of equations to 1
SET_INT_PROP(BoxOneP, NumPhases, 1); //!< The number of phases in the 1p model is 1
SET_INT_PROP(BoxOneP, NumComponents, 1); //!< The number of (pseudo) components in the 1p model is 1

//! The local residual function
SET_TYPE_PROP(BoxOneP,
              LocalResidual,
              OnePLocalResidual<TypeTag>);

//! set the heat conduction law to the dummy law
SET_TYPE_PROP(BoxOneP,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxOneP,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! the Model property
SET_TYPE_PROP(BoxOneP, Model, OnePBoxModel<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BoxOneP, RateVector, OnePRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxOneP, PrimaryVariables, OnePPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxOneP, VolumeVariables, OnePVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxOneP, FluxVariables, OnePFluxVariables<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(BoxOneP, OnePIndices, OnePIndices);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(BoxOneP, UpwindWeight, 0.5);

//! The fluid system to use by default
SET_TYPE_PROP(BoxOneP, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(BoxOneP, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

// use the conjugated gradient solver preconditioned by ILU-0.
// WARNING: CG is only expected to work for symmetric Jacobian
// matrices, which means that the either fluid used is incompressible
// and/or you are _not_ use upwinding!
SET_TYPE_PROP(BoxOneP, LinearSolver, Dumux::BoxCGILU0Solver<TypeTag> );

// disable output of a few quantities which make sense in a
// multip-hase but not in a single-phase context
SET_BOOL_PROP(BoxOneP, VtkWriteSaturations, false);
SET_BOOL_PROP(BoxOneP, VtkWriteMobilities, false);
SET_BOOL_PROP(BoxOneP, VtkWriteRelativePermeabilities, false);

// enable filter velocity output by default
SET_BOOL_PROP(BoxOneP, VtkWriteFilterVelocities, true);
// \}
} // end namepace Properties

} // end namepace Dumux

#endif
