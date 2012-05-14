// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_MPNC_PROPERTY_DEFAULTS_HH
#define DUMUX_MPNC_PROPERTY_DEFAULTS_HH

#include "mpncmodel.hh"
#include "mpncproblem.hh"
#include "mpnclocalresidual.hh"
#include "mpncfluxvariables.hh"
#include "mpncprimaryvariables.hh"
#include "mpncboundaryratevector.hh"
#include "mpncratevector.hh"
#include "mpncvolumevariables.hh"
#include "mpncnewtoncontroller.hh"
#include "mpncindices.hh"
#include "mpncproperties.hh"

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>


/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxMpNcModel
 * \file
 * \brief  Default properties for the Mp-Nc box model.
 */
namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// default property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system.
 */
SET_INT_PROP(BoxMPNC, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_INT_PROP(BoxMPNC, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_INT_PROP(BoxMPNC, NumEq,GET_PROP_TYPE(TypeTag, Indices)::NumPrimaryVars);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxMPNC, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxMPNC,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxMPNC,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

/*!
 * \brief Set the themodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(BoxMPNC, CompositionFromFugacitiesSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> type;
};

//! Use the MpNc local jacobian operator for the MpNc model
SET_TYPE_PROP(BoxMPNC,
              LocalResidual,
              MPNCLocalResidual<TypeTag>);

//! Use the MpNc specific newton controller for the MpNc model
SET_TYPE_PROP(BoxMPNC, NewtonController, MPNCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxMPNC, Model, MPNCModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxMPNC, BaseProblem, MPNCProblem<TypeTag>);

//! use an isothermal model by default
SET_BOOL_PROP(BoxMPNC, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(BoxMPNC, EnableDiffusion, false);

//! do not use smooth upwinding by default
SET_BOOL_PROP(BoxMPNC, EnableSmoothUpwinding, false);

//! the RateVector property
SET_TYPE_PROP(BoxMPNC, RateVector, MPNCRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxMPNC, BoundaryRateVector, MPNCBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxMPNC, PrimaryVariables, MPNCPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxMPNC, VolumeVariables, MPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxMPNC, FluxVariables, MPNCFluxVariables<TypeTag>);

//! truncate the newton update for the first 2 iterations of a time step
SET_INT_PROP(BoxMPNC, NewtonChoppedIterations, 2);

//! The indices required by the compositional twophase model
SET_TYPE_PROP(BoxMPNC, Indices, MPNCIndices<TypeTag, 0>);
}
}

#endif
