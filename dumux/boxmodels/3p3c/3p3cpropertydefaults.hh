// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
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
 * \ingroup ThreePThreeCModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        3p3c box model.
 */
#ifndef DUMUX_3P3C_PROPERTY_DEFAULTS_HH
#define DUMUX_3P3C_PROPERTY_DEFAULTS_HH

#include "3p3cmodel.hh"
#include "3p3cproblem.hh"
#include "3p3cnewtoncontroller.hh"
#include "3p3cprimaryvariables.hh"
#include "3p3cratevector.hh"
#include "3p3cfluxvariables.hh"
#include "3p3cvolumevariables.hh"
#include "3p3cindices.hh"
#include "3p3cproperties.hh"

#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxThreePThreeC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    
public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 3,
                  "Only fluid systems with 3 components are supported by the 3p3c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxThreePThreeC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    
public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p3c model!");
};

/*!
 * \brief Set the number of PDEs to 3
 */
SET_INT_PROP(BoxThreePThreeC, NumEq, 3);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxThreePThreeC, 
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BoxThreePThreeC,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxThreePThreeC,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! The local residual function of the conservation equations
SET_TYPE_PROP(BoxThreePThreeC, LocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//! Use the 3p3c specific newton controller for the 3p3c model
SET_TYPE_PROP(BoxThreePThreeC, NewtonController, ThreePThreeCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxThreePThreeC, Model, ThreePThreeCModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxThreePThreeC, BaseProblem, ThreePThreeCProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxThreePThreeC, PrimaryVariables, ThreePThreeCPrimaryVariables<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxThreePThreeC, RateVector, ThreePThreeCRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxThreePThreeC, VolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxThreePThreeC, FluxVariables, ThreePThreeCFluxVariables<TypeTag>);

//! The indices required by the isothermal 3p3c model
SET_TYPE_PROP(BoxThreePThreeC, ThreePThreeCIndices, ThreePThreeCIndices<TypeTag, /*PVOffset=*/0>);
SET_TYPE_PROP(BoxThreePThreeC, Indices, typename GET_PROP(TypeTag, ThreePThreeCIndices));

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxThreePThreeC, EnableSmoothUpwinding, false);

// set the model to medium chattyness by default
SET_INT_PROP(BoxThreePThreeC, ThreePThreeCVerbosity, 1);
}
}

#endif
