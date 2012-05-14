// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \ingroup OnePTwoCBoxModel
 * \file
 *
 * \brief Defines some default values for the properties of the the
 *        single-phase, two-component BOX model.
 */

#ifndef DUMUX_1P2C_PROPERTY_DEFAULTS_HH
#define DUMUX_1P2C_PROPERTY_DEFAULTS_HH

#include "1p2cproperties.hh"

#include "1p2cmodel.hh"
#include "1p2clocalresidual.hh"
#include "1p2cratevector.hh"
#include "1p2cboundaryratevector.hh"
#include "1p2cprimaryvariables.hh"
#include "1p2cvolumevariables.hh"
#include "1p2cfluxvariables.hh"
#include "1p2cindices.hh"

#include <dumux/material/fluidmatrixinteractions/dummymateriallaw.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>

namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxOnePTwoC, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(BoxOnePTwoC, NumPhases, 1); //!< The number of phases in the 1p2c model is 1
SET_INT_PROP(BoxOnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2

//! Use the 1p2c local residual function for the 1p2c model
SET_TYPE_PROP(BoxOnePTwoC, LocalResidual, OnePTwoCLocalResidual<TypeTag>);


//! set the material law to the dummy law
SET_TYPE_PROP(BoxOnePTwoC,
              MaterialLaw,
              Dumux::DummyMaterialLaw</*numPhases=*/1, typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the material law
//! from the law itself
SET_TYPE_PROP(BoxOnePTwoC,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to the dummy law
SET_TYPE_PROP(BoxOnePTwoC,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxOnePTwoC,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! define the model
SET_TYPE_PROP(BoxOnePTwoC, Model, OnePTwoCBoxModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxOnePTwoC, BaseProblem, OnePTwoCBoxProblem<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BoxOnePTwoC, RateVector, OnePTwoCRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxOnePTwoC, BoundaryRateVector, OnePTwoCBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxOnePTwoC, PrimaryVariables, OnePTwoCPrimaryVariables<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(BoxOnePTwoC, VolumeVariables, OnePTwoCVolumeVariables<TypeTag>);

//! define the FluxVariables
SET_TYPE_PROP(BoxOnePTwoC, FluxVariables, OnePTwoCFluxVariables<TypeTag>);

//! Set the indices used by the 1p2c model
SET_TYPE_PROP(BoxOnePTwoC, Indices, Dumux::OnePTwoCIndices<0>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxOnePTwoC, EnableSmoothUpwinding, false);

// disable output of a few quantities which make sense in a
// multip-hase but not in a single-phase context
SET_BOOL_PROP(BoxOnePTwoC, VtkWriteSaturations, false);
SET_BOOL_PROP(BoxOnePTwoC, VtkWriteMobilities, false);
SET_BOOL_PROP(BoxOnePTwoC, VtkWriteRelativePermeabilities, false);

// enable filter velocity output by default
SET_BOOL_PROP(BoxOnePTwoC, VtkWriteFilterVelocities, true);

}
// \}
}

#endif

