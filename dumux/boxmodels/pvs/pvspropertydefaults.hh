// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Melanie Darcis                                    *
 *   Copyright (C) 2011 by Bernd Flemisch                                    *
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
 * \ingroup PvsModel
 *
 * \brief Defines default values for most properties required by the
 *        PVS box model.
 */
#ifndef DUMUX_PVS_PROPERTY_DEFAULTS_HH
#define DUMUX_PVS_PROPERTY_DEFAULTS_HH

#include "pvsmodel.hh"
#include "pvsnewtoncontroller.hh"
#include "pvsprimaryvariables.hh"
#include "pvsratevector.hh"
#include "pvsboundaryratevector.hh"
#include "pvsvolumevariables.hh"
#include "pvsfluxvariables.hh"
#include "pvsindices.hh"
#include "pvsproperties.hh"

#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include <dumux/boxmodels/modules/velocity/boxvelocitymodules.hh>
#include <dumux/material/fluidmatrixinteractions/mp/nullmateriallaw.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 */
SET_INT_PROP(BoxPvs, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 */
SET_INT_PROP(BoxPvs, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the number of PDEs to the number of compontents
 */
SET_INT_PROP(BoxPvs, 
             NumEq,
             GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material law to the dummy law.
 */
SET_TYPE_PROP(BoxPvs,
              MaterialLaw, 
              Dumux::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases), typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxPvs,
              MaterialLawParams, 
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BoxPvs,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxPvs,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the PVS local jacobian operator for the PVS model
SET_TYPE_PROP(BoxPvs,
              LocalResidual,
              PvsLocalResidual<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(BoxPvs, VelocityModule, Dumux::BoxDarcyVelocityModule<TypeTag>);

//! Use the PVS specific newton controller for the PVS model
SET_TYPE_PROP(BoxPvs, NewtonController, PvsNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxPvs, Model, PvsModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxPvs, BaseProblem, BoxMultiPhaseProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxPvs, PrimaryVariables, PvsPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BoxPvs, RateVector, PvsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxPvs, BoundaryRateVector, PvsBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxPvs, VolumeVariables, PvsVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxPvs, FluxVariables, PvsFluxVariables<TypeTag>);

//! The indices required by the isothermal PVS model
SET_TYPE_PROP(BoxPvs, Indices, PvsIndices<TypeTag, /*PVIdx=*/0>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxPvs, EnableSmoothUpwinding, false);

// set the model to a medium verbosity
SET_INT_PROP(BoxPvs, PvsVerbosity, 1);

// disable the energy equation by default
SET_BOOL_PROP(BoxPvs, EnableEnergy, false);
}

}

#endif
