// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \ingroup DiscreteFractureVcfvModel
 *
 * \brief Defines default values for the properties required by the
 *        immiscible multi-phase model which including discrete fractures.
 */
#ifndef EWOMS_DISCRETE_FRACTURE_PROPERTY_DEFAULTS_HH
#define EWOMS_DISCRETE_FRACTURE_PROPERTY_DEFAULTS_HH

#include "discretefracturemodel.hh"
#include "discretefractureprimaryvariables.hh"
#include "discretefracturevolumevariables.hh"
#include "discretefracturefluxvariables.hh"
#include "discretefracturelocalresidual.hh"
#include "discretefractureproperties.hh"

#include <ewoms/models/immiscible/immisciblepropertydefaults.hh>

namespace Ewoms {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////

//! The class for the model
SET_TYPE_PROP(VcfvDiscreteFracture, Model, DiscreteFractureModel<TypeTag>);

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(VcfvDiscreteFracture,
              LocalResidual,
              DiscreteFractureLocalResidual<TypeTag>);

//! The type of the base base class for actual problems.
// TODO!?
//SET_TYPE_PROP(VcfvDiscreteFracture, BaseProblem, DiscreteFractureBaseProblem<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvDiscreteFracture, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvDiscreteFracture, PrimaryVariables, DiscreteFracturePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(VcfvDiscreteFracture, VolumeVariables, DiscreteFractureVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(VcfvDiscreteFracture, FluxVariables, DiscreteFractureFluxVariables<TypeTag>);

//! For the discrete fracture model, we need to use two-point flux
//! appoximation or it will converge very poorly
SET_BOOL_PROP(VcfvDiscreteFracture, UseTwoPointGradients, true);

}

}

#endif
