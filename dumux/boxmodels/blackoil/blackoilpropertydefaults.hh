// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2010-2011 by Markus Wolff                                 *
 *   Copyright (C) 2009-2011 by Melanie Darcis                               *
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
 *
 * \ingroup BlackOilBoxModel
 * \ingroup BoxProperties
 * \ingroup Properties
 *
 * \brief Defines default values for the properties required by the
 *        twophase box model.
 */
#ifndef DUMUX_BLACK_OIL_PROPERTY_DEFAULTS_HH
#define DUMUX_BLACK_OIL_PROPERTY_DEFAULTS_HH

#include "blackoilmodel.hh"
#include "blackoilindices.hh"
#include "blackoilfluxvariables.hh"
#include "blackoilprimaryvariables.hh"
#include "blackoilvolumevariables.hh"
#include "blackoilratevector.hh"
#include "blackoilboundaryratevector.hh"
#include "blackoilproperties.hh"

#include <dumux/boxmodels/modules/velocity/boxvelocitymodules.hh>
#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include <dumux/material/fluidsystems/blackoilfluidsystem.hh>
#include <dumux/material/heatconduction/dummyheatconductionlaw.hh>

namespace Dumux {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxBlackOil, NumEq, 3); //!< set the number of equations to 3
SET_INT_PROP(BoxBlackOil, NumPhases, 3); //!< The number of phases considered by the model
SET_INT_PROP(BoxBlackOil, NumComponents, 3);   //!< Number of chemical species in the system

//! Set the local residual function
SET_TYPE_PROP(BoxBlackOil,
              LocalResidual,
              BlackOilLocalResidual<TypeTag>);

//! The Model property
SET_TYPE_PROP(BoxBlackOil, Model, BlackOilModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxBlackOil, BaseProblem, BoxMultiPhaseProblem<TypeTag>);

//! The BlackOilFluidState property
SET_PROP(BoxBlackOil, BlackOilFluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::CompositionalFluidState<Scalar,
                                           FluidSystem,
                                           /*enableEnthalpy=*/false> type;
};

//! Use the Darcy relation by default
SET_TYPE_PROP(BoxBlackOil, VelocityModule, Dumux::BoxDarcyVelocityModule<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(BoxBlackOil, RateVector, BlackOilRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxBlackOil, BoundaryRateVector, BlackOilBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BoxBlackOil, PrimaryVariables, BlackOilPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxBlackOil, VolumeVariables, BlackOilVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxBlackOil, FluxVariables, BlackOilFluxVariables<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(BoxBlackOil, Indices, BlackOilIndices</*PVOffset=*/0>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxBlackOil,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(BoxBlackOil,
              HeatConductionLaw,
              Dumux::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxBlackOil,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Set the fluid system to the black-oil fluid system by default
SET_TYPE_PROP(BoxBlackOil, 
              FluidSystem, 
              Dumux::FluidSystems::BlackOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(BoxBlackOil, EnableSmoothUpwinding, false);

}

}

#endif
