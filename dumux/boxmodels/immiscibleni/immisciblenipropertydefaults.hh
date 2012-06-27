// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \ingroup ImmiscibleNIModel
 * \file
 *
 * \brief Defines default values for most properties required by the 2p2cni
 *        box model.
 */
#ifndef DUMUX_IMMISCIBLE_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_IMMISCIBLE_NI_PROPERTY_DEFAULTS_HH

#include <dumux/boxmodels/immiscible/immisciblepropertydefaults.hh>

#include "immisciblenimodel.hh"
#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include "immiscibleniindices.hh"
#include "immiscibleniboundaryratevector.hh"
#include "immisciblenilocalresidual.hh"
#include "immisciblenivolumevariables.hh"
#include "immisciblenifluxvariables.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxImmiscibleEnergy, NumEq, GET_PROP_VALUE(TypeTag, NumPhases) + 1); //!< set the number of equations

//! Use the local jacobian operator for the non-isothermal primary variable switching model
SET_TYPE_PROP(BoxImmiscibleEnergy,
              LocalResidual,
              ImmiscibleNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxImmiscibleEnergy, Model, ImmiscibleNIModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxImmiscibleEnergy, VolumeVariables, ImmiscibleNIVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxImmiscibleEnergy, FluxVariables, ImmiscibleNIFluxVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxImmiscibleEnergy, BoundaryRateVector, Dumux::ImmiscibleNIBoundaryRateVector<TypeTag>);

//! The indices required by the non-isothermal primary variable switching model
SET_TYPE_PROP(BoxImmiscibleEnergy, Indices, ImmiscibleNIIndices<TypeTag, 0>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxImmiscibleEnergy,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);


//! the FluidState property
SET_PROP(BoxImmiscibleEnergy, FluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Dumux::ImmiscibleFluidState<Scalar,
                                        FluidSystem,
                                        /*enableEnthalpy=*/true> type;
};

}

}
#endif
