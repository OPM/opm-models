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
 * \ingroup PvsNIModel
 * \file
 *
 * \brief Defines default values for most properties required by the 2p2cni
 *        box model.
 */
#ifndef DUMUX_PVS_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_PVS_NI_PROPERTY_DEFAULTS_HH

#include <dumux/boxmodels/pvs/pvspropertydefaults.hh>

#include "pvsnimodel.hh"
#include <dumux/boxmodels/common/boxmultiphaseproblem.hh>
#include "pvsniindices.hh"
#include "pvsniboundaryratevector.hh"
#include "pvsnilocalresidual.hh"
#include "pvsnivolumevariables.hh"
#include "pvsnifluxvariables.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxPvsNI, NumEq, GET_PROP_VALUE(TypeTag, NumComponents) + 1); //!< set the number of equations

//! Use the local jacobian operator for the non-isothermal primary variable switching model
SET_TYPE_PROP(BoxPvsNI,
              LocalResidual,
              PvsNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxPvsNI, Model, PvsNIModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(BoxPvsNI, BaseProblem, BoxMultiPhaseProblem<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxPvsNI, VolumeVariables, PvsNIVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxPvsNI, FluxVariables, PvsNIFluxVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxPvsNI, BoundaryRateVector, Dumux::PvsNIBoundaryRateVector<TypeTag>);

//! The indices required by the non-isothermal primary variable switching model
SET_TYPE_PROP(BoxPvsNI, Indices, PvsNIIndices<TypeTag, 0>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(BoxPvsNI,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

}

}
#endif
