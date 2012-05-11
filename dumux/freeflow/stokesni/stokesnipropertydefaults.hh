// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
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
 * \ingroup BoxStokesNIModel
 *
 * \file
 *
 * \brief Sets default properties for the non-isothermal compositional
 *        Stokes box model.
 */
#ifndef DUMUX_STOKES_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_STOKES_NI_PROPERTY_DEFAULTS_HH


#include "stokesnifluxvariables.hh"
#include "stokesniindices.hh"
#include "stokesnilocalresidual.hh"
#include "stokesnimodel.hh"
#include "stokesnivolumevariables.hh"
#include "stokesniboundaryratevector.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(BoxStokesNI, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static constexpr int value = GridView::dimensionworld + FluidSystem::numComponents + 1;
};

//! Use the stokesni local jacobian operator for the compositional stokes model
SET_TYPE_PROP(BoxStokesNI,
              LocalResidual,
              StokesNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxStokesNI, Model, StokesNIModel<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BoxStokesNI, BoundaryRateVector, StokesNIBoundaryRateVector<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxStokesNI, VolumeVariables, StokesNIVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxStokesNI, FluxVariables, StokesNIFluxVariables<TypeTag>);

// the indices for the StokesNI model
SET_TYPE_PROP(BoxStokesNI,
              StokesNIIndices,
              StokesNIIndices<TypeTag>);
SET_TYPE_PROP(BoxStokesNI,
              Indices,
              typename GET_PROP_TYPE(TypeTag, StokesNIIndices));
}
}
#endif
