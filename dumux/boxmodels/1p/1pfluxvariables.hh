// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Katherina Baber
 *   Copyright (C) 2008-2009 by Onur Dogan                                   *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        the flux of the fluid over a face of a finite volume for the one-phase model.
 *
 *        This means pressure and temperature gradients, phase densities at
 *           the integration point, etc.
 */
#ifndef DUMUX_1P_FLUX_VARIABLES_HH
#define DUMUX_1P_FLUX_VARIABLES_HH

#include "1pproperties.hh"

#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the flux of the fluid over a face of a
 *        finite volume for the one-phase model.
 */
template <class TypeTag>
class OnePFluxVariables : public BoxMultiPhaseFluxVariables<TypeTag>
{
};

} // end namepace

#endif
