// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::BlackOilFluxVariables
 */
#ifndef EWOMS_BLACK_OIL_FLUX_VARIABLES_HH
#define EWOMS_BLACK_OIL_FLUX_VARIABLES_HH

#include "blackoilproperties.hh"

#include <ewoms/disc/vcfv/vcfvmultiphasefluxvariables.hh>

namespace Ewoms {

/*!
 * \ingroup BlackOilVcfvModel
 * \ingroup VCFVFluxVariables
 *
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the black-oil model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class BlackOilFluxVariables : public VcfvMultiPhaseFluxVariables<TypeTag>
{
};

} // namespace Ewoms

#endif
