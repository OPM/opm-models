// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
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
 * \copydoc Dumux::RichardsFluxVariables
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#include "richardsproperties.hh"

#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup BoxFluxVariables
 *
 * \brief Calculates and stores the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsFluxVariables : public BoxMultiPhaseFluxVariables<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    /*!
     * \copydoc BoxMultiPhaseFluxVariables::usePhase
     */
    bool usePhase(int phaseIdx)
    { return phaseIdx == GET_PROP_VALUE(TypeTag, LiquidPhaseIndex); }
};

} // end namepace

#endif
