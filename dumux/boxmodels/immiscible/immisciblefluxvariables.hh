// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Bernd Flemisch                                    *
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
 * \copydoc Dumux::ImmiscibleFluxVariables
 */
#ifndef DUMUX_IMMISCIBLE_FLUX_VARIABLES_HH
#define DUMUX_IMMISCIBLE_FLUX_VARIABLES_HH

#include "immiscibleproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup ImmiscibleBoxModel
 * \ingroup BoxFluxVariables
 *
 * \brief This class provides the data all quantities that are required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the immiscible multi-phase model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class ImmiscibleFluxVariables
    : public BoxMultiPhaseFluxVariables<TypeTag>
    , public BoxMultiPhaseEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef BoxMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;
    typedef BoxMultiPhaseEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc BoxMultiPhaseFluxVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc BoxMultiPhaseFluxVariables::updateBoundary()
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context, 
                        int bfIdx, 
                        int timeIdx, 
                        const FluidState &fs, 
                        typename FluidSystem::ParameterCache &paramCache)
    {
        MultiPhaseFluxVariables::updateBoundary(context, 
                                                bfIdx, 
                                                timeIdx, 
                                                fs, 
                                                paramCache);
        EnergyFluxVariables::updateBoundary_(context, bfIdx, timeIdx, fs);
    }
};

} // end namepace

#endif
