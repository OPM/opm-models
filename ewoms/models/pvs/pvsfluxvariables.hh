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
 * \copydoc Ewoms::PvsFluxVariables
 */
#ifndef EWOMS_PVS_FLUX_VARIABLES_HH
#define EWOMS_PVS_FLUX_VARIABLES_HH

#include "pvsproperties.hh"

#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/models/modules/diffusion/vcfvdiffusionmodule.hh>
#include <ewoms/disc/vcfv/vcfvmultiphasefluxvariables.hh>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 * \ingroup VCFVFluxVariables
 *
 * \brief Contains all data which is required to calculate all fluxes
 *        over a finite-volume face for the primary variable switching
 *        model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class PvsFluxVariables
    : public VcfvMultiPhaseFluxVariables<TypeTag>
    , public VcfvEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
    , public VcfvDiffusionFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion)>
{
    typedef VcfvMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;



    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef VcfvDiffusionFluxVariables<TypeTag, enableDiffusion> DiffusionFluxVariables;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef VcfvEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::update
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);
        DiffusionFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::updateBoundary
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context,
                        int bfIdx,
                        int timeIdx,
                        const FluidState &fluidState,
                        typename FluidSystem::ParameterCache &paramCache)
    {
        MultiPhaseFluxVariables::updateBoundary(context,
                                                bfIdx,
                                                timeIdx,
                                                fluidState,
                                                paramCache);
        DiffusionFluxVariables::updateBoundary_(context, bfIdx, timeIdx, fluidState);
        EnergyFluxVariables::updateBoundary_(context, bfIdx, timeIdx, fluidState);
    }
};

} // namespace Ewoms

#endif
