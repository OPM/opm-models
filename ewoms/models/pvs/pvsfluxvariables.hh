/*
  Copyright (C) 2009-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::PvsFluxVariables
 */
#ifndef EWOMS_PVS_FLUX_VARIABLES_HH
#define EWOMS_PVS_FLUX_VARIABLES_HH

#include "pvsproperties.hh"

#include <ewoms/models/common/multiphasebasefluxvariables.hh>
#include <ewoms/models/common/energymodule.hh>
#include <ewoms/models/common/diffusionmodule.hh>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 * \ingroup FluxVariables
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
    : public MultiPhaseBaseFluxVariables<TypeTag>
    , public EnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
    , public DiffusionFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion)>
{
    typedef MultiPhaseBaseFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef Ewoms::DiffusionFluxVariables<TypeTag, enableDiffusion> DiffusionFluxVariables;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc ParentType::update
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);
        DiffusionFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ParentType::updateBoundary
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context, int bfIdx, int timeIdx,
                        const FluidState &fluidState,
                        typename FluidSystem::ParameterCache &paramCache)
    {
        ParentType::updateBoundary(context, bfIdx, timeIdx,
                                                fluidState, paramCache);
        DiffusionFluxVariables::updateBoundary_(context, bfIdx, timeIdx,
                                                fluidState);
        EnergyFluxVariables::updateBoundary_(context, bfIdx, timeIdx, fluidState);
    }
};

} // namespace Ewoms

#endif
