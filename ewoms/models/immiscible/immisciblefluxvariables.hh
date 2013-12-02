// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \copydoc Ewoms::ImmiscibleFluxVariables
 */
#ifndef EWOMS_IMMISCIBLE_FLUX_VARIABLES_HH
#define EWOMS_IMMISCIBLE_FLUX_VARIABLES_HH

#include "immiscibleproperties.hh"

#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/disc/vcfv/vcfvmultiphasefluxvariables.hh>

namespace Ewoms {

/*!
 * \ingroup ImmiscibleVcfvModel
 * \ingroup VCFVFluxVariables
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
    : public VcfvMultiPhaseFluxVariables<TypeTag>,
      public VcfvEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef VcfvMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;
    typedef VcfvEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::registerParameters()
     */
    static void registerParameters()
    { MultiPhaseFluxVariables::registerParameters(); }

    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::updateBoundary()
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context, int bfIdx, int timeIdx,
                        const FluidState &fluidState,
                        typename FluidSystem::ParameterCache &paramCache)
    {
        MultiPhaseFluxVariables::updateBoundary(context, bfIdx, timeIdx,
                                                fluidState, paramCache);
        EnergyFluxVariables::updateBoundary_(context, bfIdx, timeIdx, fluidState);
    }
};

} // namespace Ewoms

#endif
