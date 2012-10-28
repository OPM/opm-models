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
 * \copydoc Dumux::PvsFluxVariables
 */
#ifndef DUMUX_PVS_FLUX_VARIABLES_HH
#define DUMUX_PVS_FLUX_VARIABLES_HH

#include "pvsproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/modules/diffusion/boxdiffusionmodule.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup PvsModel
 * \ingroup BoxFluxVariables
 *
 * \brief Contains all data which is required to calculate all fluxes
 *        of components over a face of a finite volume for the
 *        compositional multi-phase primary variable switching box
 *        model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class PvsFluxVariables
    : public BoxMultiPhaseFluxVariables<TypeTag>
    , public BoxMultiPhaseEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
    , public BoxDiffusionFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion)>
{
    typedef BoxMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef BoxDiffusionFluxVariables<TypeTag, enableDiffusion> DiffusionFluxVariables;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef BoxMultiPhaseEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc BoxMultiPhaseFluxVariables::update
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);
        DiffusionFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc BoxMultiPhaseFluxVariables::updateBoundary
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

} // end namepace

#endif
