// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Dumux::NcpFluxVariables
 */
#ifndef DUMUX_NCP_FLUX_VARIABLES_HH
#define DUMUX_NCP_FLUX_VARIABLES_HH

#include "diffusion/ncpfluxvariables.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup NcpModel
 * \ingroup BoxFluxVariables
 *
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the compositional NCP model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class NcpFluxVariables
    : public BoxMultiPhaseFluxVariables<TypeTag>
    , public BoxMultiPhaseEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    typedef BoxMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };   
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef NcpFluxVariablesDiffusion<TypeTag, enableDiffusion> FluxVariablesDiffusion;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef BoxMultiPhaseEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc BoxMultiPhaseFluxVariables::update
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);

        // update the flux data of the diffusion module (i.e. with or
        // without diffusion)
        diffusionVars_.update(elemCtx, scvfIdx, timeIdx);
        
        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc BoxMultiPhaseFluxVariables::updateBoundary
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

    ////////////////////////////////////////////////
    // forward calls to the diffusion module
    //! \cond 0
    Scalar porousDiffCoeffL(int compIdx) const
    { return diffusionVars_.porousDiffCoeffL(compIdx); }

    Scalar porousDiffCoeffG(int compIIdx, int compJIdx) const
    { return diffusionVars_.porousDiffCoeffG(compIIdx, compJIdx); }

    const Scalar moleFraction(int phaseIdx, int compIdx) const
    { return diffusionVars_.moleFraction(phaseIdx, compIdx); }

    const DimVector &moleFracGrad(int phaseIdx,
                               int compIdx) const
    { return diffusionVars_.moleFracGrad(phaseIdx, compIdx); }
    //! \endcond
    // end of forward calls to the diffusion module
    ////////////////////////////////////////////////

private:
    // data for the diffusion and the energy modules
    FluxVariablesDiffusion diffusionVars_;
};

} // end namepace

#endif
