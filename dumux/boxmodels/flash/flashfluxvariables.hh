// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::FlashFluxVariables
 */
#ifndef DUMUX_FLASH_FLUX_VARIABLES_HH
#define DUMUX_FLASH_FLUX_VARIABLES_HH

#include "flashproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup FlashModel
 * \ingroup BoxFluxVariables
 *
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the compositional multi-phase model based on
 *        flash calculations.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class FlashFluxVariables
    : public BoxMultiPhaseFluxVariables<TypeTag>
    , public BoxMultiPhaseEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    typedef BoxMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef BoxMultiPhaseEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \copydoc BoxMultiPhaseFluxVariables::update
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);

        calculateGradients_(elemCtx, scvfIdx, timeIdx);

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

    /*!
     * \brief Returns the effective molecular diffusion coefficient of
     *        a component in the porous medium at the integration point.
     *
     * \copydetails Doxygen::phaseIdxParam
     * \copydetails Doxygen::compIdxParam
     */
    Scalar porousDiffCoeff(int phaseIdx, int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return porousDiffCoeff_[phaseIdx];
    }

    /*!
     * \brief Returns the molar density of a phase at the integration
     *        point [mol/m^3].
     *
     * \copydetails Doxygen::phaseIdxParam
     */
    Scalar molarDensity(int phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return molarDensity_[phaseIdx];
    }

    /*!
     * \brief Returns the mole fraction of a component in a phase at
     *        the integration point [-].
     *
     * \copydetails Doxygen::phaseIdxParam
     * \copydetails Doxygen::compIdxParam
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return moleFrac_[phaseIdx][compIdx];
    }

    /*!
     * \brief Returns the gradient of the mole fraction of a component
     *        in a phase at the integration point [1/m].
     *
     * \copydetails Doxygen::phaseIdxParam
     * \copydetails Doxygen::compIdxParam
     */
    const DimVector &moleFracGrad(int phaseIdx, int compIdx) const
    { return moleFracGrad_[phaseIdx][compIdx]; }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx,
                             int timeIdx)
    {
        // reset all quantities to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            molarDensity_[phaseIdx] = Scalar(0);
            porousDiffCoeff_[phaseIdx] = Scalar(0);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFracGrad_[phaseIdx][compIdx] = Scalar(0);
                moleFrac_[phaseIdx][compIdx] = Scalar(0);
            }
        }

        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];

        // calculate gradients
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx)
        {
            // FE gradient at vertex
            const DimVector &feGrad = scvf.grad[scvIdx];
            Scalar shapeValue = scvf.shapeValue[scvIdx];
            const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
            const auto &fluidState = volVars.fluidState();
            DimVector tmp;

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                molarDensity_[phaseIdx] += shapeValue * fluidState.molarDensity(phaseIdx);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    tmp = feGrad;
                    tmp *= fluidState.moleFraction(phaseIdx, compIdx);
                    moleFracGrad_[phaseIdx][compIdx] += tmp;

                    moleFrac_[phaseIdx][compIdx] +=
                        shapeValue * fluidState.moleFraction(phaseIdx, compIdx);
                }
            }

            tmp = feGrad;
        }

#if 0
        const auto &volVarsI = elemCtx.volVars(this->insideIdx(), timeIdx);
        const auto &volVarsJ = elemCtx.volVars(this->outsideIdx(), timeIdx);
        const auto &fsI = volVarsI.fluidState();
        const auto &fsJ = volVarsJ.fluidState();

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient
            Scalar tauI =
                1.0/(volVarsI.porosity() * volVarsI.porosity())
                *
                std::pow(std::max(0.0, volVarsI.porosity() * fsI.saturation(phaseIdx)),
                         7.0/3);
            Scalar tauJ =
                1.0/(volVarsJ.porosity() * volVarsJ.porosity())
                *
                std::pow(std::max(0.0, volVarsJ.porosity() * fsJ.saturation(phaseIdx)),
                         7.0/3);

            // Diffusion coefficient in the porous medium
            // -> harmonic mean
            porousDiffCoeff_[phaseIdx] =
                harmonicMean(volVarsI.porosity()
                             * fsI.saturation(phaseIdx)
                             * tauI
                             * volVarsI.diffCoeff(phaseIdx)
                             ,
                             volVarsJ.porosity()
                             * fsJ.saturation(phaseIdx)
                             * tauJ
                             * volVarsJ.diffCoeff(phaseIdx));
        }
#endif
    }

    Scalar porousDiffCoeff_[numPhases];
    Scalar molarDensity_[numPhases];
    Scalar moleFrac_[numPhases][numComponents];
    DimVector moleFracGrad_[numPhases][numComponents];
};

} // end namepace

#endif
