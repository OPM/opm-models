/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2011 by Philipp Nuske

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
 * \copydoc Ewoms::NcpLocalResidual
 */
#ifndef EWOMS_NCP_LOCAL_RESIDUAL_HH
#define EWOMS_NCP_LOCAL_RESIDUAL_HH

#include "ncpproperties.hh"

#include <ewoms/models/common/diffusionmodule.hh>
#include <ewoms/models/common/energymodule.hh>

namespace Ewoms {
/*!
 * \ingroup NcpModel
 *
 * \brief Details needed to calculate the local residual in the
 *        compositional multi-phase NCP-model .
 */
template <class TypeTag>
class NcpLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef Ewoms::DiffusionModule<TypeTag, enableDiffusion> DiffusionModule;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

    typedef Dune::BlockVector<EqVector> LocalBlockVector;

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::addPhaseStorage
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int dofIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto &fluidState = intQuants.fluidState();

        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            int eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                fluidState.molarity(phaseIdx, compIdx)
                * fluidState.saturation(phaseIdx)
                * intQuants.porosity();
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx), phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeStorage
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        storage = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, dofIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidHeatStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx));
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux, const ElementContext &elemCtx,
                     int scvfIdx, int timeIdx) const
    {
        flux = 0.0;
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        const auto &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const auto &evalPointExtQuants = elemCtx.evalPointExtensiveQuantities(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream and the finite volume of the current phase
            const IntensiveQuantities &up =
                elemCtx.intensiveQuantities(evalPointExtQuants.upstreamIndex(phaseIdx), timeIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] +=
                    extQuants.volumeFlux(phaseIdx) * up.fluidState().molarity(phaseIdx, compIdx);

                Valgrind::CheckDefined(flux[eqIdx]);
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        DiffusionModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     *
     * By default, this method only asks the problem to specify a
     * source term.
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }


    /*!
     * \brief Set the values of the constraint volumes of the current element.
     */
    void evalConstraints_(LocalBlockVector &residual,
                          LocalBlockVector &storageTerm,
                          const ElementContext &elemCtx, int timeIdx) const
    {
        // set the auxiliary functions
        for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                residual[dofIdx][ncp0EqIdx + phaseIdx] =
                    phaseNcp_(elemCtx, dofIdx, timeIdx, phaseIdx);
            }
        }

        // overwrite the constraint equations
        ParentType::evalConstraints_(residual, storageTerm, elemCtx, timeIdx);
    }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     */
    Scalar phaseNcp_(const ElementContext &elemCtx,
                     int dofIdx,
                     int timeIdx,
                     int phaseIdx) const
    {
        const auto &fluidStateEval = elemCtx.evalPointIntensiveQuantities(dofIdx, timeIdx).fluidState();
        const auto &fluidState = elemCtx.intensiveQuantities(dofIdx, timeIdx).fluidState();

        Scalar aEval = phaseNotPresentIneq_(fluidStateEval, phaseIdx);
        Scalar bEval = phasePresentIneq_(fluidStateEval, phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq_(fluidState, phaseIdx);
        return phaseNotPresentIneq_(fluidState, phaseIdx);
    }

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     */
    template <class FluidState>
    Scalar phasePresentIneq_(const FluidState &fluidState, int phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     */
    template <class FluidState>
    Scalar phaseNotPresentIneq_(const FluidState &fluidState, int phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Scalar a = 1;
        for (int i = 0; i < numComponents; ++i)
            a -= fluidState.moleFraction(phaseIdx, i);
        return a;
    }
};

} // namespace Ewoms

#endif
