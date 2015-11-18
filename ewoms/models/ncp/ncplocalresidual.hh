// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    typedef Ewoms::DiffusionModule<TypeTag, enableDiffusion> DiffusionModule;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

    typedef Dune::FieldVector<Evaluation, numEq> EvalEqVector;
    typedef Dune::BlockVector<EvalEqVector> ElemEvalEqVector;
    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::addPhaseStorage
     */
    template <class LhsEval>
    void addPhaseStorage(Dune::FieldVector<LhsEval, numEq> &storage,
                         const ElementContext &elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx,
                         unsigned phaseIdx) const
    {
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto &fluidState = intQuants.fluidState();

        // compute storage term of all components within all phases
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            unsigned eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                Toolbox::template toLhs<LhsEval>(fluidState.molarity(phaseIdx, compIdx))
                * Toolbox::template toLhs<LhsEval>(fluidState.saturation(phaseIdx))
                * Toolbox::template toLhs<LhsEval>(intQuants.porosity());
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx), phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext &elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        storage = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, dofIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidHeatStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx));
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux, const ElementContext &elemCtx,
                     unsigned scvfIdx, unsigned timeIdx) const
    {
        flux = Toolbox::createConstant(0.0);
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          unsigned scvfIdx, unsigned timeIdx) const
    {
        const auto &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned interiorIdx = extQuants.interiorIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream and the downstream DOFs
            // of the current phase
            unsigned upIdx = extQuants.upstreamIndex(phaseIdx);
            const IntensiveQuantities &up = elemCtx.intensiveQuantities(upIdx, timeIdx);

            // this is a bit hacky because it is specific to the element-centered
            // finite volume scheme. (N.B. that if finite differences are used to
            // linearize the system of equations, it does not matter.)
            if (upIdx == interiorIdx) {
                Evaluation tmp =
                    up.fluidState().molarDensity(phaseIdx)
                    * extQuants.volumeFlux(phaseIdx);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*up.fluidState().moleFraction(phaseIdx, compIdx);
                }
            }
            else {
                Evaluation tmp =
                    Toolbox::value(up.fluidState().molarDensity(phaseIdx))
                    * extQuants.volumeFlux(phaseIdx);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    flux[conti0EqIdx + compIdx] +=
                        tmp*Toolbox::value(up.fluidState().moleFraction(phaseIdx, compIdx));
                }
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          unsigned scvfIdx, unsigned timeIdx) const
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
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }


    /*!
     * \brief Set the values of the constraint volumes of the current element.
     */
    void evalConstraints_(ElemEvalEqVector& residual,
                          ElemEvalEqVector& storageTerm,
                          const ElementContext& elemCtx, unsigned timeIdx) const
    {
        // set the auxiliary functions
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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
    Evaluation phaseNcp_(const ElementContext &elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx,
                         unsigned phaseIdx) const
    {
        const auto &fluidState = elemCtx.intensiveQuantities(dofIdx, timeIdx).fluidState();

        const Evaluation& a = phaseNotPresentIneq_(fluidState, phaseIdx);
        const Evaluation& b = phasePresentIneq_(fluidState, phaseIdx);
        return Toolbox::min(a, b);
    }

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     */
    template <class FluidState>
    Evaluation phasePresentIneq_(const FluidState &fluidState, unsigned phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     */
    template <class FluidState>
    Evaluation phaseNotPresentIneq_(const FluidState &fluidState, unsigned phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Evaluation a = 1;
        for (unsigned i = 0; i < numComponents; ++i)
            a -= fluidState.moleFraction(phaseIdx, i);
        return a;
    }
};

} // namespace Ewoms

#endif
