// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Bernd Flemisch                                    *
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
#ifndef DUMUX_NCP_LOCAL_RESIDUAL_HH
#define DUMUX_NCP_LOCAL_RESIDUAL_HH

#include "ncpfluxvariables.hh"
#include "diffusion/ncpdiffusion.hh"

#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>
#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/common/math.hh>

namespace Dumux
{
/*!
 * \ingroup NcpModel
 * \ingroup BoxLocalResidual
 * \brief Compositional NCP-model specific details needed to
 *        approximately calculate the local defect in the box scheme.
 *
 * This class is used to fill the gaps in BoxLocalResidual for
 * M-phase, N-component flow using NCPs as the model equations.
 */
template<class TypeTag>
class NcpLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    friend class BoxLocalResidual<TypeTag>;

protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        numEq = GET_PROP_VALUE(TypeTag, NumEq),

        enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy),
        ncp0EqIdx = Indices::ncp0EqIdx,
        conti0EqIdx = Indices::conti0EqIdx
    };

    typedef Dumux::BoxConstraintsContext<TypeTag> ConstraintsContext;
    typedef Dune::BlockVector<EqVector> LocalBlockVector;
    typedef BoxMultiPhaseEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Evaluate the storage term [kg/m^3] in a single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int scvIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();
        
        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            int eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                fs.molarity(phaseIdx, compIdx)
                * fs.saturation(phaseIdx)
                * volVars.porosity();
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.volVars(scvIdx, timeIdx), phaseIdx);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The number of moles of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        storage = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, scvIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidHeatStorage(storage, elemCtx.volVars(scvIdx, timeIdx));
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0.0;
        computeAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        computeDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the finite volume of the
            // current phase
            const VolumeVariables &up = elemCtx.volVars(evalPointFluxVars.upstreamIndex(phaseIdx), timeIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] +=
                    fluxVars.volumeFlux(phaseIdx)
                    * fluxVars.upstreamWeight(phaseIdx)
                    * up.fluidState().molarity(phaseIdx, compIdx);

                Valgrind::CheckDefined(flux[eqIdx]);
            }
        }
        
        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current sub control volume face
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
#if 0
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;

        // add diffusive flux of gas component in liquid phase
        Scalar tmp = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int compIdx = 1; // HACK
            tmp = - (fluxVars.moleFracGrad(phaseIdx, compIdx)*normal);
            tmp *=
                fluxVars.porousDiffCoeff(phaseIdx, compIdx) *
                fluxVars.molarDensity(phaseIdx);

            flux[conti0EqIdx + compIdx] += tmp;
            flux[conti0EqIdx + (1 - compIdx)] -= tmp;
        }
#endif

        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Set the values of the constraint volumes of the current element.
     */
    void evalConstraints_(LocalBlockVector &residual,
                          LocalBlockVector &storageTerm,
                          const ElementContext &elemCtx,
                          int timeIdx) const
    {
        // set the auxiliary functions
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                residual[scvIdx][ncp0EqIdx + phaseIdx] =
                    phaseNcp_(elemCtx, scvIdx, timeIdx, phaseIdx);
            }
        }

        // overwrite the constraint equations
        ParentType::evalConstraints_(residual, storageTerm, elemCtx, timeIdx);
    }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     */
    Scalar phaseNcp_(const ElementContext &elemCtx,
                     int scvIdx,
                     int timeIdx,
                     int phaseIdx) const
    {
        const auto &fsEval = elemCtx.evalPointVolVars(scvIdx, timeIdx).fluidState();
        const auto &fs = elemCtx.volVars(scvIdx, timeIdx).fluidState();

        Scalar aEval = phaseNotPresentIneq_(fsEval, phaseIdx);
        Scalar bEval = phasePresentIneq_(fsEval, phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq_(fs, phaseIdx);
        return phaseNotPresentIneq_(fs, phaseIdx);
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

} // end namepace

#endif
