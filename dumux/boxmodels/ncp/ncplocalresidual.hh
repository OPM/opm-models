// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
#include "energy/ncplocalresidualenergy.hh"

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

    typedef NcpLocalResidualEnergy<TypeTag, enableEnergy> EnergyResid;
    typedef Dumux::BoxConstraintsContext<TypeTag> ConstraintsContext;

    typedef Dune::BlockVector<EqVector> LocalBlockVector;

public:
    /*!
     * \brief Evaluate the amount moles within a sub-control volume in
     *        a phase.
     *
     * The result should be averaged over the volume.
     */
    void addPhaseStorageMass(EqVector &storage,
                             const ElementContext &elemCtx,
                             int scvIdx,
                             int timeIdx,
                             int phaseIdx) const
    {
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);

        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            storage[compIdx] +=
                volVars.fluidState().saturation(phaseIdx)
                * volVars.fluidState().molarity(phaseIdx, compIdx)
                * volVars.porosity();
        }        
    }

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VolumeVariables &volVars =
            elemCtx.volVars(scvIdx, timeIdx);

        storage = 0;

        // compute mass storage terms
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorageMass(storage, elemCtx, scvIdx, timeIdx, phaseIdx);
        Valgrind::CheckDefined(storage);

        // compute energy storage term
        EnergyResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
    }

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within all sub-control volumes of an
     *        element.
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int timeIdx,
                         int phaseIdx) const
    {
        // calculate the phase storage for all sub-control volumes
        for (int scvIdx=0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            EqVector tmp(0.0);
            
            // compute mass and energy storage terms in terms of
            // averaged quantities
            addPhaseStorageMass(tmp,
                                elemCtx,
                                scvIdx,
                                timeIdx,
                                phaseIdx);
            EnergyResid::addPhaseStorage(storage,
                                         elemCtx.volVars(scvIdx, /*timeIdx=*/0),
                                         phaseIdx);
            
            // multiply with volume of sub-control volume
            tmp *=
                elemCtx.volVars(scvIdx, /*timeIdx=*/0).extrusionFactor() *
                elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[scvIdx].volume;

            // Add the storage of the current SCV to the total storage
            storage += tmp;
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    void addAdvectiveMassFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx,
                              int phaseIdx) const
    {
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &evalFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);
        const auto &up = elemCtx.volVars(evalFluxVars.upstreamIdx(phaseIdx), timeIdx);

        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            flux[conti0EqIdx + compIdx] =
                up.fluidState().molarity(phaseIdx, compIdx)
                * fluxVars.filterVelocityNormal(phaseIdx);
        }
    }


    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    void addDiffusiveMassFlux(EqVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx,
                              int phaseIdx) const
    {
#if 0
        if (!enableDiffusion) {
            flux = 0.0;
            return;
        }

        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const VolumeVariables &volVarsI = elemCtx.volVars(fluxVars.insideIdx(), timeIdx);
        const VolumeVariables &volVarsJ = elemCtx.volVars(fluxVars.outsideIdx(), timeIdx);
        if (volVarsI.fluidState().saturation(phaseIdx) < 1e-4 ||
            volVarsJ.fluidState().saturation(phaseIdx) < 1e-4)
        {
            return; // phase is not present in one of the finite volumes
        }

        // approximate the total concentration of the phase at the
        // integration point by the arithmetic mean of the
        // concentration of the sub-control volumes
        Scalar molarDensity;
        molarDensity = volVarsI.fluidState().molarDensity(phaseIdx);
        molarDensity += volVarsJ.fluidState().molarDensity(phaseIdx);
        molarDensity /= 2;

        Diffusion::flux(flux, elemCtx, scvfIdx, timeIdx, phaseIdx, molarDensity);
#endif
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

        RateVector tmp;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            tmp = 0.0;
            addAdvectiveMassFlux(tmp,
                                 elemCtx,
                                 scvfIdx,
                                 timeIdx,
                                 phaseIdx);
            flux += tmp;
            Valgrind::CheckDefined(flux);
            addDiffusiveMassFlux(flux,
                                 elemCtx,
                                 scvfIdx,
                                 timeIdx,
                                 phaseIdx);
            Valgrind::CheckDefined(flux);

            // energy transport in fluid phases
            EnergyResid::addPhaseEnthalpyFlux(flux,
                                              elemCtx,
                                              scvfIdx,
                                              timeIdx,
                                              phaseIdx,
                                              tmp);
            Valgrind::CheckDefined(flux);
        }

        // energy transport in fluid phases
        EnergyResid::addHeatConduction(flux,
                                       elemCtx,
                                       scvfIdx,
                                       timeIdx);
        Valgrind::CheckDefined(flux);
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
