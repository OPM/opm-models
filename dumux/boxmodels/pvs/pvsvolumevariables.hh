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
 * \copydoc Dumux::PvsVolumeVariables
 */
#ifndef DUMUX_PVS_VOLUME_VARIABLES_HH
#define DUMUX_PVS_VOLUME_VARIABLES_HH

#include "pvsproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/modules/diffusion/boxdiffusionmodule.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/common/math.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup PvsModel
 * \ingroup BoxVolumeVariables
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase primary
 *        variable switching box model.
 */
template <class TypeTag>
class PvsVolumeVariables
    : public BoxVolumeVariables<TypeTag>
    , public BoxDiffusionVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public BoxMultiPhaseEnergyVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;

    enum { switch0Idx = Indices::switch0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;
    typedef BoxDiffusionVolumeVariables<TypeTag, enableDiffusion> DiffusionVolumeVariables;
    typedef BoxMultiPhaseEnergyVolumeVariables<TypeTag, enableEnergy> EnergyVolumeVariables;

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \copydoc ImmiscibleVolumeVariables::update
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx,
                           scvIdx,
                           timeIdx);
        EnergyVolumeVariables::updateTemperatures_(fluidState_, elemCtx, scvIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        const auto &problem = elemCtx.problem();

        /////////////
        // set the saturations
        /////////////
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fluidState_.setSaturation(phaseIdx, priVars.explicitSaturationValue(phaseIdx));
            sumSat += fluidState_.saturation(phaseIdx);
        }
        fluidState_.setSaturation(priVars.implicitSaturationIdx(), 1.0 - sumSat);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // set the absolute phase pressures in the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx,
                                   priVars[pressure0Idx] + (pC[phaseIdx] - pC[0]));

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;
        int lowestPresentPhaseIdx = priVars.lowestPresentPhaseIdx();
        int numNonPresentPhases = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!priVars.phaseIsPresent(phaseIdx))
                ++ numNonPresentPhases;
        }

        // now comes the tricky part: calculate phase compositions
        if (numNonPresentPhases == numPhases - 1) {
            // only one phase is present, i.e. the primary variables
            // contain the complete composition of the phase
            Scalar sumx = 0;
            for (int compIdx = 1; compIdx < numComponents; ++compIdx) {
                Scalar x = priVars[switch0Idx + compIdx - 1];
                fluidState_.setMoleFraction(lowestPresentPhaseIdx, compIdx, x);
                sumx += x;
            }

            // set the mole fraction of the first component
            fluidState_.setMoleFraction(lowestPresentPhaseIdx, 0, 1 - sumx);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             lowestPresentPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        else {
            // create the auxiliary constraints
            int numAuxConstraints = numComponents + numNonPresentPhases - numPhases;
            MMPCAuxConstraint<Scalar> auxConstraints[numComponents];

            int auxIdx = 0;
            int switchIdx = 0;
            for (; switchIdx < numPhases - 1; ++switchIdx) {
                int compIdx = switchIdx + 1;
                int switchPhaseIdx = switchIdx;
                if (switchIdx >= lowestPresentPhaseIdx)
                    switchPhaseIdx += 1;

                if (!priVars.phaseIsPresent(switchPhaseIdx)) {
                    auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx, priVars[switch0Idx + switchIdx]);
                    ++auxIdx;
                }
            }

            for (; auxIdx < numAuxConstraints; ++auxIdx, ++switchIdx) {
                int phaseIdx = priVars.lowestPresentPhaseIdx();
                int compIdx = numPhases - numNonPresentPhases + auxIdx;

                auxConstraints[auxIdx].set(phaseIdx, compIdx, priVars[switch0Idx + switchIdx]);
                std::cout << "auxConstraints[" << auxIdx << "] = (" << phaseIdx << ", " << compIdx << ", " << priVars[switch0Idx + switchIdx] << ")\n";
            }

            // both phases are present, i.e. phase compositions are a
            // result of the the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_,
                                                 paramCache,
                                                 priVars.phasePresence(),
                                                 auxConstraints,
                                                 numAuxConstraints,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);

        }

        /////////////
        // calculate the remaining quantities
        /////////////

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // porosity
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, scvIdx, timeIdx);

        // update the quantities specific for the velocity model
        VelocityVolumeVariables::update_(elemCtx, scvIdx, timeIdx);

        // energy related quantities
        EnergyVolumeVariables::update_(fluidState_, paramCache, elemCtx, scvIdx, timeIdx);

        // update the diffusion specific quantities of the volume variables
        DiffusionVolumeVariables::update_(fluidState_, paramCache, elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleVolumeVariables::fluidState
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::intrinsicPermeability
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::relativePermeability
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::mobility
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState().viscosity(phaseIdx); }

    /*!
     * \copydoc ImmiscibleVolumeVariables::porosity
     */
    Scalar porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Scalar porosity_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
};

} // end namepace

#endif
