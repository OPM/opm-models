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
 * \copydoc Ewoms::PvsIntensiveQuantities
 */
#ifndef EWOMS_PVS_INTENSIVE_QUANTITIES_HH
#define EWOMS_PVS_INTENSIVE_QUANTITIES_HH

#include "pvsproperties.hh"

#include <ewoms/models/common/energymodule.hh>
#include <ewoms/models/common/diffusionmodule.hh>

#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/constraintsolvers/MiscibleMultiPhaseComposition.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>

namespace Ewoms {
/*!
 * \ingroup PvsModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the compositional multi-phase primary
 *        variable switching model.
 */
template <class TypeTag>
class PvsIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public DiffusionIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
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

    typedef Opm::MiscibleMultiPhaseComposition<Scalar, FluidSystem>
    MiscibleMultiPhaseComposition;
    typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem>
    ComputeFromReferencePhase;

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename VelocityModule::VelocityIntensiveQuantities VelocityIntensiveQuantities;
    typedef Ewoms::DiffusionIntensiveQuantities<TypeTag, enableDiffusion> DiffusionIntensiveQuantities;
    typedef Ewoms::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;

public:
    //! The type of the object returned by the fluidState() method
    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::update
     */
    void update(const ElementContext &elemCtx,
                int dofIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx,
                           dofIdx,
                           timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto &problem = elemCtx.problem();

        /////////////
        // set the saturations
        /////////////
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setSaturation(phaseIdx,
                                      priVars.explicitSaturationValue(phaseIdx));
            Valgrind::CheckDefined(fluidState_.saturation(phaseIdx));
            sumSat += fluidState_.saturation(phaseIdx);
        }
        Valgrind::CheckDefined(priVars.implicitSaturationIdx());
        Valgrind::CheckDefined(sumSat);
        fluidState_.setSaturation(priVars.implicitSaturationIdx(), 1.0 - sumSat);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams &materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // set the absolute phase pressures in the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, priVars[pressure0Idx]
                                              + (pC[phaseIdx] - pC[0]));

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;
        int lowestPresentPhaseIdx = priVars.lowestPresentPhaseIdx();
        int numNonPresentPhases = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!priVars.phaseIsPresent(phaseIdx))
                ++numNonPresentPhases;
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
            ComputeFromReferencePhase::solve(fluidState_, paramCache,
                                             lowestPresentPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setEnthalpy=*/false);
        }
        else {
            // create the auxiliary constraints
            int numAuxConstraints = numComponents + numNonPresentPhases
                                    - numPhases;
            Opm::MMPCAuxConstraint<Scalar> auxConstraints[numComponents];

            int auxIdx = 0;
            int switchIdx = 0;
            for (; switchIdx < numPhases - 1; ++switchIdx) {
                int compIdx = switchIdx + 1;
                int switchPhaseIdx = switchIdx;
                if (switchIdx >= lowestPresentPhaseIdx)
                    switchPhaseIdx += 1;

                if (!priVars.phaseIsPresent(switchPhaseIdx)) {
                    auxConstraints[auxIdx].set(lowestPresentPhaseIdx,
                                               compIdx,
                                               priVars[switch0Idx + switchIdx]);
                    ++auxIdx;
                }
            }

            for (; auxIdx < numAuxConstraints; ++auxIdx, ++switchIdx) {
                int compIdx = numPhases - numNonPresentPhases + auxIdx;
                auxConstraints[auxIdx].set(lowestPresentPhaseIdx, compIdx,
                                           priVars[switch0Idx + switchIdx]);
            }

            // both phases are present, i.e. phase compositions are a
            // result of the the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_, paramCache,
                                                 priVars.phasePresence(),
                                                 auxConstraints,
                                                 numAuxConstraints,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);
        }

#ifndef NDEBUG
        // make valgrind happy and set the enthalpies to NaN
        if (!enableEnergy) {
            Scalar myNan = std::numeric_limits<Scalar>::quiet_NaN();
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState_.setEnthalpy(phaseIdx, myNan);
        }
#endif

        /////////////
        // calculate the remaining quantities
        /////////////

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        VelocityIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        fluidState_.checkDefined();
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    Scalar mobility(int phaseIdx) const
    {
        return relativePermeability(phaseIdx) / fluidState().viscosity(phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    Scalar porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Scalar porosity_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
};

} // namespace Ewoms

#endif
