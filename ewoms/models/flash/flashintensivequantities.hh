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
 * \copydoc Ewoms::FlashIntensiveQuantities
 */
#ifndef EWOMS_FLASH_INTENSIVE_QUANTITIES_HH
#define EWOMS_FLASH_INTENSIVE_QUANTITIES_HH

#include "flashproperties.hh"
#include "flashindices.hh"

#include <ewoms/models/common/energymodule.hh>
#include <ewoms/models/common/diffusionmodule.hh>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the intensive quantities of the flash-based compositional multi-phase model
 */
template <class TypeTag>
class FlashIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public DiffusionIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;

    // primary variable indices
    enum { cTot0Idx = Indices::cTot0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FlashSolver) FlashSolver;

    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Ewoms::DiffusionIntensiveQuantities<TypeTag, enableDiffusion> DiffusionIntensiveQuantities;
    typedef Ewoms::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;

public:
    //! The type of the object returned by the fluidState() method
    typedef Opm::CompositionalFluidState<Evaluation, FluidSystem, enableEnergy> FluidState;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    {
        ParentType::update(elemCtx,
                           dofIdx,
                           timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto &problem = elemCtx.problem();
        Scalar flashTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, FlashTolerance);

        // extract the total molar densities of the components
        ComponentVector cTotal;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            cTotal[compIdx] = priVars.makeEvaluation(cTot0Idx + compIdx, timeIdx);

        typename FluidSystem::ParameterCache paramCache;
        const auto *hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        if (hint) {
            // use the same fluid state as the one of the hint, but
            // make sure that we don't overwrite the temperature
            // specified by the primary variables
            const Evaluation& T = fluidState_.temperature(/*phaseIdx=*/0);
            fluidState_.assign(hint->fluidState());
            fluidState_.setTemperature(T);
        }
        else
            FlashSolver::guessInitial(fluidState_, paramCache, cTotal);

        // compute the phase compositions, densities and pressures
        const MaterialLawParams &materialParams =
            problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        FlashSolver::template solve<MaterialLaw>(fluidState_,
                                                 paramCache,
                                                 materialParams,
                                                 cTotal,
                                                 flashTolerance);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // set the phase viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        /////////////
        // calculate the remaining quantities
        /////////////

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);
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
    const Evaluation& relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Ewoms

#endif
