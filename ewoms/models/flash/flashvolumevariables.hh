// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2013 by Andreas Lauser                               *
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
 * \copydoc Ewoms::FlashVolumeVariables
 */
#ifndef EWOMS_FLASH_VOLUME_VARIABLES_HH
#define EWOMS_FLASH_VOLUME_VARIABLES_HH

#include "flashproperties.hh"
#include "flashindices.hh"

#include <ewoms/disc/vcfv/vcfvvolumevariables.hh>
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/models/modules/diffusion/vcfvdiffusionmodule.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/common/math.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 * \ingroup VcfvVolumeVariables
 *
 * \brief Contains the quantities which are constant within a finite
 *        volume for the flash-based compositional multi-phase model.
 */
template <class TypeTag>
class FlashVolumeVariables
    : public VcfvVolumeVariables<TypeTag>
    , public VcfvDiffusionVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion) >
    , public VcfvEnergyVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef VcfvVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // primary variable indices
    enum { cTot0Idx = Indices::cTot0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FlashSolver) FlashSolver;

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;
    typedef VcfvDiffusionVolumeVariables<TypeTag, enableDiffusion> DiffusionVolumeVariables;
    typedef VcfvEnergyVolumeVariables<TypeTag, enableEnergy> EnergyVolumeVariables;

public:
    //! The type of the object returned by the fluidState() method
    typedef Ewoms::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/enableEnergy> FluidState;

    /*!
     * \copydoc VcfvVolumeVariables::update
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
        Scalar flashTolerance = GET_PARAM(TypeTag, Scalar, FlashTolerance);
        if (flashTolerance <= 0) {
            // make the tolerance of the flash solver 10 times smaller
            // than the epsilon value used by the newton solver to
            // calculate the partial derivatives
            const auto &model = elemCtx.model();
            flashTolerance =
                model.localJacobian().baseEpsilon()
                / (100*18e-3); // assume the molar weight of water
        }

        // extract the total molar densities of the components
        ComponentVector cTotal;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            cTotal[compIdx] = priVars[cTot0Idx + compIdx];

        typename FluidSystem::ParameterCache paramCache;
        const auto *hint = elemCtx.hint(scvIdx, timeIdx);
        if (hint) {
            // use the same fluid state as the one of the hint, but
            // make sure that we don't overwrite the temperature
            // specified by the primary variables
            Scalar T = fluidState_.temperature(/*phaseIdx=*/0);
            fluidState_.assign(hint->fluidState());
            fluidState_.setTemperature(T);
        }
        else
            FlashSolver::guessInitial(fluidState_, paramCache, cTotal);

        // compute the phase compositions, densities and pressures
        const MaterialLawParams &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        FlashSolver::template solve<MaterialLaw>(fluidState_, paramCache, materialParams, cTotal, flashTolerance);

        // set the phase viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
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
