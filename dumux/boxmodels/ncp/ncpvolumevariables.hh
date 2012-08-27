// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \brief Contains the secondary variables (Quantities which are
 *        constant within a finite volume) of the M-phase, N-component
 *        model.
 */
#ifndef DUMUX_NCP_VOLUME_VARIABLES_HH
#define DUMUX_NCP_VOLUME_VARIABLES_HH

#include "diffusion/ncpvolumevariables.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>
#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux {
/*!
 * \ingroup NcpModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the M-phase, N-component model.
 */
template <class TypeTag>
class NcpVolumeVariables
    : public BoxVolumeVariables<TypeTag>
    , public NcpVolumeVariablesDiffusion<TypeTag, GET_PROP_VALUE(TypeTag, EnableDiffusion)>
    , public BoxMultiPhaseEnergyVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, CompositionFromFugacitiesSolver) CompositionFromFugacitiesSolver;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { fugacity00Idx = Indices::fugacity00Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/enableEnergy> FluidState;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef NcpVolumeVariablesDiffusion<TypeTag, enableDiffusion> DiffusionVolumeVariables;
    typedef BoxMultiPhaseEnergyVolumeVariables<TypeTag, enableEnergy> EnergyVolumeVariables;
    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;

public:
    NcpVolumeVariables()
    { }

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx,
                           scvIdx,
                           timeIdx);
        ParentType::checkDefined();


        typename FluidSystem::ParameterCache paramCache;
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        /////////////
        // set the phase saturations
        /////////////
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            sumSat += priVars[saturation0Idx + phaseIdx];
            fluidState_.setSaturation(phaseIdx, priVars[saturation0Idx + phaseIdx]);
        }
        fluidState_.setSaturation(numPhases - 1, 1.0 - sumSat);
        Valgrind::CheckDefined(sumSat);


        /////////////
        // set the fluid phase temperatures
        /////////////
        EnergyVolumeVariables::updateTemperatures_(fluidState_, elemCtx, scvIdx, timeIdx);

        /////////////
        // set the phase pressures
        /////////////

        // retrieve capillary pressure parameters
        const auto &problem = elemCtx.problem();
        const MaterialLawParams &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        // calculate capillary pressures
        Scalar capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState_);
        // add to the pressure of the first fluid phase
        Scalar pressure0 = priVars[pressure0Idx];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            fluidState_.setPressure(phaseIdx, pressure0 + (capPress[phaseIdx] - capPress[0]));

        /////////////
        // set the fluid compositions
        /////////////
        const auto *hint = elemCtx.hint(scvIdx, timeIdx);

        // calculate phase compositions
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // initial guess
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar xIJ = 1.0/numComponents;
                if (hint)
                    // use the hint for the initial mole fraction!
                    xIJ = hint->fluidState().moleFraction(phaseIdx, compIdx);

                // set initial guess of the component's mole fraction
                fluidState_.setMoleFraction(phaseIdx, compIdx, xIJ);
            }

            ComponentVector fug;
            // retrieve component fugacities
            Scalar p0 = fluidState_.pressure(0);
            Scalar pAlpha = fluidState_.pressure(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fug[compIdx] = priVars[fugacity00Idx + compIdx]/p0 * pAlpha;

            // calculate the phase composition from the component
            // fugacities
            CompositionFromFugacitiesSolver::solve(fluidState_, paramCache, phaseIdx, fug);
        }

        /////////////
        // Porosity
        /////////////

        // porosity
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        /////////////
        // Phase mobilities
        /////////////

        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams,
                                            fluidState_);

        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
        }

        /////////////
        // diffusion
        /////////////

        // update the diffusion part of the volume data
        DiffusionVolumeVariables::update(fluidState_,
                                         paramCache,
                                         elemCtx,
                                         scvIdx,
                                         timeIdx);
        DiffusionVolumeVariables::checkDefined();

        /////////////
        // energy
        /////////////

        // energy related quantities
        EnergyVolumeVariables::update_(fluidState_, paramCache, elemCtx, scvIdx, timeIdx);

        fluidState_.checkDefined();

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, scvIdx, timeIdx);

        // update the quantities specific for the velocity model
        VelocityVolumeVariables::update_(elemCtx, scvIdx, timeIdx);

        checkDefined();
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the intrinsic permeability tensor for the sub-control volume
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume.
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        ParentType::checkDefined();

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(relativePermeability_);

        fluidState_.checkDefined();
#endif
    }

protected:
    FluidState fluidState_;
    Scalar porosity_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
};

} // end namepace

#endif
