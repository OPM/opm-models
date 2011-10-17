// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
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
/*!
 * \file
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include "2p2cproperties.hh"
#include "2p2cindices.hh"

#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class TwoPTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };
    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum {
    };

    // present phases
    enum {
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl
    };

    // primary variable indices
    enum {
        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension};

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int historyIdx)
    {
        ParentType::update(elemCtx,
                           scvIdx,
                           historyIdx);

        completeFluidState(fluidState_, elemCtx, scvIdx, historyIdx);

        /////////////
        // calculate the remaining quantities
        /////////////
        const auto &spatialParams = elemCtx.problem().spatialParameters();
        const MaterialLawParams &materialParams =
            spatialParams.materialLawParams(elemCtx, scvIdx);
        
        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        // energy related quantities
        asImp_().updateEnergy_(paramCache, elemCtx, scvIdx, historyIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // relative permeabilities
            Scalar kr;
            Scalar Sw = fluidState_.saturation(/*phaseIdx=*/0);
            if (phaseIdx == 0)
                kr = MaterialLaw::krw(materialParams, Sw);
            else // ATTENTION: krn requires the liquid saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, Sw);
            relativePermeability_[phaseIdx] = kr;
            Valgrind::CheckDefined(relativePermeability_[phaseIdx]);

            // binary diffusion coefficents
            diffCoeff_[phaseIdx] =
                FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                        paramCache,
                                                        phaseIdx,
                                                        /*compIdx=*/0,
                                                        /*compIdx=*/1);
            Valgrind::CheckDefined(diffCoeff_[phaseIdx]);
        }

        // porosity
        porosity_ = spatialParams.porosity(elemCtx, scvIdx);
        Valgrind::CheckDefined(porosity_);
    }

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        Implementation::updateTemperature_(fluidState, elemCtx, scvIdx, historyIdx);
        
        const auto &priVars = elemCtx.primaryVars(scvIdx, historyIdx);
        const auto &problem = elemCtx.problem();
        const auto &model = elemCtx.model();
        const auto &spatialParams = problem.spatialParameters();
        int globalVertIdx = model.dofMapper().map(elemCtx.element(), scvIdx, dim);
        int phasePresence = model.phasePresence(globalVertIdx, historyIdx);

        /////////////
        // set the saturations
        /////////////
        Scalar Sw;
        if (phasePresence == gPhaseOnly)
            Sw = 0.0;
        else if (phasePresence == lPhaseOnly) {
            Sw = 1.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == plSg)
                Sw = 1.0 - priVars[switchIdx];
            else if (formulation == pgSl)
                Sw = priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        fluidState.setSaturation(/*phaseIdx=*/0, Sw);
        fluidState.setSaturation(/*phaseIdx=*/1, 1 - Sw);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams &materialParams =
            spatialParams.materialLawParams(elemCtx, scvIdx);
        Scalar pC = MaterialLaw::pC(materialParams, Sw);

        if (formulation == plSg) {
            fluidState.setPressure(/*phaseIdx=*/0, priVars[pressureIdx]);
            fluidState.setPressure(/*phaseIdx=*/1, priVars[pressureIdx] + pC);
        }
        else if (formulation == pgSl) {
            fluidState.setPressure(/*phaseIdx=*/1, priVars[pressureIdx]);
            fluidState.setPressure(/*phaseIdx=*/0, priVars[pressureIdx] - pC);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;
        // now comes the tricky part: calculate phase compositions
        if (phasePresence == bothPhases) {
            // both phases are present, phase compositions are a
            // result of the the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState,
                                                 paramCache,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);

        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, i.e. gas phase
            // composition is stored explicitly.

            // extract _mass_ fractions in the gas phase
            Scalar massFractionG[numComponents];
            massFractionG[/*compIdx=*/0] = priVars[switchIdx];
            massFractionG[/*compIdx=*/1] = 1 - massFractionG[/*compIdx=*/0];

            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(/*compIdx=*/0);
            Scalar M2 = FluidSystem::molarMass(/*compIdx=*/1);
            Scalar X2 = massFractionG[/*compIdx=*/1];
            Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions and set the fluid state
            fluidState.setMoleFraction(/*phaseIdx=*/1, /*compIdx=*/0, massFractionG[/*compIdx=*/0]*avgMolarMass/M1);
            fluidState.setMoleFraction(/*phaseIdx=*/1, /*compIdx=*/1, massFractionG[/*compIdx=*/1]*avgMolarMass/M2);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             /*phaseIdx=*/1,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == lPhaseOnly) {
            // only the liquid phase is present, i.e. liquid phase
            // composition is stored explicitly.

            // extract _mass_ fractions in the gas phase
            Scalar massFractionL[numComponents];
            massFractionL[/*compIdx=*/1] = priVars[switchIdx];
            massFractionL[/*compIdx=*/0] = 1 - massFractionL[/*compIdx=*/1];

            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(/*compIdx=*/0);
            Scalar M2 = FluidSystem::molarMass(/*compIdx=*/1);
            Scalar X2 = massFractionL[/*compIdx=*/1];
            Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions and set the fluid state
            fluidState.setMoleFraction(/*phaseIdx=*/0, /*compIdx=*/0, massFractionL[/*compIdx=*/0]*avgMolarMass/M1);
            fluidState.setMoleFraction(/*phaseIdx=*/0, /*compIdx=*/1, massFractionL[/*compIdx=*/1]*avgMolarMass/M2);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             /*phaseIdx=*/0,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }


    /*!
     * \brief Returns the relative permeability of a given phase
     *        within the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState().viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(/*phaseIdx=*/1) - fluidState_.pressure(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase
     */
    Scalar diffCoeff(int phaseIdx) const
    { return diffCoeff_[phaseIdx]; }

protected:
    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx));
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(typename FluidSystem::ParameterCache &paramCache,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int historyIdx)
    { }

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar relativePermeability_[numPhases]; //!< Relative permeability within the control volume
    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
