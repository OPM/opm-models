// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
#ifndef DUMUX_3P3C_VOLUME_VARIABLES_HH
#define DUMUX_3P3C_VOLUME_VARIABLES_HH

#include "3p3cproperties.hh"

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/common/math.hh>

#include <dune/common/collectivecommunication.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePThreeCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class ThreePThreeCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    // constraint solvers
    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wCompIdx = Indices::wCompIdx,
        aCompIdx = Indices::aCompIdx,
        cCompIdx = Indices::cCompIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        pressure0Idx = Indices::pressure0Idx
    };

    // present phases
    enum {
        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly
    };


public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;


    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param primaryVars The primary variables
     * \param problem The problem
     * \param element The element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const ElementContext &elemCtx,
                int spaceIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx, spaceIdx, timeIdx);

        Implementation::updateTemperature_(fluidState_, elemCtx, spaceIdx, timeIdx);

        const auto &problem = elemCtx.problem();
        const auto &primaryVars = elemCtx.primaryVars(spaceIdx, timeIdx);
        int phasePresence = primaryVars.phasePresence();
        Valgrind::CheckDefined(primaryVars);
        
        // first, set the the saturations
        if (phasePresence == threePhases)
        {
            fluidState_.setSaturation(wPhaseIdx, primaryVars[switch1Idx]);
            fluidState_.setSaturation(nPhaseIdx, primaryVars[switch2Idx]);
            fluidState_.setSaturation(gPhaseIdx, 
                                      1.0
                                      - fluidState_.saturation(wPhaseIdx) 
                                      - fluidState_.saturation(nPhaseIdx));
        }
        else if (phasePresence == wPhaseOnly)
        {
            fluidState_.setSaturation(wPhaseIdx, 1.0);
            fluidState_.setSaturation(nPhaseIdx, 0.0);
            fluidState_.setSaturation(gPhaseIdx, 0.0);
        }
        else if (phasePresence == gnPhaseOnly)
        {
            fluidState_.setSaturation(wPhaseIdx, 0.0);
            fluidState_.setSaturation(nPhaseIdx, primaryVars[switch2Idx]);
            fluidState_.setSaturation(gPhaseIdx, 1.0 - primaryVars[switch2Idx]);
        }
        else if (phasePresence == wnPhaseOnly)
        {
            fluidState_.setSaturation(nPhaseIdx, primaryVars[switch2Idx]);
            fluidState_.setSaturation(wPhaseIdx, 1.0 - primaryVars[switch2Idx]);
            fluidState_.setSaturation(gPhaseIdx, 0.0);
        }
        else if (phasePresence == gPhaseOnly)
        {
            fluidState_.setSaturation(wPhaseIdx, 0.0);
            fluidState_.setSaturation(nPhaseIdx, 0.0);
            fluidState_.setSaturation(gPhaseIdx, 1.0);
        }
        else if (phasePresence == wgPhaseOnly)
        {
            fluidState_.setSaturation(wPhaseIdx, primaryVars[switch1Idx]);
            fluidState_.setSaturation(nPhaseIdx, 0.0);
            fluidState_.setSaturation(gPhaseIdx, 1.0 - primaryVars[switch1Idx]);
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(fluidState_.saturation(gPhaseIdx));
        Valgrind::CheckDefined(fluidState_.saturation(wPhaseIdx));
        Valgrind::CheckDefined(fluidState_.saturation(nPhaseIdx));

        // next, calculate the pressures
        Scalar p0 = primaryVars[pressure0Idx];
        
        // calculate capillary pressures
        const auto &matParams = problem.materialLawParams(elemCtx, spaceIdx, timeIdx);

        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, matParams, fluidState_);
        Valgrind::CheckDefined(pC);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            fluidState_.setPressure(phaseIdx, p0 + pC[phaseIdx] - pC[/*phaseIdx=*/0]);

        Scalar pg = fluidState_.pressure(gPhaseIdx);
        Scalar pw = fluidState_.pressure(wPhaseIdx);
        Scalar pn = fluidState_.pressure(nPhaseIdx);
        Valgrind::CheckDefined(pg);
        Valgrind::CheckDefined(pw);
        Valgrind::CheckDefined(pn);

        // now comes the tricky part: calculate phase composition
        typename FluidSystem::ParameterCache paramCache;
        if (phasePresence == threePhases) {
            // all phases are present, phase compositions are a
            // result of the the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            MiscibleMultiPhaseComposition::solve(fluidState_,
                                                 paramCache,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wPhaseOnly) {
            // only the water phase is present, water phase composition is
            // stored explicitly.

            // extract mole fractions in the water phase
            Scalar xwa = primaryVars[switch1Idx];
            Scalar xwc = primaryVars[switch2Idx];
            Scalar xww = 1 - xwa - xwc;

            // write water mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, aCompIdx, xwa);
            fluidState_.setMoleFraction(wPhaseIdx, cCompIdx, xwc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == gnPhaseOnly) {
            // only gas and NAPL phases are present
            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase

            // calculate and set the required fugacity coefficients 
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState_, paramCache, nPhaseIdx, cCompIdx);
            fluidState_.setFugacityCoefficient(nPhaseIdx, cCompIdx, phi);

            // we have all (partly hypothetical) phase pressures
            // and temperature and the mole fraction of water in
            // the gas phase
            Scalar partPressNAPL = fluidState_.fugacityCoefficient(nPhaseIdx, cCompIdx)*pn;

            Scalar xgw = primaryVars[switch1Idx];
            Scalar xgc = partPressNAPL/pg;
            Scalar xga = 1.-xgw-xgc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wnPhaseOnly) {
            // only water and NAPL phases are present

            // calculate and set the required fugacity coefficients 
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState_, paramCache, nPhaseIdx, cCompIdx);
            fluidState_.setFugacityCoefficient(nPhaseIdx, cCompIdx, phi);
            phi = FluidSystem::fugacityCoefficient(fluidState_, paramCache, wPhaseIdx, cCompIdx);
            fluidState_.setFugacityCoefficient(wPhaseIdx, cCompIdx, phi);

            Scalar pPartialC = fluidState_.fugacityCoefficient(nPhaseIdx,cCompIdx)*pn;
            Scalar henryC = fluidState_.fugacityCoefficient(wPhaseIdx,cCompIdx)*pw;

            Scalar xwa = primaryVars[switch1Idx];
            Scalar xwc = pPartialC/henryC;
            Scalar xww = 1.-xwa-xwc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            fluidState_.setMoleFraction(wPhaseIdx, aCompIdx, xwa);
            fluidState_.setMoleFraction(wPhaseIdx, cCompIdx, xwc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             wPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, gas phase composition is
            // stored explicitly here below.
            const Scalar xgw = primaryVars[switch1Idx];
            const Scalar xgc = primaryVars[switch2Idx];
            Scalar xga = 1 - xgw - xgc;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == wgPhaseOnly) {
            // only water and gas phases are present

            // calculate the required fugacity coefficients
            Scalar phi = FluidSystem::fugacityCoefficient(fluidState_, paramCache, wPhaseIdx, wCompIdx);
            fluidState_.setFugacityCoefficient(wPhaseIdx, wCompIdx, phi);

            Scalar xgc = primaryVars[switch2Idx];
            Scalar partPressH2O = fluidState_.fugacityCoefficient(wPhaseIdx, wCompIdx)*pw;

            Scalar xgw = partPressH2O/pg;
            Scalar xga = 1.-xgc-xgw;

            // write mole fractions in the fluid state
            fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            fluidState_.setMoleFraction(gPhaseIdx, aCompIdx, xga);
            fluidState_.setMoleFraction(gPhaseIdx, cCompIdx, xgc);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState_,
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else
            assert(false); // unhandled phase state

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       paramCache,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);
        }

        MaterialLaw::relativePermeabilities(relativePermeability_, matParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        /* ATTENTION: The conversion to effective diffusion parameters
         *            for the porous media happens in the flux
         *            variables!
         */

        // diffusivity coefficents
        diffusionCoefficient_[gPhaseIdx][wCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              gPhaseIdx,
                                              wCompIdx);
        diffusionCoefficient_[gPhaseIdx][cCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              gPhaseIdx,
                                              cCompIdx);
        diffusionCoefficient_[gPhaseIdx][aCompIdx] = 0.0; // dummy, should not be used !

        diffusionCoefficient_[wPhaseIdx][aCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              aCompIdx);
        diffusionCoefficient_[wPhaseIdx][cCompIdx] =
            FluidSystem::diffusionCoefficient(fluidState_,
                                              paramCache,
                                              wPhaseIdx,
                                              cCompIdx);
        diffusionCoefficient_[wPhaseIdx][wCompIdx] = 0.0; // dummy, should not be used !

        /* no diffusion in NAPL phase considered  at the moment */
        diffusionCoefficient_[nPhaseIdx][cCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][wCompIdx] = 0.0;
        diffusionCoefficient_[nPhaseIdx][aCompIdx] = 0.0;

        Valgrind::CheckDefined(diffusionCoefficient_);

        // porosity
        porosity_ = problem.porosity(elemCtx, spaceIdx, timeIdx);
        Valgrind::CheckDefined(porosity_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(paramCache, elemCtx, spaceIdx, timeIdx);
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
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the molecular diffusion coefficent of a component in a phase;
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    { return diffusionCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    {}                                    
    
    /*!
     * \brief Set the enthalpy rate per second of a rate vector, .
     */
    static void setEnthalpyRate(RateVector &rateVec, Scalar rate)
    { }

    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class FluidState>
    static void setEnthalpyRate(RateVector &v,
                                const FluidState &fluidState, 
                                int phaseIdx, 
                                Scalar volume)
    { }

protected:
    static void updateTemperature_(FluidState &fs, const ElementContext &elemCtx, int spaceIdx, int timeIdx)
    { fs.setTemperature(elemCtx.problem().temperature(elemCtx, spaceIdx, timeIdx)); }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    template <class ParameterCache>
    void updateEnergy_(const ParameterCache &paramCache, 
                       const ElementContext &elemCtx,
                       int spaceIdx,
                       int timeIdx)
    {}

    Scalar porosity_;
    Scalar relativePermeability_[numPhases];
    Scalar diffusionCoefficient_[numPhases][numComponents];

    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
