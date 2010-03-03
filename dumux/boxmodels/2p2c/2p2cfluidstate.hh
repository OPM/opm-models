/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Calcultes the phase state from the primary variables in the
 *        2pNc model.
 */
#ifndef DUMUX_2P2C_PHASE_STATE_HH
#define DUMUX_2P2C_PHASE_STATE_HH

#include "2p2cproperties.hh"

#include <dumux/new_material/fluidstate.hh>

namespace Dune
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        2pNc model.
 */
template <class TypeTag>
class TwoPTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                             TwoPTwoCFluidState<TypeTag> >
{
    typedef TwoPTwoCFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params                       MaterialLawParams;

    enum {
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,
        
        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx,
    };

    // present phases
    enum {
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum { 
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),
        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl
    };

public:   
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { numSolvents = 1 };

    /*!
     * \brief Update the phase state from the primary variables.
     */
    void update(const PrimaryVarVector &primaryVars,
                const MaterialLawParams &pcParams,
                Scalar temperature,
                int phasePresence)
    {
        Valgrind::CheckDefined(primaryVars);

        temperature_ = temperature;
        
        // extract non-wetting phase pressure
        if (phasePresence == gPhaseOnly)
            Sg_ = 1.0;
        else if (phasePresence == lPhaseOnly) {
            Sg_ = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == plSg)
                Sg_ = primaryVars[switchIdx];
            else if (formulation == pgSl)
                Sg_ = 1.0 - primaryVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(Sg_);

        // calculate capillary pressure
        Scalar pC = MaterialLaw::pC(pcParams, 1 - Sg_);

        // extract the pressures
        if (formulation == plSg) {
            phasePressure_[lPhaseIdx] = primaryVars[pressureIdx];
            phasePressure_[gPhaseIdx] = phasePressure_[lPhaseIdx] + pC;
        }
        else if (formulation == pgSl) {
            phasePressure_[gPhaseIdx] = primaryVars[pressureIdx];
            phasePressure_[lPhaseIdx] = phasePressure_[gPhaseIdx] - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        Valgrind::CheckDefined(phasePressure_);
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(Sg_);

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases) {
            // both phases are present, phase composition results from
            // the gas <-> liquid equilibrium.
            Scalar pn = phasePressure_[gPhaseIdx];
            Scalar beta0 = FluidSystem::degasPressure(lCompIdx, temperature_, pn);
            Scalar beta1 = FluidSystem::degasPressure(gCompIdx, temperature_, pn);

            // calculate _mole_ (!) fractions in wetting phase,
            // mis-use the concentration_ array temporary storage
            concentration_[lPhaseIdx][lCompIdx] = (pn - beta1)/(beta0 - beta1);
            concentration_[lPhaseIdx][gCompIdx] = 1 - concentration_[lPhaseIdx][lCompIdx];

            // calculate the partial pressures
            partialPressure_[lCompIdx] = beta0*concentration_[lPhaseIdx][lCompIdx];
            partialPressure_[gCompIdx] = beta1*concentration_[lPhaseIdx][gCompIdx];

            // calculate the average molar mass of the wetting
            // phase. keep in mind that concentration_ still stores
            // the mole fractions at this point
            avgMolarMass_[lPhaseIdx] = 
                concentration_[lPhaseIdx][lCompIdx] * FluidSystem::molarMass(lCompIdx) +
                concentration_[lPhaseIdx][gCompIdx] * FluidSystem::molarMass(gCompIdx);
            
            // temporarily set totalConcentration_ to 1, so that the
            // moleFrac() and massFrac() methods return correct values
            // which is a requisite of calculating the phase density
            totalConcentration_[lPhaseIdx] = 1.0;
            // ask the fluid system to calculate the density of the
            // wetting phase
            Scalar lPhaseDensity = 
                FluidSystem::phaseDensity(lPhaseIdx,
                                          temperature_,
                                          phasePressure_[lPhaseIdx],
                                          *this);
            // calculate the real total concentration from the phase density
            totalConcentration_[lPhaseIdx] = lPhaseDensity / avgMolarMass_[lPhaseIdx];

            // finally undo the mis-use of concentration_ by
            // converting the mole fractions into concentrations
            concentration_[lPhaseIdx][lCompIdx] = 
                concentration_[lPhaseIdx][lCompIdx] * totalConcentration_[lPhaseIdx];
            concentration_[lPhaseIdx][gCompIdx] = 
                concentration_[lPhaseIdx][gCompIdx] * totalConcentration_[lPhaseIdx];
            
            // now we deal with the gas phase. we already have the
            // partial pressures, we now let the fluid system convert
            // these to concentrations. Generally, this involves to
            // resolve the gas' equation of state to the specific
            // volume.
            concentration_[gPhaseIdx][lCompIdx] = 
                FluidSystem::componentDensity(gPhaseIdx,
                                              lCompIdx,
                                              temperature_, 
                                              partialPressure_[lCompIdx])
                /
                FluidSystem::molarMass(lCompIdx);
            concentration_[gPhaseIdx][gCompIdx] = 
                FluidSystem::componentDensity(gPhaseIdx,
                                              gCompIdx,
                                              temperature_, 
                                              partialPressure_[gCompIdx])
                /
                FluidSystem::molarMass(gCompIdx);

            // Using the component concentrations, we can calculate
            // the total concentration in the gas phase
            totalConcentration_[gPhaseIdx] = 
                concentration_[gPhaseIdx][lCompIdx] +
                concentration_[gPhaseIdx][gCompIdx];
            // the total concentration now allows us to calculate the
            // average molar mass in the gas phase
            avgMolarMass_[gPhaseIdx] = 
                FluidSystem::molarMass(lCompIdx)*moleFrac(gPhaseIdx, lCompIdx) +
                FluidSystem::molarMass(gCompIdx)*moleFrac(gPhaseIdx, gCompIdx);
            
            // make sure everything was initialized from a defined
            // memory
            Valgrind::CheckDefined(partialPressure_);
            Valgrind::CheckDefined(concentration_);
            Valgrind::CheckDefined(totalConcentration_);
            Valgrind::CheckDefined(avgMolarMass_);
        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, gas phase composition is
            // stored explicitly.

            // extract _mass_ (!) fraction in the gas phase, misuse
            // the concentration_ array for temporary storage
            concentration_[gPhaseIdx][lCompIdx] = primaryVars[switchIdx];
            concentration_[gPhaseIdx][gCompIdx] = 1 - concentration_[gPhaseIdx][lCompIdx];

            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = concentration_[gPhaseIdx][gCompIdx]; // mass fraction of solvent in gas
            avgMolarMass_[gPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions
            concentration_[gPhaseIdx][lCompIdx] *= avgMolarMass_[gPhaseIdx]/M1;
            concentration_[gPhaseIdx][gCompIdx] *= avgMolarMass_[gPhaseIdx]/M2;
            
            // temporary set the total concentration 1.0 to make the
            // massFrac() and moleFrac() methods return correct values
            totalConcentration_[gPhaseIdx] = 1.0;

            // ask the fluid system to calculate the partial pressures
            // of the components. The FluidSystem calls
            // setPartialPressure().
            FluidSystem::computePartialPressures(temperature_,
                                                 phasePressure_[gPhaseIdx],
                                                 *this);

            // Ask the fluid system for the total gas density
            Scalar gPhaseDensity = 
                FluidSystem::phaseDensity(gPhaseIdx,
                                          temperature_,
                                          phasePressure_[gPhaseIdx],
                                          *this);
            
            // convert mole fractions to concentrations
            concentration_[gPhaseIdx][lCompIdx] *= gPhaseDensity/avgMolarMass_[gPhaseIdx];
            concentration_[gPhaseIdx][gCompIdx] *= gPhaseDensity/avgMolarMass_[gPhaseIdx];

            // calculate the real total concentration
            totalConcentration_[gPhaseIdx]  = concentration_[gPhaseIdx][lCompIdx];
            totalConcentration_[gPhaseIdx] += concentration_[gPhaseIdx][gCompIdx];

            ////////////////////////////////////////////////////////
            // calculate values for the non-existing liquid phase. we
            // need this for example for meaningful pressure potential
            // gradients in the case where liquid is present in one cell
            // and not pressent in the neighboring one.
            ////////////////////////////////////////////////////////
            Scalar beta0 = FluidSystem::degasPressure(lCompIdx,
                                                      temperature_,
                                                      phasePressure_[gPhaseIdx]);
            Scalar beta1 = FluidSystem::degasPressure(gCompIdx,
                                                      temperature_,
                                                      phasePressure_[gPhaseIdx]);
            concentration_[lPhaseIdx][lCompIdx] =
                partialPressure_[lCompIdx] / beta0;
            concentration_[lPhaseIdx][gCompIdx] =
                partialPressure_[gCompIdx] / beta1;
            // average molar mass in liquid "phase". keep in mind that
            // concentration_[lPhase] stores the mole fractions at
            // this point.
            avgMolarMass_[lPhaseIdx]  = concentration_[lPhaseIdx][lCompIdx]*FluidSystem::molarMass(lCompIdx);
            avgMolarMass_[lPhaseIdx] += concentration_[lPhaseIdx][gCompIdx]*FluidSystem::molarMass(gCompIdx);
            avgMolarMass_[lPhaseIdx] /=
                concentration_[lPhaseIdx][lCompIdx] + 
                concentration_[lPhaseIdx][gCompIdx];

            // set total concentration of liquid phase to 1, so that
            // massFrac() and moleFrac() return meaningful values, and
            // ask the fluid system to return the total density of the
            // liquid phase
            totalConcentration_[lPhaseIdx] = 1.0;
            Scalar lPhaseDensity = FluidSystem::phaseDensity(lPhaseIdx,
                                                             temperature_,
                                                             phasePressure_[lPhaseIdx],
                                                             *this);
            // convert totalConcentration_ to the real value
            totalConcentration_[lPhaseIdx] = lPhaseDensity/avgMolarMass_[lPhaseIdx];

            // calculate the "real" concentrations in the non-existing
            // liquid phase
            concentration_[lPhaseIdx][lCompIdx] *= totalConcentration_[lPhaseIdx];
            concentration_[lPhaseIdx][gCompIdx] *= totalConcentration_[lPhaseIdx];

            // make sure everything was initialized from a defined
            // memory
            Valgrind::CheckDefined(partialPressure_);
            Valgrind::CheckDefined(concentration_);
            Valgrind::CheckDefined(totalConcentration_);
            Valgrind::CheckDefined(avgMolarMass_);
        }
        else if (phasePresence == lPhaseOnly) {
            // only the liquid phase is present, i.e. liquid phase
            // composition is stored explicitly.

            // extract _mass_ (!) fraction in the liquid phase, mis-use
            // the concentration_ array for temporary storage
            concentration_[lPhaseIdx][gCompIdx] = primaryVars[switchIdx];
            concentration_[lPhaseIdx][lCompIdx] = 1 - concentration_[lPhaseIdx][gCompIdx];

            // calculate average molar mass of the liquid phase
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = concentration_[lPhaseIdx][gCompIdx]; // mass fraction of gas in liquid
            avgMolarMass_[lPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass fractions to mole fractions
            concentration_[lPhaseIdx][lCompIdx] *= avgMolarMass_[lPhaseIdx]/M1;
            concentration_[lPhaseIdx][gCompIdx] *= avgMolarMass_[lPhaseIdx]/M2;

            // set the total concentration to 1 to make the moleFrac()
            // and massFrac() methods return meaningful values
            totalConcentration_[lPhaseIdx] = 1.0;

            // ask the fluid system for the liquid density
            Scalar lPhaseDensity = 
                FluidSystem::phaseDensity(lPhaseIdx,
                                          temperature_,
                                          phasePressure_[lPhaseIdx],
                                          *this);
            
            // convert mole fractions to concentrations
            concentration_[lPhaseIdx][lCompIdx] *= lPhaseDensity/avgMolarMass_[lPhaseIdx];
            concentration_[lPhaseIdx][gCompIdx] *= lPhaseDensity/avgMolarMass_[lPhaseIdx];
            
            // calculate the real total concentration
            totalConcentration_[lPhaseIdx]  = concentration_[lPhaseIdx][lCompIdx];
            totalConcentration_[lPhaseIdx] += concentration_[lPhaseIdx][gCompIdx];
            
            ////////////////////////////////////////////////////////
            // calculate values for the non-existing gas phase. we
            // need this for example for meaningful pressure potential
            // gradients in the case where gas is present in one cell
            // and not pressent in the neighboring one.
            ////////////////////////////////////////////////////////
            
            // calculate partial pressures.
            Scalar beta0 = FluidSystem::degasPressure(lCompIdx, temperature_, phasePressure_[gPhaseIdx]);
            Scalar beta1 = FluidSystem::degasPressure(gCompIdx, temperature_, phasePressure_[gPhaseIdx]);
            partialPressure_[lCompIdx] = beta0*moleFrac(lPhaseIdx, lCompIdx);
            partialPressure_[gCompIdx] = beta1*moleFrac(lPhaseIdx, gCompIdx);
            
            // calculate concentrations in the (non-existing) gas
            // phase by inverting each partial pressure
            // individually. this assumes that dalton's law holds.
            concentration_[gPhaseIdx][lCompIdx] = 
                FluidSystem::componentDensity(gPhaseIdx,
                                              lCompIdx,
                                              temperature_, 
                                              partialPressure_[lCompIdx])
                /
                FluidSystem::molarMass(lCompIdx);
            concentration_[gPhaseIdx][gCompIdx] =
                FluidSystem::componentDensity(gPhaseIdx,
                                              gCompIdx,
                                              temperature_, 
                                              partialPressure_[gCompIdx])
                /
                FluidSystem::molarMass(gCompIdx);

            // calculate the real total concentration in the gas phase
            totalConcentration_[gPhaseIdx]  = concentration_[gPhaseIdx][lCompIdx];
            totalConcentration_[gPhaseIdx] += concentration_[gPhaseIdx][gCompIdx];
            
            // average molar mass in gas phase.
            avgMolarMass_[gPhaseIdx]  = moleFrac(gPhaseIdx, lCompIdx)*FluidSystem::molarMass(lCompIdx);
            avgMolarMass_[gPhaseIdx] += moleFrac(gPhaseIdx, gCompIdx)*FluidSystem::molarMass(gCompIdx);

            // make sure everything was initialized from a defined
            // memory
            Valgrind::CheckDefined(partialPressure_);
            Valgrind::CheckDefined(concentration_);
            Valgrind::CheckDefined(totalConcentration_);
            Valgrind::CheckDefined(avgMolarMass_);
        }
    }

public:
    /*!
     * \brief Returns the saturation of a phase.
     */
    Scalar saturation(int phaseIdx) const
    { 
        if (phaseIdx == lPhaseIdx)
            return Scalar(1.0) - Sg_;
        else
            return Sg_;
    };
    
    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        return
            concentration_[phaseIdx][compIdx]/
            totalConcentration_[phaseIdx];
    }
    
    /*!
     * \brief Returns the total concentration of a phase [mol / m^3].
     *
     * This is equivalent to the sum of all component concentrations.
     */
    Scalar totalConcentration(int phaseIdx) const
    { return totalConcentration_[phaseIdx]; };

    /*!
     * \brief Returns the concentration of a component in a phase [mol / m^3].
     */
    Scalar concentration(int phaseIdx, int compIdx) const
    { return concentration_[phaseIdx][compIdx]; };
    
    /*!
     * \brief Returns the mass fraction of a component in a phase.
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    { 
        return
            moleFrac(phaseIdx, compIdx) *
            FluidSystem::molarMass(compIdx)/avgMolarMass_[phaseIdx];
    }
    
    /*!
     * \brief Returns the density of a phase [kg / m^3].
     */
    Scalar density(int phaseIdx) const
    { return totalConcentration_[phaseIdx]*avgMolarMass_[phaseIdx]; }
    
    /*!
     * \brief Returns mean molar mass of a phase [kg / mol].
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return avgMolarMass_[phaseIdx]; };
    
    /*!
     * \brief Returns the partial pressure of a component in the gas phase [Pa].
     */
    Scalar partialPressure(int compIdx) const
    { return partialPressure_[compIdx]; }

    /*!
     * \brief Set the partial pressure of a component in the gas phase
     *        [Pa].
     *
     * This method is required in order to use
     * FluidSystem::computePartialPressures().
     */
    void setPartialPressure(int compIdx, Scalar value)
    { partialPressure_[compIdx] = value; }
    
    /*!
     * \brief Returns the pressure of a fluid phase [Pa].
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure [Pa]
     */
    Scalar capillaryPressure() const
    { return phasePressure_[gPhaseIdx] - phasePressure_[lPhaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids [K].
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

    /*!
     * \brief Return henry coefficent or vapor pressure, depending on wether
     *        the component is the solvent or a solute.
     */
    Scalar beta(int compIdx) const
    { return partialPressure_[compIdx]/moleFrac(lPhaseIdx, compIdx); }

public:
    Scalar partialPressure_[numComponents];
    Scalar concentration_[numPhases][numComponents];
    Scalar totalConcentration_[numPhases];
    Scalar avgMolarMass_[numPhases];
    Scalar phasePressure_[numComponents];
    Scalar temperature_;
    Scalar Sg_;
};

} // end namepace

#endif
