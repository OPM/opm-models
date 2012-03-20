// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Represents the primary variables used in the 2-phase,
 *        2-component box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_2P2C_PRIMARY_VARIABLES_HH
#define DUMUX_2P2C_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "2p2cindices.hh"
#include "2p2cproperties.hh"

#include <iostream>

namespace Dumux
{
/*!
 * \ingroup 2P2CModel
 *
 * \brief Represents the primary variables used in the 2-phase,
 *        2-component box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class TwoPTwoCPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;
    typedef TwoPTwoCPrimaryVariables<TypeTag> ThisType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { switch0Idx = Indices::switch0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) EnergyModule;

    typedef Dumux::NcpFlash<Scalar, FluidSystem> NcpFlash;

public:
    /*!
     * \brief Default constructor
     */
    TwoPTwoCPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetDefined(*this);

        wasSwitched_ = false;
        phasePresence_ = 0;
        lowestPresentPhaseIdx_ = -1;
    };

    /*!
     * \brief Constructor with assignment from scalar
     */
    explicit TwoPTwoCPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::CheckDefined(value);
        Valgrind::SetDefined(*this);

        phasePresence_ = 0;
        lowestPresentPhaseIdx_ = -1;
        wasSwitched_ = false;
    };

    /*!
     * \brief Copy constructor
     */
    TwoPTwoCPrimaryVariables(const TwoPTwoCPrimaryVariables &value)
        : ParentType(value)
    {
        Valgrind::SetDefined(*this);

        phasePresence_ = value.phasePresence_;
        lowestPresentPhaseIdx_ = value.lowestPresentPhaseIdx_;
        wasSwitched_ = false;
    };

    /*!
     * \brief Set the primary variables from an arbitrary fluid state
     *        in a mass conservative way.
     *
     * If an energy equation is included, the fluid temperatures are
     * the same as the one given in the fluid state, *not* the
     * enthalpy.
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     * \param matParams The capillary pressure law parameters
     * \param isInEquilibrium If true, the fluid state expresses
     *                        thermodynamic equilibrium assuming the
     *                        relations expressed by the fluid
     *                        system. This implies that in addition to
     *                        the quantities mentioned above, the
     *                        fugacities are also defined.
     */
    template <class FluidState>
    void assignMassConservative(const FluidState &fluidState,
                                const MaterialLawParams &matParams,
                                bool isInEquilibrium = false)
    {
#ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            assert(fluidState.temperature(0) == fluidState.temperature(phaseIdx));
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium
        typename FluidSystem::ParameterCache paramCache;
        Dumux::CompositionalFluidState<Scalar, FluidSystem> fsFlash;
             
        // use the externally given fluid state as initial value for
        // the flash calculation
        fsFlash.assign(fluidState);
        
        // calculate the phase densities
        paramCache.updateAll(fsFlash);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fsFlash.setDensity(phaseIdx, 
                               FluidSystem::density(fsFlash,
                                                    paramCache,
                                                    phaseIdx));
        }

        // calculate the "global molarities" 
        ComponentVector globalMolarities(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    fsFlash.saturation(phaseIdx)*fsFlash.molarity(phaseIdx, compIdx);
            }
        }

        // run the flash calculation
        //NcpFlash::guessInitial(fsFlash, paramCache, globalMolarities);       
        NcpFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams, globalMolarities);
        
        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \brief Return the fluid phases which are present in a given
     *        control volume.
     */
    short phasePresence() const
    { return phasePresence_; };

    /*!
     * \brief Set which fluid phases are present in a given control
     *        volume.
     */
    void setPhasePresence(short value)
    { 
        phasePresence_ = value;
        
        lowestPresentPhaseIdx_ = -1;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIsPresent(phaseIdx)) {
                lowestPresentPhaseIdx_ = phaseIdx;
                break;
            }
        }
    };

    /*!
     * \brief Set whether a given indivividual phase should be present
     *        or not.
     */
    void setPhasePresent(int phaseIdx, bool yesno = true)
    { 
        if (yesno) setPhasePresence(phasePresence_ | (1 << phaseIdx));
        else setPhasePresence(phasePresence_ & ~(1 << phaseIdx));
    }

    /*!
     * \brief Returns the phase with is determined by the closure
     *        condition of saturation.
     */
    unsigned char implicitSaturationIdx() const 
    { return lowestPresentPhaseIdx_; }
    
    /*!
     * \brief Returns true iff a phase is present for a given phase presence.
     */
    static bool phaseIsPresent(int phaseIdx, short phasePresence)
    { return phasePresence & (1 << phaseIdx); }

    /*!
     * \brief Returns true iff a phase is present for the current phase presence.
     */
    bool phaseIsPresent(int phaseIdx) const
    { return phasePresence_ & (1 << phaseIdx); }

    /*!
     * \brief Returns whether the primary variables where switched in the last iteration
     */
    bool wasSwitched() const
    { 
        return wasSwitched_;
    }

    /*!
     * \brief Set whether the primary variables where switched in the last iteration
     */
    void setSwitched(bool value)
    { 
        Valgrind::CheckDefined(value);
        wasSwitched_ = value;
    }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const ThisType &value)
    { 
        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;
        lowestPresentPhaseIdx_ = value.lowestPresentPhaseIdx_;
        wasSwitched_ = value.wasSwitched_;

        return *this;
    }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const Scalar value)
    { 
        ParentType::operator=(value);

        phasePresence_ = 0;
        lowestPresentPhaseIdx_ = -1;
        wasSwitched_ = false;
        return *this;
    }
  
    /*!
     * \brief Returns an explcitly stored saturation for a given phase.
     * 
     * (or 0 if the saturation is not explicitly stored.)
     */
    Scalar explicitSaturationValue(int phaseIdx) const
    {
        if (phaseIdx == lowestPresentPhaseIdx_)
            // the saturation of the present phase with the lowest is
            // not stored explicitly
            return 0.0;
        else if (!phaseIsPresent(phaseIdx)) 
            // non-present phases have saturation 0
            return 0.0;

        return (*this)[switch0Idx + phaseIdx - 1];
    };

    /*!
     * \brief Assign the primary variables "naively" from a fluid state.
     *
     * \attention Some mass might get lost/added if the fluid state is
     *            not in thermodynamic equilibrium!
     */
    template <class FluidState> 
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // set the pressure of the first phase
        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        Valgrind::CheckDefined((*this)[pressure0Idx]);

        // determine the phase presence.
        phasePresence_ = 0;
        lowestPresentPhaseIdx_ = -1;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // determine wether the phase is be present or not. There
            // are three cases where a phase is present:
            // 
            // a) the current saturation is larger than 0
            // b) the phase appeared in the current iteration and the
            //    current saturation is smaller than the switch tolerance           
            // c) the phase was not present in the previous iteration
            //    but the sum of the mole fractions in the phase is
            //    above the switch tolerance
            if (fluidState.saturation(phaseIdx) > 0) {
                // case a)
                phasePresence_ |= (1 << phaseIdx);
            }
            else {
                // case c)
                Scalar sumMoleFrac = 0;
                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);

                if (sumMoleFrac > 1.0) {
                    phasePresence_ |= (1 << phaseIdx);
                }
            }
            
            if (lowestPresentPhaseIdx_ < 0 && phaseIsPresent(phaseIdx))
                lowestPresentPhaseIdx_ = phaseIdx;
        }
        
        // assert that some phase is present
        assert(phasePresence_ != 0);

        // set the primary variables which correspond to mole
        // fractions of the present phase which has the lowest index.
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (switchIdx >= lowestPresentPhaseIdx_)
                ++ phaseIdx;

            if (phaseIsPresent(phaseIdx))
            {
                (*this)[switch0Idx + switchIdx] = 
                    fluidState.saturation(phaseIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
            else {
                (*this)[switch0Idx + switchIdx] = 
                    fluidState.moleFraction(lowestPresentPhaseIdx_, compIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
        }
    }

    /*!
     * \brief Assign the primary variables "naively" from a fluid state.
     *
     * \attention Some mass might get lost/added if the fluid state is
     *            not in thermodynamic equilibrium!
     */
    template <class FluidState> 
    void assignNaive(const FluidState &fluidState, short oldPhasePresence)
    {
        Scalar switchTol = 0.0;
        if (wasSwitched()) {
            switchTol = 0.02;
        }

        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        // set the pressure of the first phase
        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        Valgrind::CheckDefined((*this)[pressure0Idx]);

        // determine the phase presence.
        phasePresence_ = 0;
        lowestPresentPhaseIdx_ = -1;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // determine wether the phase is be present or not. There
            // are three cases where a phase is present:
            // 
            // a) the current saturation is larger than 0
            // b) the phase appeared in the current iteration and the
            //    current saturation is smaller than the switch tolerance           
            // c) the phase was not present in the previous iteration
            //    but the sum of the mole fractions in the phase is
            //    above the switch tolerance
            if (fluidState.saturation(phaseIdx) > 0) {
                // case a)
                phasePresence_ |= (1 << phaseIdx);
            }
            else if (phaseIsPresent(phaseIdx, oldPhasePresence) &&
                     fluidState.saturation(phaseIdx) > -switchTol) {
                // case b)
                phasePresence_ |= (1 << phaseIdx);
            }
            else if (!phaseIsPresent(phaseIdx, oldPhasePresence)) {
                // case c)
                Scalar sumMoleFrac = 0;
                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);

                if (sumMoleFrac > 1.0 + switchTol) {
                    phasePresence_ |= (1 << phaseIdx);
                }
            }
            
            if (lowestPresentPhaseIdx_ < 0 && phaseIsPresent(phaseIdx))
                lowestPresentPhaseIdx_ = phaseIdx;
        }
        
        // assert that some phase is present
        assert(phasePresence_ != 0);
      
        // set the primary variables which correspond to mole
        // fractions of the present phase which has the lowest index.
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (switchIdx >= lowestPresentPhaseIdx_)
                ++ phaseIdx;

            if (phaseIsPresent(phaseIdx))
            {
                if (!phaseIsPresent(phaseIdx, oldPhasePresence))
                    (*this)[switch0Idx + switchIdx] = 0.001;
                else if (!phaseIsPresent(lowestPresentPhaseIdx_, oldPhasePresence))
                    (*this)[switch0Idx + switchIdx] = 
                        fluidState.saturation(phaseIdx) - 0.001;
                else
                    (*this)[switch0Idx + switchIdx] = 
                        fluidState.saturation(phaseIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
            else {
                (*this)[switch0Idx + switchIdx] = 
                    fluidState.moleFraction(lowestPresentPhaseIdx_, compIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
        }
    }

    /*!
     * \brief Prints the names of the primary variables and their values.
     */
    void print(std::ostream &os = std::cout) const
    {
        os << "(p_" << FluidSystem::phaseName(0)
           << " = " << this->operator[](pressure0Idx);
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++ switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (phaseIdx >= lowestPresentPhaseIdx_)
                ++ phaseIdx; // skip the saturation of the present
                             // phase with the lowest index

            if  (phaseIsPresent(phaseIdx)) {
                os << ", S_" << FluidSystem::phaseName(phaseIdx)
                   << " = " << (*this)[switch0Idx + switchIdx];
            }
            else {
                os << ", x_" << FluidSystem::phaseName(lowestPresentPhaseIdx_)
                   << "^" << FluidSystem::componentName(compIdx)
                   << " = " << (*this)[switch0Idx + switchIdx];
            }
        };
        os << ")";
        os << ", phase presence: " << static_cast<int>(phasePresence_);
    };
    
protected:
    unsigned char phasePresence_;
    char lowestPresentPhaseIdx_;
    bool wasSwitched_;
};

} // end namepace

#endif
