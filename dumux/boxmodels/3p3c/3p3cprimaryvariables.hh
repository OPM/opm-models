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
 * \brief Represents the primary variables used in the 3-phase,
 *        3-component box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_3P3C_PRIMARY_VARIABLES_HH
#define DUMUX_3P3C_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "3p3cindices.hh"
#include "3p3cproperties.hh"

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
class ThreePThreeCPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;
    typedef ThreePThreeCPrimaryVariables<TypeTag> ThisType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCIndices) Indices;

    // primary variable indices
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { switch1Idx = Indices::switch1Idx };
    enum { switch2Idx = Indices::switch2Idx };

    enum { wPhaseIdx = Indices::wPhaseIdx };
    enum { nPhaseIdx = Indices::nPhaseIdx };
    enum { gPhaseIdx = Indices::gPhaseIdx };

    enum { wCompIdx = Indices::wCompIdx };
    enum { aCompIdx = Indices::aCompIdx };
    enum { cCompIdx = Indices::cCompIdx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) EnergyModule;

    typedef Dumux::NcpFlash<Scalar, FluidSystem> NcpFlash;

public:
    /*!
     * \brief Default constructor
     */
    ThreePThreeCPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetDefined(*this);
        Valgrind::SetUndefined(*static_cast<ParentType*>(this));

        wasSwitched_ = false;
        phasePresence_ = 0;
    };

    /*!
     * \brief Constructor with assignment from scalar
     */
    explicit ThreePThreeCPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::SetDefined(*this);
        Valgrind::CheckDefined(value);

        phasePresence_ = 0;
        wasSwitched_ = false;
    };

    /*!
     * \brief Copy constructor
     */
    ThreePThreeCPrimaryVariables(const ThreePThreeCPrimaryVariables &value)
        : ParentType(value)
    {
        Valgrind::SetDefined(*this);

        phasePresence_ = value.phasePresence_;
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
    { phasePresence_ = value; };

    /*!
     * \brief Set whether a given indivividual phase should be present
     *        or not.
     */
    void setPhasePresent(int phaseIdx, bool yesno = true)
    { 
        if (yesno) phasePresence_ |= (1 << phaseIdx);
        else phasePresence_ &= ~(1 << phaseIdx);
    }

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
    { return wasSwitched_; }

    /*!
     * \brief Set whether the primary variables where switched in the last iteration
     */
    void setSwitched(bool value)
    { wasSwitched_ = value; }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const Implementation &value)
    { 
        Valgrind::SetDefined(*this);

        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;
        wasSwitched_ = value.wasSwitched_;

        return *this;
    }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const Scalar value)
    {
        Valgrind::SetDefined(*this);
        ParentType::operator=(value);

        phasePresence_ = 0;
        wasSwitched_ = false;
        return *this;
    }

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

        determinePhasePresence_(fluidState);

        // set the pressure
        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);

        if (phasePresence_ == Indices::threePhases) {
            (*this)[switch1Idx] = fluidState.saturation(wPhaseIdx);
            (*this)[switch2Idx] = fluidState.saturation(nPhaseIdx);
        }
        else if (phasePresence_ == Indices::wPhaseOnly) {
            (*this)[switch1Idx] = fluidState.moleFraction(wPhaseIdx, aCompIdx);
            (*this)[switch2Idx] = fluidState.moleFraction(wPhaseIdx, cCompIdx);
        }
        else if (phasePresence_ == Indices::gnPhaseOnly) {
            (*this)[switch1Idx] = fluidState.moleFraction(gPhaseIdx, wCompIdx);
            (*this)[switch2Idx] = fluidState.saturation(nPhaseIdx);
        }
        else if (phasePresence_ == Indices::wnPhaseOnly) {
            (*this)[switch1Idx] = fluidState.moleFraction(wPhaseIdx, aCompIdx);
            (*this)[switch2Idx] = fluidState.saturation(nPhaseIdx);
        }
        else if (phasePresence_ == Indices::gPhaseOnly) {
            (*this)[switch1Idx] = fluidState.moleFraction(gPhaseIdx, wCompIdx);
            (*this)[switch2Idx] = fluidState.moleFraction(gPhaseIdx, cCompIdx);
        }
        else if (phasePresence_ == Indices::wgPhaseOnly) {
            (*this)[switch1Idx] = fluidState.saturation(wPhaseIdx);
            (*this)[switch2Idx] = fluidState.moleFraction(gPhaseIdx, cCompIdx);
        }
        else
            assert(false);
    }
           
protected:
    template <class FluidState> 
    void determinePhasePresence_(const FluidState &fluidState)
    {
        // determine the phase presence.
        phasePresence_ = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // NCP-like condition
            Scalar a = 1;
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                a -= fluidState.moleFraction(phaseIdx, compIdx);
            Scalar b = fluidState.saturation(phaseIdx);
            
            if (b > a) {
                phasePresence_ |= (1 << phaseIdx);
            }
        }
    }

    unsigned char phasePresence_;
    bool wasSwitched_;
};

} // end namepace

#endif
