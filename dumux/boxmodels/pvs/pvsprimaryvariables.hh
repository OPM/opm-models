// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                                    *
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
 * \brief Represents the primary variables used in the primary variable switching compositional
 *         box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_PVS_PRIMARY_VARIABLES_HH
#define DUMUX_PVS_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "pvsindices.hh"
#include "pvsproperties.hh"

#include <iostream>

namespace Dumux
{
/*!
 * \ingroup PvsModel
 *
 * \brief Represents the primary variables used in the primary variable switching compositional
 *         box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class PvsPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;
    typedef PvsPrimaryVariables<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

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
    PvsPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetDefined(*this);
    }

    /*!
     * \brief Constructor with assignment from scalar
     */
    explicit PvsPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::CheckDefined(value);
        Valgrind::SetDefined(*this);

        phasePresence_ = 0;
    }

    /*!
     * \brief Copy constructor
     */
    PvsPrimaryVariables(const PvsPrimaryVariables &value)
        : ParentType(value)
    {
        Valgrind::SetDefined(*this);

        phasePresence_ = value.phasePresence_;
    }

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
    { return phasePresence_; }

    /*!
     * \brief Set which fluid phases are present in a given control
     *        volume.
     */
    void setPhasePresence(short value)
    { phasePresence_ = value; }

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
    { return lowestPresentPhaseIdx(); }
    
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
     * \brief Returns the phase with the lowest index that is present.
     */
    int lowestPresentPhaseIdx() const
    { return ffs(phasePresence_) - 1; }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const PrimaryVariables &value)
    { 
        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;

        return *this;
    }

    /*!
     * \brief Assignment operator
     */
    ThisType &operator=(const Scalar value)
    { 
        ParentType::operator=(value);

        phasePresence_ = 0;
        return *this;
    }

    /*!
     * \brief Returns an explcitly stored saturation for a given phase.
     * 
     * (or 0 if the saturation is not explicitly stored.)
     */
    Scalar explicitSaturationValue(int phaseIdx) const
    {
        if (!phaseIsPresent(phaseIdx) || phaseIdx == lowestPresentPhaseIdx()) 
            // non-present phases have saturation 0
            return 0.0;

        return (*this)[switch0Idx + phaseIdx - 1];
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

        // determine the phase presence.
        phasePresence_ = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // use a NCP condition to determine if the phase is
            // present or not
            Scalar a = 1;
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                a -= fluidState.moleFraction(phaseIdx, compIdx);
            }
            Scalar b = fluidState.saturation(phaseIdx);
            
            if (b > a)
                phasePresence_ |= (1 << phaseIdx);
        }
        
        // assert that some phase is present
        assert(phasePresence_ != 0);

        // set the primary variables which correspond to mole
        // fractions of the present phase which has the lowest index.
        int lowestPhaseIdx = lowestPresentPhaseIdx();
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (switchIdx >= lowestPhaseIdx)
                ++ phaseIdx;

            if (phaseIsPresent(phaseIdx))
            {
                (*this)[switch0Idx + switchIdx] = 
                    fluidState.saturation(phaseIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
            else {
                (*this)[switch0Idx + switchIdx] = 
                    fluidState.moleFraction(lowestPhaseIdx, compIdx);
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
        int lowestPhaseIdx = lowestPresentPhaseIdx();
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++ switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (phaseIdx >= lowestPhaseIdx)
                ++ phaseIdx; // skip the saturation of the present
                             // phase with the lowest index

            if  (phaseIsPresent(phaseIdx)) {
                os << ", S_" << FluidSystem::phaseName(phaseIdx)
                   << " = " << (*this)[switch0Idx + switchIdx];
            }
            else {
                os << ", x_" << FluidSystem::phaseName(lowestPhaseIdx)
                   << "^" << FluidSystem::componentName(compIdx)
                   << " = " << (*this)[switch0Idx + switchIdx];
            }
        };
        os << ")";
        os << ", phase presence: " << static_cast<int>(phasePresence_);
    }
    
protected:
    unsigned char phasePresence_;
};

} // end namepace

#endif
