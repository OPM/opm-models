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

namespace Dumux
{
/*!
 * \ingroup 2P2CModel
 *
 * \brief Represents the primary variables used in the M-phase,
 *        N-component box model.
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
    enum { pressureIdx = Indices::pressureIdx };
    enum { switchIdx = Indices::switchIdx };

    // phase presence
    enum {
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    enum {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        
        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl
    };

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

        phasePresence_ = -123;
        Valgrind::SetUndefined(*static_cast<ParentType*>(this));
        Valgrind::SetUndefined(phasePresence_);
        wasSwitched_ = false;
    };

    /*!
     * \brief Constructor with assignment from scalar
     */
    TwoPTwoCPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::CheckDefined(value);

        phasePresence_ = -123;
        Valgrind::SetDefined(*this);
        Valgrind::SetUndefined(phasePresence_);
        wasSwitched_ = false;
    };

    /*!
     * \brief Copy constructor
     */
    TwoPTwoCPrimaryVariables(const TwoPTwoCPrimaryVariables &value)
        : ParentType(value)
    {
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
    int phasePresence() const
    { return phasePresence_; };

    /*!
     * \brief Set which fluid phases are present in a given control
     *        volume.
     */
    void setPhasePresence(int value)
    { phasePresence_ = value; };

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
    ThisType &operator=(const ThisType &value)
    { 
        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;
        wasSwitched_ = false;
        return *this;
    }
    
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        if (formulation == plSg)
            (*this)[pressureIdx] = fluidState.pressure(/*phaseIdx=*/0);
        else if (formulation == pgSl)
            (*this)[pressureIdx] = fluidState.pressure(/*phaseIdx=*/1);
        else {
            // invalid formulation
            assert(false);
        }
        
        Scalar Sl = fluidState.saturation(/*phaseIdx=*/0);
        Scalar Sg = fluidState.saturation(/*phaseIdx=*/1);
        if (Sg > 0 && Sl > 0) {
            phasePresence_ = bothPhases;
            if (formulation == plSg)
                (*this)[switchIdx] = Sg;
            else if (formulation == pgSl)
                (*this)[switchIdx] = Sl;
            else {
                // invalid formulation
                assert(false);
            }
        }
        else if (Sg > 0) { 
            phasePresence_ = gPhaseOnly;
            (*this)[switchIdx] = fluidState.massFraction(/*phaseIdx=*/1, /*compIdx=*/0);
        }
        else { 
            // only liquid present
            assert(Sl > 0);

            phasePresence_ = lPhaseOnly;
            (*this)[switchIdx] = fluidState.massFraction(/*phaseIdx=*/0, /*compIdx=*/1);
        }
    }
    
protected:
    short phasePresence_;
    bool wasSwitched_;
};

} // end namepace

#endif
