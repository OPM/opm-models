/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::PvsPrimaryVariables
 */
#ifndef EWOMS_PVS_PRIMARY_VARIABLES_HH
#define EWOMS_PVS_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <ewoms/models/common/energymodule.hh>
#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include "pvsindices.hh"
#include "pvsproperties.hh"

#include <iostream>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 *
 * \brief Represents the primary variables used in the primary
 *        variable switching compositional model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class PvsPrimaryVariables
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq)>
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
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Opm::NcpFlash<Scalar, FluidSystem> NcpFlash;

public:
    PvsPrimaryVariables() : ParentType()
    { Valgrind::SetDefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    explicit PvsPrimaryVariables(Scalar value) : ParentType(value)
    {
        Valgrind::CheckDefined(value);
        Valgrind::SetDefined(*this);

        phasePresence_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables &)
     */
    PvsPrimaryVariables(const PvsPrimaryVariables &value) : ParentType(value)
    {
        Valgrind::SetDefined(*this);

        phasePresence_ = value.phasePresence_;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
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
        Opm::CompositionalFluidState<Scalar, FluidSystem> fsFlash;

        // use the externally given fluid state as initial value for
        // the flash calculation
        fsFlash.assign(fluidState);

        // calculate the phase densities
        paramCache.updateAll(fsFlash);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fsFlash.setDensity(phaseIdx, FluidSystem::density(fsFlash, paramCache,
                                                              phaseIdx));
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx]
                    += fsFlash.saturation(phaseIdx)
                       * fsFlash.molarity(phaseIdx, compIdx);
            }
        }

        // run the flash calculation
        // NcpFlash::guessInitial(fsFlash, paramCache, globalMolarities);
        NcpFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams,
                                              globalMolarities);

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
     *
     * \param value The new phase presence. The phase with index i is
     *              present if the i-th bit of \c value is 1.
     */
    void setPhasePresence(short value)
    { phasePresence_ = value; }

    /*!
     * \brief Set whether a given indivividual phase should be present
     *        or not.
     *
     * \param phaseIdx The index of the phase which's presence ought to be set
     *or reset.
     * \param yesno If true, the presence of the phase is set, else it is reset
     */
    void setPhasePresent(int phaseIdx, bool yesno = true)
    {
        if (yesno)
            setPhasePresence(phasePresence_ | (1 << phaseIdx));
        else
            setPhasePresence(phasePresence_ & ~(1 << phaseIdx));
    }

    /*!
     * \brief Returns the index of the phase with's its saturation is
     *        determined by the closure condition of saturation.
     */
    unsigned char implicitSaturationIdx() const
    { return lowestPresentPhaseIdx(); }

    /*!
     * \brief Returns true iff a phase is present for a given phase
     *        presence.
     *
     * \param phaseIdx The index of the phase which's presence is
     *                 queried.
     * \param phasePresence The bit-map of present phases.
     */
    static bool phaseIsPresent(int phaseIdx, short phasePresence)
    { return phasePresence & (1 << phaseIdx); }

    /*!
     * \brief Returns true iff a phase is present for the current
     *        phase presence.
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    bool phaseIsPresent(int phaseIdx) const
    { return phasePresence_ & (1 << phaseIdx); }

    /*!
     * \brief Returns the phase with the lowest index that is present.
     */
    int lowestPresentPhaseIdx() const
    { return ffs(phasePresence_) - 1; }

    /*!
     * \brief Assignment operator from an other primary variables object
     */
    ThisType &operator=(const PrimaryVariables &value)
    {
        ParentType::operator=(value);
        phasePresence_ = value.phasePresence_;

        return *this;
    }

    /*!
     * \brief Assignment operator from a scalar value
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
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    Scalar explicitSaturationValue(int phaseIdx) const
    {
        if (!phaseIsPresent(phaseIdx) || phaseIdx == lowestPresentPhaseIdx())
            // non-present phases have saturation 0
            return 0.0;

        return (*this)[switch0Idx + phaseIdx - 1];
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
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
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // use a NCP condition to determine if the phase is
            // present or not
            Scalar a = 1;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                a -= fluidState.moleFraction(phaseIdx, compIdx);
            }
            Scalar b = fluidState.saturation(phaseIdx);

            if (b > a)
                phasePresence_ |= (1 << phaseIdx);
        }

        // some phase must be present
        if (phasePresence_ == 0)
            OPM_THROW(Opm::NumericalProblem,
                      "Phase state was 0, i.e., no fluid is present");

        // set the primary variables which correspond to mole
        // fractions of the present phase which has the lowest index.
        int lowestPhaseIdx = lowestPresentPhaseIdx();
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (switchIdx >= lowestPhaseIdx)
                ++phaseIdx;

            if (phaseIsPresent(phaseIdx)) {
                (*this)[switch0Idx + switchIdx]
                    = fluidState.saturation(phaseIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
            else {
                (*this)[switch0Idx + switchIdx]
                    = fluidState.moleFraction(lowestPhaseIdx, compIdx);
                Valgrind::CheckDefined((*this)[switch0Idx + switchIdx]);
            }
        }

        // set the mole fractions in of the remaining components in
        // the phase with the lowest index
        for (int compIdx = numPhases - 1; compIdx < numComponents - 1;
             ++compIdx) {
            (*this)[switch0Idx + compIdx]
                = fluidState.moleFraction(lowestPhaseIdx, compIdx + 1);
            Valgrind::CheckDefined((*this)[switch0Idx + compIdx]);
        }
    }

    /*!
     * \copydoc FlashPrimaryVariables::print
     */
    void print(std::ostream &os = std::cout) const
    {
        os << "(p_" << FluidSystem::phaseName(0) << " = "
           << this->operator[](pressure0Idx);
        int lowestPhaseIdx = lowestPresentPhaseIdx();
        for (int switchIdx = 0; switchIdx < numPhases - 1; ++switchIdx) {
            int phaseIdx = switchIdx;
            int compIdx = switchIdx + 1;
            if (phaseIdx >= lowestPhaseIdx)
                ++phaseIdx; // skip the saturation of the present
                            // phase with the lowest index

            if (phaseIsPresent(phaseIdx)) {
                os << ", S_" << FluidSystem::phaseName(phaseIdx) << " = "
                   << (*this)[switch0Idx + switchIdx];
            }
            else {
                os << ", x_" << FluidSystem::phaseName(lowestPhaseIdx) << "^"
                   << FluidSystem::componentName(compIdx) << " = "
                   << (*this)[switch0Idx + switchIdx];
            }
        };
        for (int compIdx = numPhases - 1; compIdx < numComponents - 1;
             ++compIdx) {
            os << ", x_" << FluidSystem::phaseName(lowestPhaseIdx) << "^"
               << FluidSystem::componentName(compIdx + 1) << " = "
               << (*this)[switch0Idx + compIdx];
        }
        os << ")";
        os << ", phase presence: " << static_cast<int>(phasePresence_);
    }

private:
    unsigned char phasePresence_;
};

} // namespace Ewoms

#endif
