// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::BlackOilFluidState
 */
#ifndef EWOMS_BLACK_OIL_FLUID_STATE_HH
#define EWOMS_BLACK_OIL_FLUID_STATE_HH

#include "blackoilproperties.hh"

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 *
 * \brief Implements a "taylor-made" fluid state class for the black-oil model.
 *
 * I.e., it uses exactly the same quantities which are used by the ECL blackoil
 * model. Further quantities are computed "on the fly" and are accessing them is thus
 * relatively slow.
 */
template <class TypeTag>
class BlackOilFluidState
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

public:
    typedef Evaluation Scalar;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        Valgrind::CheckDefined(pvtRegionIdx_);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Valgrind::CheckDefined(saturation_[phaseIdx]);
            Valgrind::CheckDefined(pressure_[phaseIdx]);
            Valgrind::CheckDefined(invB_[phaseIdx]);
        }

        Valgrind::CheckDefined(Rs_);
        Valgrind::CheckDefined(Rv_);
#endif // NDEBUG
    }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        assert(false); // not yet implemented
    }

    void setPvtRegionIndex(unsigned short newPvtRegionIdx)
    { pvtRegionIdx_ = newPvtRegionIdx; }

    void setPressure(int phaseIdx, const Evaluation& p)
    { pressure_[phaseIdx] = p; }

    void setSaturation(int phaseIdx, const Evaluation& S)
    { saturation_[phaseIdx] = S; }

    void setInvB(int phaseIdx, const Evaluation& b)
    { invB_[phaseIdx] = b; }

    void setDensity(int phaseIdx, const Evaluation& rho)
    { density_[phaseIdx] = rho; }

    void setRs(const Evaluation& newRs)
    { Rs_ = newRs; }

    void setRv(const Evaluation& newRv)
    { Rv_ = newRv; }

    void setInvB(unsigned phaseIdx, const Evaluation& newb)
    { invB_[phaseIdx] = newb; }

    const Evaluation& pressure(unsigned phaseIdx) const
    { return pressure_[phaseIdx]; }

    const Evaluation& saturation(unsigned phaseIdx) const
    { return saturation_[phaseIdx]; }

    const Evaluation& temperature(unsigned phaseIdx) const
    { return temperature_; }

    const Evaluation& invB(unsigned phaseIdx) const
    { return invB_[phaseIdx]; }

    const Evaluation& Rs() const
    { return Rs_; }

    const Evaluation& Rv() const
    { return Rv_; }

    unsigned short pvtRegionIndex() const
    { return pvtRegionIdx_; };

    //////
    // slow methods
    //////
    Evaluation density(unsigned phaseIdx) const
    { return density_[phaseIdx]; }

    Evaluation molarDensity(unsigned phaseIdx) const
    {
        const auto& rho = density(phaseIdx);

        if (phaseIdx == waterPhaseIdx)
            return rho/FluidSystem::molarMass(waterCompIdx, pvtRegionIdx_);

        return
            rho*(moleFraction(phaseIdx, gasCompIdx)/FluidSystem::molarMass(gasCompIdx, pvtRegionIdx_)
                 + moleFraction(phaseIdx, oilCompIdx)/FluidSystem::molarMass(oilCompIdx, pvtRegionIdx_));

    }

    Evaluation molarVolume(unsigned phaseIdx) const
    { return 1.0/molarDensity(phaseIdx); }

    Evaluation viscosity(unsigned phaseIdx) const
    { return FluidSystem::viscosity(*this, phaseIdx, pvtRegionIdx_); }

    Evaluation enthalpy(unsigned phaseIdx) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The black-oil model does not support energy conservation yet.");
    }

    Evaluation internalEnergy(unsigned phaseIdx) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The black-oil model does not support energy conservation yet.");
    }

    Evaluation massFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_);
            }
            break;
        }

        OPM_THROW(std::logic_error,
                  "Invalid phase or component index!");
    }

    Evaluation moleFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_),
                                                          pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_),
                                                    pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_),
                                                    pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_),
                                                          pvtRegionIdx_);
            }
            break;
        }

        OPM_THROW(std::logic_error,
                  "Invalid phase or component index!");
    }

    Evaluation molarity(unsigned phaseIdx, unsigned compIdx) const
    { return moleFraction(phaseIdx, compIdx)*molarDensity(phaseIdx); }

    Evaluation averageMolarMass(unsigned phaseIdx) const
    {
        Evaluation result(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            result += FluidSystem::molarMass(compIdx, pvtRegionIdx_)*moleFraction(phaseIdx, compIdx);
        return result;
    }

    Evaluation fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return FluidSystem::fugacityCoefficient(*this, phaseIdx, compIdx, pvtRegionIdx_); }

    Evaluation fugacity(unsigned phaseIdx, unsigned compIdx) const
    {
        return
            fugacityCoefficient(phaseIdx, compIdx)
            *moleFraction(phaseIdx, compIdx)
            *pressure(phaseIdx);
    }

private:
    static const Evaluation temperature_;
    std::array<Evaluation, numPhases> pressure_;
    std::array<Evaluation, numPhases> saturation_;
    std::array<Evaluation, numPhases> invB_;
    std::array<Evaluation, numPhases> density_;
    Evaluation Rs_;
    Evaluation Rv_;
    unsigned short pvtRegionIdx_;
};

template <class TypeTag>
const typename BlackOilFluidState<TypeTag>::Evaluation BlackOilFluidState<TypeTag>::temperature_ =
    BlackOilFluidState<TypeTag>::FluidSystem::surfaceTemperature;

} // namespace Ewoms

#endif
