// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
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
 * \brief Contains the energy part of volume variables of the M-phase,
 *        N-component model.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_ENERGY_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_ENERGY_HH

#include <dumux/boxmodels/mpnc/mpncproperties.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux {

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the two-phase, N-component model.
 *
 * This is the dummy class for the isothermal case. Note that we're
 * only isothermal in the sense that the temperature at a location and
 * a time is specified outside of the model!
 */
template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer /*=don't care*/>
class MPNCVolumeVariablesEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    //typedef typename GET_PROP_TYPE(TypeTag, MPNCEnergyIndices) EnergyIndices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const ElementContext &elemCtx,
                            int scvIdx,
                            int timeIdx) const
    {
        Scalar T = elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx);
        fs.setTemperature(T);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     *
     * Since we are isothermal, we don't need to do anything!
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        // the phase temperatures where already set by the base volume
        // variables!
    }


    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    {}                                    
    
    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class FluidState>
    static void setEnthalpyRate(RateVector &v,
                                const FluidState &fluidState, 
                                int phaseIdx, 
                                Scalar volume)
    { };


    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    static std::string primaryVarName(int pvIdx)
    {
        return "";
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    static std::string eqName(int eqIdx)
    {
        return "";
    };

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolid() const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "heatCapacitySolid() called with the energy equation being disabled");
    };

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    { }
};

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the two-phase, N-component model.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;

    typedef typename GET_PROP_TYPE(TypeTag, MPNCIndices) Indices;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    enum { temperatureIdx = Indices::temperatureIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    enum { numEnergyEqs = Indices::NumPrimaryEnergyVars};

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const ElementContext &elemCtx,
                            int scvIdx,
                            int timeIdx) const
    {
        // retrieve temperature from solution vector
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        fs.setTemperature(priVars[temperatureIdx]);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        const auto &problem = elemCtx.problem();

        Valgrind::SetUndefined(*this);

        // heat capacities of the fluids plus the porous medium
        heatCapacitySolid_ = problem.heatCapacitySolid(elemCtx, scvIdx, timeIdx);
        Valgrind::CheckDefined(heatCapacitySolid_);

        const auto &heatCondParams = problem.heatConductionParams(elemCtx, scvIdx, timeIdx);
        heatConductivity_ = HeatConductionLaw::heatConductivity(heatCondParams, fs);
        Valgrind::CheckDefined(heatConductivity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
            fs.setEnthalpy(phaseIdx, h);
        }
    }


    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    {
        priVars[temperatureIdx] = fs.temperature(/*phaseIdx=*/0);
    }

    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class FluidState>
    static void setEnthalpyRate(RateVector &rateVec,
                                const FluidState &fluidState, 
                                int phaseIdx, 
                                Scalar volume)
    {
        rateVec[energyEqIdx] = 
            volume
            * fluidState.density(phaseIdx)
            * fluidState.enthalpy(phaseIdx);
    };

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    static std::string primaryVarName(int pvIdx)
    {
        if (pvIdx == Indices::temperatureIdx)
            return "temperature";
        return "";
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    static std::string eqName(int eqIdx)
    {
        if (eqIdx == Indices::energyEqIdx)
            return "energy";
        return "";
    };

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolid() const
    { return heatCapacitySolid_; };

    /*!
     * \brief Returns the lumped thermal conductivity [W/(K m)] of the
     *        medium.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; };

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(heatCapacitySolid_);
        Valgrind::CheckDefined(heatConductivity_);
    };

protected:
    Scalar heatCapacitySolid_;
    Scalar heatConductivity_;
};

} // end namepace

#endif
