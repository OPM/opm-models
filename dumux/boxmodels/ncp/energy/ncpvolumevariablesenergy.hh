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
 * \brief Contains the energy part of volume variables of the
 *        compositional NCP model.
 */
#ifndef DUMUX_NCP_VOLUME_VARIABLES_ENERGY_HH
#define DUMUX_NCP_VOLUME_VARIABLES_ENERGY_HH

#include <dumux/boxmodels/ncp/ncpproperties.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <string>

namespace Dumux {

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the compositional NCP model.
 *
 * This is the dummy class for the isothermal case. Note that we're
 * only isothermal in the sense that the temperature at a location and
 * a time is specified outside of the model!
 */
template <class TypeTag, bool enableEnergy/*=false*/>
class NcpVolumeVariablesEnergy
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    template <class FluidState>
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
    template <class FluidState>
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
     * \brief Returns the name of a primary variable or an empty
     *        string if the specified primary variable index does not belong to
     *        the energy module.
     */
    static std::string primaryVarName(int pvIdx)
    { return ""; };

    /*!
     * \brief Returns the name of an equation or an empty
     *        string if the specified equation index does not belong to
     *        the energy module.
     */
    static std::string eqName(int eqIdx)
    { return ""; };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     */
    static Scalar primaryVarWeight(const Model &model, int globalVertexIdx, int pvIdx)
    { return -1; }

    /*!
     * \brief Returns the relative weight of a equation.
     */
    static Scalar eqWeight(const Model &model, int globalVertexIdx, int eqIdx)
    { return -1; }

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


    /*!
     * \brief Set the enthalpy rate on a boundary rate vector
     */
    template <class FluidState>
    static void enthalpyBoundaryFlux(BoundaryRateVector &rate,
                                     const FluxVariables &fluxVars,
                                     const VolumeVariables &insideVolVars,
                                     const FluidState &fs,
                                     const typename FluidSystem::ParameterCache &paramCache,
                                     int phaseIdx,
                                     Scalar density)
    { }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolid() const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "heatCapacitySolid() called with the energy equation being disabled");
    }

    /*!
     * \brief Returns the lumped heat conducitvity [W/(K m)] of
     *        the porous medium and the fluid phases.
     */
    Scalar heatConductivity() const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "heatCapacitySolid() called with the energy equation being disabled");
    }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    { }
};

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the compositional NCP model.
 */
template <class TypeTag>
class NcpVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { temperatureIdx = Indices::temperatureIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    template <class FluidState>
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
    template <class FluidState>
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
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     */
    void globalPhaseStorage(EqVector &dest, int phaseIdx)
    {
        dest = 0;
        EqVector tmp;

        ElementContext elemCtx(this->problem_());
        auto elemIt = this->gridView_().template begin<0>();
        const auto &elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            const auto &fvElemGeom = elemCtx.fvElemGeom(/*timeIdx=*/0);
            
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                if (elemIt->partitionType() != Dune::InteriorEntity)
                    continue; // ignore ghost and overlap elements

                tmp = 0;
                this->localResidual().addPhaseStorage(tmp, 
                                                      elemCtx,
                                                      scvIdx,
                                                      /*timeIdx=*/0,
                                                      phaseIdx);
                tmp *= 
                    fvElemGeom.subContVol[scvIdx].volume
                    * elemCtx.volVars(scvIdx, /*timeIdx=*/0).extrusionFactor();
                dest += tmp;
            }
        };

        dest = this->gridView_().comm().sum(dest);
    }


    /*!
     * \brief Returns the name of a primary variable or an empty
     *        string if the specified primary variable index does not belong to
     *        the energy module.
     */
    static std::string primaryVarName(int pvIdx)
    { 
        if (pvIdx == temperatureIdx)
            return "temperature";
        return "";
    };

    /*!
     * \brief Returns the name of an equation or an empty
     *        string if the specified equation index does not belong to
     *        the energy module.
     */
    static std::string eqName(int eqIdx)
    { 
        if (eqIdx == energyEqIdx)
            return "energy";
        return "";
    };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     */
    static Scalar primaryVarWeight(const Model &model, int globalVertexIdx, int pvIdx)
    {
        if (pvIdx != temperatureIdx)
            return -1;

        // make the weight of the temperature primary variable inversly proportional to its value
        return std::max(1.0/1000, 1.0/model.solution(/*timeIdx=*/0)[globalVertexIdx][temperatureIdx]);
    }

    /*!
     * \brief Returns the relative weight of a equation.
     */
    static Scalar eqWeight(const Model &model, int globalVertexIdx, int eqIdx)
    {
        if (eqIdx != energyEqIdx)
            return -1;

        // approximate heat capacity of 1kg of air
        return 1.0/1.0035e3;
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
     * \brief Set the enthalpy rate per second of a rate vector, .
     */
    static void setEnthalpyRate(RateVector &rateVec, Scalar rate)
    {
        rateVec[energyEqIdx] = rate;
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
    }

    /*!
     * \brief Set the enthalpy rate on a boundary rate vector
     */
    template <class FluidState>
    static void enthalpyBoundaryFlux(BoundaryRateVector &rate,
                                     const FluxVariables &fluxVars,
                                     const VolumeVariables &insideVolVars,
                                     const FluidState &fs,
                                     const typename FluidSystem::ParameterCache &paramCache,
                                     int phaseIdx,
                                     Scalar density)
    {
        Scalar enthalpy =
            (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
            ? FluidSystem::enthalpy(fs, paramCache, phaseIdx)
            : insideVolVars.fluidState().enthalpy(phaseIdx);
        
        rate[energyEqIdx] += 
            fluxVars.volumeFlux(phaseIdx)
            * enthalpy
            * density;
    }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    Scalar heatCapacitySolid() const
    { return heatCapacitySolid_; }

    /*!
     * \brief Returns the lumped thermal conductivity [W/(K m)] of the
     *        medium.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(heatCapacitySolid_);
        Valgrind::CheckDefined(heatConductivity_);
    }

protected:
    Scalar heatCapacitySolid_;
    Scalar heatConductivity_;
};

} // end namepace

#endif
