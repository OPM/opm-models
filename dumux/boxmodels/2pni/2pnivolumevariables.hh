// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
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
 *        finite volume in the non-isothermal two-phase
 *        model.
 */
#ifndef DUMUX_2PNI_VOLUME_VARIABLES_HH
#define DUMUX_2PNI_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/2p/2pvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
template <class TypeTag>
class TwoPNIVolumeVariables : public TwoPVolumeVariables<TypeTag>
{
    typedef TwoPVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { temperatureIdx = Indices::temperatureIdx };
    enum { energyEqIdx = Indices::energyEqIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    typedef typename ParentType::FluidState FluidState;

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacitySolid() const
    { return heatCapacitySolid_; };

    /*!
     * \brief Returns the total conductivity capacity
     *        \f$\mathrm{[W/m^2 / (K/m)]}\f$ of the rock matrix in the
     *        sub-control volume.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; };

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    {
        priVars[temperatureIdx] = fs.temperature(/*phaseIdx=*/0);
#ifndef NDEBUG
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(fs.temperature(/*phaseIdx=*/0) == fs.temperature(phaseIdx));
        }
#endif
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

protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent..
    friend class TwoPVolumeVariables<TypeTag>;

    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int timeIdx)
    {
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        fluidState.setTemperature(priVars[temperatureIdx]);
    }

    template <class ParameterCache>
    static void updateEnthalpy_(FluidState &fluidState,
                                const ParameterCache &paramCache,
                                const ElementContext &elemCtx,
                                int scvIdx,
                                int timeIdx)
    {
        // compute and set the internal energies of the fluid phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);

            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx)
    {
        // compute and set the heat capacity of the solid phase
        const auto &spatialParams = elemCtx.problem().spatialParameters();
        const auto &heatCondParams = spatialParams.heatConducionParams(elemCtx, scvIdx, timeIdx);

        heatCapacitySolid_ = spatialParams.heatCapacitySolid(elemCtx, scvIdx, timeIdx);
        heatConductivity_ = HeatConductionLaw::heatConductivity(heatCondParams, this->fluidState());

        Valgrind::CheckDefined(heatCapacitySolid_);
        Valgrind::CheckDefined(heatConductivity_);
    }
    Scalar heatCapacitySolid_;
    Scalar heatConductivity_;
};

} // end namepace

#endif
