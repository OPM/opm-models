// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the black-oil model.
 */
#ifndef DUMUX_BLACK_OIL_VOLUME_VARIABLES_HH
#define DUMUX_BLACK_OIL_VOLUME_VARIABLES_HH

#include "blackoilproperties.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup BlackOilBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        
        pressure0Idx = Indices::pressure0Idx,
        saturation0Idx = Indices::saturation0Idx,

        wCompIdx = FluidSystem::wCompIdx,
        oCompIdx = FluidSystem::oCompIdx,
        gCompIdx = FluidSystem::gCompIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        oPhaseIdx = FluidSystem::oPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx
    };

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx,
                           scvIdx,
                           timeIdx);

        fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx));

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &problem = elemCtx.problem();
        const typename MaterialLaw::Params &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        // update the saturations
        Scalar sumSat = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++ phaseIdx) {
            fluidState_.setSaturation(phaseIdx, priVars[saturation0Idx + phaseIdx]);
            sumSat += priVars[saturation0Idx + phaseIdx];
        }
        fluidState_.setSaturation(numPhases - 1, 1 - sumSat);
        
        // update the pressures
        Scalar p0 = priVars[0];
        Scalar pC[numPhases];
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fluidState_.setPressure(phaseIdx, p0 + (pC[phaseIdx] - pC[0]));
        }
        
        // update phase compositions. first, set everything to 0, then
        // make the gas/water phases consist of only the gas/water
        // components and calculate the composition of the liquid oil
        // phase from the gas formation factor plus the gas/oil
        // formation volume factors and the reference densities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        // set composition of gas and water phases
        fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, 1.0);
        fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);

        // retrieve the relevant black-oil parameters from the fluid
        // system.
        Scalar Bg = FluidSystem::gasFormationVolumeFactor(fluidState_.pressure(oPhaseIdx));
        Scalar Bo = FluidSystem::oilFormationVolumeFactor(fluidState_.pressure(oPhaseIdx));
        Scalar Rs = FluidSystem::gasFormationFactor(fluidState_.pressure(oPhaseIdx));
        Scalar rhoo = FluidSystem::surfaceDensity(oPhaseIdx)/Bo;
        Scalar rhogref = FluidSystem::surfaceDensity(gPhaseIdx);
        Scalar MG = FluidSystem::molarMass(gPhaseIdx);
        Scalar MO = FluidSystem::molarMass(oPhaseIdx);

        // calculate composition of oil phase in terms of mass
        // fractions.
        Scalar XoG = Rs*rhogref / rhoo;
        Scalar XoO = 1 - XoG;

        // convert to mole fractions
        Scalar avgMolarMass = MO*MG/(MG + XoO*(MO - MG));
        Scalar xoG = XoG*avgMolarMass/MG;
        Scalar xoO = 1 - XoG;
        
        // finally set the oil-phase composition. yeah!
        fluidState_.setMoleFraction(oPhaseIdx, gCompIdx, xoG);
        fluidState_.setMoleFraction(oPhaseIdx, oCompIdx, xoO);

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
        }

        // set the gas and oil densities
        fluidState_.setDensity(oPhaseIdx, rhoo);
        fluidState_.setDensity(wPhaseIdx, FluidSystem::density(fluidState_, paramCache, wPhaseIdx));
        fluidState_.setDensity(gPhaseIdx, rhogref/Bg);
        
        // set the water density
        Scalar rhow = FluidSystem::density(fluidState_, paramCache, wPhaseIdx);
        fluidState_.setViscosity(wPhaseIdx, rhow);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // retrieve the porosity from the problem
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the relative permeability of a given phase
     *        within the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState().viscosity(phaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }


protected:
    FluidState fluidState_;
    Scalar porosity_;
    Scalar relativePermeability_[numPhases];
};

}

#endif
