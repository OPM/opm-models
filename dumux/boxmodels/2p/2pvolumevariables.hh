// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
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
 *        finite volume in the two-phase model.
 */
#ifndef DUMUX_2P_VOLUME_VARIABLES_HH
#define DUMUX_2P_VOLUME_VARIABLES_HH

#include "2pproperties.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPFluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
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

        asImp_().updateTemperature_(elemCtx, scvIdx, timeIdx);

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &problem = elemCtx.problem();
        const typename MaterialLaw::Params &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        Scalar Sn = priVars[Indices::saturationIdx];
        fluidState_.setSaturation(nPhaseIdx, Sn);
        fluidState_.setSaturation(wPhaseIdx, 1 - Sn);
        
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        
        Scalar pW = priVars[Indices::pressureIdx];
        fluidState_.setPressure(wPhaseIdx, pW);
        fluidState_.setPressure(nPhaseIdx,
                                pW + (pC[nPhaseIdx] - pC[wPhaseIdx]));

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);
        }

        asImp_().updateEnergy_(paramCache,
                               elemCtx,
                               scvIdx,
                               timeIdx);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // porosity
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(paramCache, elemCtx, scvIdx, timeIdx);
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

protected:
    void updateTemperature_(const ElementContext &elemCtx,
                            int scvIdx,
                            int timeIdx)
    { fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx)); }

    template<class ParameterCache>
    void updateEnergy_(const ParameterCache &paramCache,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx)
    { }

    FluidState fluidState_;
    Scalar porosity_;
    Scalar relativePermeability_[numPhases];

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
