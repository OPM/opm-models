// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan                                        *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Quantities required by the one-phase box model defined on a vertex.
 */
#ifndef DUMUX_1P_VOLUME_VARIABLES_HH
#define DUMUX_1P_VOLUME_VARIABLES_HH

#include "1pproperties.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>

#include <dumux/material/fluidstates/immisciblefluidstate.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class OnePVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, OnePIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

public:
    //! Type of the fluid state
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx, scvIdx, timeIdx);

        completeFluidState(fluidState_, elemCtx, scvIdx, timeIdx);

        // porosity
        const auto &problem = elemCtx.problem();
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(elemCtx, scvIdx, timeIdx);
    };

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int timeIdx)
    {
        Implementation::updateTemperature_(fluidState, elemCtx, scvIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        fluidState.setPressure(/*phaseIdx=*/0, priVars[Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
        }

        Implementation::updateEnthalpy_(fluidState,
                                        paramCache,
                                        elemCtx,
                                        scvIdx,
                                        timeIdx);
    }

    /*!
     * \brief Returns the thermodynamic state of the fluid.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the relative permeability of the fluid []
     *
     * This is always 1 for single phase flow.
     */
    Scalar relativePermeability(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return 1.0;
    };

    /*!
     * \brief Returns the mobility of the fluid [1 / (Pa s)]
     */
    Scalar mobility(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx);
    };

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
    {
    }
    
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
    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int timeIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx));
    }

    template<class ParameterCache>
    static void updateEnthalpy_(FluidState &fluidState,
                                const ParameterCache &paramCache,
                                const ElementContext &elemCtx,
                                int scvIdx,
                                int timeIdx)
    { }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx)
    { }

    FluidState fluidState_;
    Scalar porosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
