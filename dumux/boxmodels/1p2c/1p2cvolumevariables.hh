// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Quantities required by the single-phase, two-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1P2C_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "1p2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx
    };


    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \todo implement dispersivity
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx, scvIdx, timeIdx);

        asImp_().updateTemperature_(elemCtx, scvIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        fluidState_.setPressure(/*phaseIdx=*/0, priVars[pressureIdx]);

        Scalar x1 = priVars[x1Idx]; //mole or mass fraction of component 1
        fluidState_.setMoleFraction(/*phaseIdx=*/0, /*compIdx=*/0, 1 - x1);
        fluidState_.setMoleFraction(/*phaseIdx=*/0, /*compIdx=*/1, x1);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, /*phaseIdx=*/0);

        Scalar value;
        value = FluidSystem::density(fluidState_, paramCache, /*phaseIdx=*/0);
        fluidState_.setDensity(/*phaseIdx=*/0, value);
        value = FluidSystem::viscosity(fluidState_, paramCache, /*phaseIdx=*/0);
        fluidState_.setViscosity(/*phaseIdx=*/0, value);

        const auto &problem = elemCtx.problem();
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);
        tortuosity_ = problem.tortuosity(elemCtx, scvIdx, timeIdx);
        dispersivity_ = 0; // problem.dispersivity(elemCtx, scvIdx, timeIdx);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        paramCache.updatePhase(fluidState_, /*phaseIdx=*/0);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             /*phaseIdx=*/0,
                                                             /*comp0Idx=*/0,
                                                             /*comp1Idx=*/1);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(tortuosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(paramCache, elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the binary diffusion coefficient in the fluid
     */
    Scalar diffCoeff() const
    { return diffCoeff_; }

    /*!
     * \brief Returns the tortuosity of the streamlines of the fluid.
     */
    Scalar tortuosity() const
    { return tortuosity_; }

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
    }

    /*!
     * \brief Returns the mobility of the fluid [1 / (Pa s)]
     */
    Scalar mobility(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx);
    }

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

    Scalar porosity_;    //!< Effective porosity within the control volume
    Scalar tortuosity_;
    Scalar dispersivity_;
    Scalar diffCoeff_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}// end namepace

#endif
