// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Felix Bode                                        *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"

#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup BoxVolumeVariables
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model
 */
template <class TypeTag>
class RichardsVolumeVariables 
    : public BoxVolumeVariables<TypeTag>
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pwIdx = Indices::pwIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = Indices::wPhaseIdx };
    enum { nPhaseIdx = Indices::nPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables as a vector for the finite
     *                volume.
     * \param problem The physical problem which needs to be solved.
     * \param element The DUNE Codim<0> enitity which intersects
     *                the control volume of the box method
     * \param elemGeom The element's finite volume geometry
     * \param scvIdx The local index of the sub control volume inside the element
     * \param isOldSol Specifies whether the solution is from
     *                 the previous time step or from the current one
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        assert(!FluidSystem::isLiquid(nPhaseIdx));

        ParentType::update(elemCtx, scvIdx, timeIdx);

        fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx));

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &problem = elemCtx.problem();
        const typename MaterialLaw::Params &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        /////////
        // calculate the pressures
        /////////
                    
        // first, we have to find the minimum capillary pressure (i.e. Sw = 0)
        fluidState_.setSaturation(wPhaseIdx, 1.0);
        fluidState_.setSaturation(nPhaseIdx, 0.0);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
                    
        // non-wetting pressure can be larger than the
        // reference pressure if the medium is fully
        // saturated by the wetting phase
        Scalar pW = priVars[pwIdx];
        Scalar pN = std::max(elemCtx.problem().referencePressure(elemCtx, scvIdx, /*timeIdx=*/0),
                             pW + (pC[nPhaseIdx] - pC[wPhaseIdx]));
                    
        /////////
        // calculate the saturations
        /////////
        fluidState_.setPressure(wPhaseIdx, pW);
        fluidState_.setPressure(nPhaseIdx, pN);

        PhaseVector sat;
        MaterialLaw::saturations(sat, materialParams, fluidState_);
        fluidState_.setSaturation(wPhaseIdx, sat[wPhaseIdx]);
        fluidState_.setSaturation(nPhaseIdx, 1.0 - sat[wPhaseIdx]);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        // compute and set the wetting phase viscosity
        Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, wPhaseIdx);
        fluidState_.setViscosity(wPhaseIdx, mu);
        fluidState_.setViscosity(nPhaseIdx, 1e-20);

        // compute and set the wetting phase density
        Scalar rho = FluidSystem::density(fluidState_, paramCache, wPhaseIdx);
        fluidState_.setDensity(wPhaseIdx, rho);
        fluidState_.setDensity(nPhaseIdx, 1e-20);

        //////////
        // specify the other parameters
        //////////
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, scvIdx, timeIdx);
        
        // update the quantities specific for the velocity model
        VelocityVolumeVariables::update_(elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \brief Returns a reference to the fluid state for the volume
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the average porosity [] within the control volume.
     *
     * The porosity is defined as the ratio of the pore space to the
     * total volume, i.e. \f[ \Phi := \frac{V_{pore}}{V_{pore} + V_{rock}} \f]
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the intrinsic permeability tensor for the sub-control volume
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
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
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the non-wetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    { }

    /*!
     * \brief Set the enthalpy rate per second of a rate vector, .
     */
    static void setEnthalpyRate(RateVector &rateVec, Scalar rate)
    {}
    
    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class FluidState>
    static void setEnthalpyRate(RateVector &rateVec,
                                const FluidState &fluidState, 
                                int phaseIdx, 
                                Scalar volume)
    { }

protected:
    FluidState fluidState_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
    Scalar porosity_;
};

}

#endif
