// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2012 by Andreas Lauser

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
 * \copydoc Ewoms::RichardsVolumeVariables
 */
#ifndef EWOMS_RICHARDS_VOLUME_VARIABLES_HH
#define EWOMS_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <ewoms/disc/vcfv/vcfvvolumevariables.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 * \ingroup VcfvVolumeVariables
 *
 * \brief Volume averaged quantities required by the Richards model.
 */
template <class TypeTag>
class RichardsVolumeVariables
    : public VcfvVolumeVariables<TypeTag>,
      public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef VcfvVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pressureWIdx = Indices::pressureWIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { nPhaseIdx = 1 - wPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    //! The type returned by the fluidState() method
    typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \copydoc VcfvVolumeVariables::update
     */
    void update(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        assert(!FluidSystem::isLiquid(nPhaseIdx));

        ParentType::update(elemCtx, scvIdx, timeIdx);

        fluidState_.setTemperature(
            elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx));

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &problem = elemCtx.problem();
        const typename MaterialLaw::Params &materialParams
            = problem.materialLawParams(elemCtx, scvIdx, timeIdx);
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
        Scalar pW = priVars[pressureWIdx];
        Scalar pN = std::max(elemCtx.problem().referencePressure(elemCtx, scvIdx,
                                                                 /*timeIdx=*/0),
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
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, scvIdx, timeIdx);

        // update the quantities specific for the velocity model
        VelocityVolumeVariables::update_(elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleVolumeVariables::fluidState
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::porosity
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::intrinsicPermeability
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::relativePermeability
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::mobility
     */
    Scalar mobility(int phaseIdx) const
    {
        return relativePermeability(phaseIdx) / fluidState().viscosity(phaseIdx);
    }

private:
    FluidState fluidState_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
    Scalar porosity_;
};

} // namespace Ewoms

#endif
