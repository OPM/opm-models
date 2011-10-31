// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"

#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/boxmodels/common/boxvolumevariables.hh>

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
class RichardsVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    enum {
        pwIdx = Indices::pwIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

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
                int historyIdx)
    {
        assert(!FluidSystem::isLiquid(nPhaseIdx));

        ParentType::update(elemCtx, scvIdx, historyIdx);

        completeFluidState(fluidState_, elemCtx, scvIdx, historyIdx);

        //////////
        // specify the other parameters
        //////////
        const auto &spatialParams = elemCtx.problem().spatialParameters();
        const MaterialLawParams &matParams =
            spatialParams.materialLawParams(elemCtx, scvIdx);
        relativePermeabilityWetting_ =
            MaterialLaw::krw(matParams,
                             fluidState_.saturation(wPhaseIdx));

        porosity_ = spatialParams.porosity(elemCtx, scvIdx);

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(elemCtx, scvIdx, historyIdx);
    }

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        Implementation::updateTemperature_(fluidState, elemCtx, scvIdx, historyIdx);

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &spatialParams = elemCtx.problem().spatialParameters();
        const typename MaterialLaw::Params &materialParams =
            spatialParams.materialLawParams(elemCtx, scvIdx);
        const auto &priVars = elemCtx.primaryVars(scvIdx, historyIdx);

        // pressures
        Scalar minPc = MaterialLaw::pC(materialParams, 1.0);
        Scalar pnRef = elemCtx.problem().referencePressure(elemCtx, scvIdx);
        fluidState.setPressure(wPhaseIdx, priVars[pwIdx]);
        fluidState.setPressure(nPhaseIdx, std::max(pnRef, priVars[pwIdx] + minPc));

        // saturations
        Scalar Sw = MaterialLaw::Sw(materialParams,
                                    fluidState.pressure(nPhaseIdx)
                                    - fluidState.pressure(wPhaseIdx));
        fluidState.setSaturation(wPhaseIdx, Sw);
        fluidState.setSaturation(nPhaseIdx, 1 - Sw);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        int phaseIdx = wPhaseIdx;
        // compute and set the viscosity
        Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, mu);

        // compute and set the density
        Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, rho);

        Implementation::updateEnthalpy_(fluidState,
                                        paramCache,
                                        elemCtx,
                                        scvIdx,
                                        historyIdx);
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
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return relativePermeabilityWetting_;
        else
            return 1.0;
    }

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

protected:
    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx));
    };

    template<class ParameterCache>
    static void updateEnthalpy_(FluidState &fluidState,
                                const ParameterCache &paramCache,
                                const ElementContext &elemCtx,
                                int scvIdx,
                                int historyIdx)
    { }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const ElementContext &elemCtx,
                       int scvIdx,
                       int historyIdx)
    { }

    FluidState fluidState_;
    Scalar relativePermeabilityWetting_;
    Scalar porosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
