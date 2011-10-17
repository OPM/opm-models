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
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int historyIdx)
    {
        ParentType::update(elemCtx,
                           scvIdx,
                           historyIdx);

        completeFluidState(fluidState_, elemCtx, scvIdx, historyIdx);

        const SpatialParameters &spatialParams = 
            elemCtx.problem().spatialParameters();

        // material law parameters
        const MaterialLawParams &materialParams =
            spatialParams.materialLawParams(elemCtx, scvIdx);

        relativePermeability_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx));
        relativePermeability_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx));


        // porosity
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

        int formulation = GET_PROP_VALUE(TypeTag, Formulation);
        if (formulation == Indices::pwSn) {
            Scalar Sn = priVars[Indices::saturationIdx];
            fluidState.setSaturation(nPhaseIdx, Sn);
            fluidState.setSaturation(wPhaseIdx, 1 - Sn);

            Scalar pW = priVars[Indices::pressureIdx];
            fluidState.setPressure(wPhaseIdx, pW);
            fluidState.setPressure(nPhaseIdx,
                                   pW + MaterialLaw::pC(materialParams, 1 - Sn));
        }
        else if (formulation == Indices::pnSw) {
            Scalar Sw = priVars[Indices::saturationIdx];
            fluidState.setSaturation(wPhaseIdx, Sw);
            fluidState.setSaturation(nPhaseIdx, 1 - Sw);

            Scalar pN = priVars[Indices::pressureIdx];
            fluidState.setPressure(nPhaseIdx, pN);
            fluidState.setPressure(wPhaseIdx,
                                   pN - MaterialLaw::pC(materialParams, Sw));
        }

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
                                        historyIdx);
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
    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx));
    }

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
