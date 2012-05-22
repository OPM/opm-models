// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief Contains the supplemental quantities, which are constant within a
 *        finite volume in the non-isothermal compositional Stokes box model.
 */
#ifndef DUMUX_STOKES_NI_VOLUME_VARIABLES_HH
#define DUMUX_STOKES_NI_VOLUME_VARIABLES_HH

#include <dumux/freeflow/stokes/stokesvolumevariables.hh>
#include "stokesniproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokesNIModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-component Stokes
 *        box model.
 */
template <class TypeTag>
class StokesNIVolumeVariables : public StokesVolumeVariables<TypeTag>
{
    typedef StokesVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { temperatureIdx = Indices::temperatureIdx };

public:
    /*!
     * \copydoc BoxVolumeVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvIdx, timeIdx);

        // set the heat conductivity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(this->fluidState_, phaseIdx);
        heatConductivity_ = FluidSystem::thermalConductivity(this->fluidState_,
                                                             paramCache,
                                                             phaseIdx);
    }

    /*!
     * \brief Returns the heat conductivity of the fluid phase in
     *        the sub-control volume.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; }


protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent...
    friend class StokesVolumeVariables<TypeTag>;

    template<class ParameterCache>
    void updateEnergy_(const ParameterCache& paramCache,
                       const ElementContext &elemCtx,
                       int scvIdx, int timeIdx)
    {
        Scalar h = FluidSystem::enthalpy(this->fluidState_, paramCache, phaseIdx);
        this->fluidState_.setEnthalpy(phaseIdx, h);
    }

    void updateTemperature_(const ElementContext &elemCtx,
                            int scvIdx, int timeIdx)
    {
        Scalar T = elemCtx.primaryVars(scvIdx, timeIdx)[temperatureIdx];
        this->fluidState_.setTemperature(T);
    }

    Scalar heatConductivity_;
};

} // end namespace

#endif
