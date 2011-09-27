// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase model.
 */
#ifndef DUMUX_2PNI_VOLUME_VARIABLES_HH
#define DUMUX_2PNI_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/2p/2pvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase model.
 */
template <class TypeTag>
class TwoPNIVolumeVariables : public TwoPVolumeVariables<TypeTag>
{
    typedef TwoPVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    typedef typename ParentType::FluidState FluidState;
    
    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     *
     */
    Scalar internalEnergy(int phaseIdx) const
    { return this->fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     *  \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return this->fluidState_.enthalpy(phaseIdx); };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent..
    friend class TwoPVolumeVariables<TypeTag>;
    
    /*!
     * \brief Update the temperature for a given control volume.
     */
    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int historyIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx));
    }

    template <class ParameterCache>
    static void updateEnthalpy_(FluidState &fluidState,
                                const ParameterCache &paramCache,
                                const ElementContext &elemCtx,
                                int scvIdx,
                                int historyIdx)
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
                       int historyIdx)
    {
        // compute and set the heat capacity of the solid phase
        heatCapacity_ = elemCtx.problem().spatialParameters().heatCapacity(elemCtx, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);
    }

    Scalar heatCapacity_;
};

} // end namepace

#endif
