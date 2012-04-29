// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 *        finite volume in the Stokes box model.
 */
#ifndef DUMUX_STOKES_VOLUME_VARIABLES_HH
#define DUMUX_STOKES_VOLUME_VARIABLES_HH

#include "stokesproperties.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the Stokes box model.
 */
template <class TypeTag>
class StokesVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, StokesIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };
    enum { pressureIdx = Indices::pressureIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*!
     * \copydoc BoxVolumeVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvIdx, timeIdx);

        asImp_().updateTemperature_(elemCtx, scvIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        fluidState_.setPressure(phaseIdx, priVars[pressureIdx]);

        // create NullParameterCache and do dummy update
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        fluidState_.setDensity(phaseIdx,
                               FluidSystem::density(fluidState_,
                                                    paramCache,
                                                    phaseIdx));
        fluidState_.setViscosity(phaseIdx,
                                 FluidSystem::viscosity(fluidState_,
                                                        paramCache,
                                                        phaseIdx));
        
        // compute and set the energy related quantities
        asImp_().updateEnergy_(paramCache, elemCtx, scvIdx, timeIdx);

        // the effective velocity of the control volume
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            velocity_[dimIdx] = priVars[Indices::velocity0Idx + dimIdx];

        // the gravitational acceleration applying to the material
        // inside the volume
        gravity_ = elemCtx.problem().gravity();
    }

    /*!
     * \brief Update the gradients for the sub-control volumes.
     */
    void updateScvGradients(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        // calculate the pressure gradient at the SCV using finite
        // element gradients
        pressureGrad_ = 0.0;
        for (int i = 0; i < elemCtx.numScv(); ++i) {
            const auto &feGrad = elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx].gradCenter[i];
            Vector tmp(feGrad);
            tmp *= elemCtx.volVars(i, timeIdx).fluidState().pressure(phaseIdx);
            
            pressureGrad_ += tmp;
        }            
    }


    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }
    
    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ of
     *        the fluid within the sub-control volume.
     */
    Scalar molarDensity() const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the velocity vector in the sub-control volume.
     */
    const Vector &velocity() const
    { return velocity_; }

    /*!
     * \brief Returns the pressure gradient in the sub-control volume.
     */
    const Vector &pressureGradient() const
    { return pressureGrad_; }

    /*!
     * \brief Returns the gravitational acceleration vector in the
     *        sub-control volume.
     */
    const Vector &gravity() const
    { return gravity_; } 

protected:
    template<class ParameterCache>
    void updateEnergy_(const ParameterCache &paramCache,
                       const ElementContext &elemCtx,
                       int scvIdx, int timeIdx)
    { }

    void updateTemperature_(const ElementContext &elemCtx,
                            int scvIdx, int timeIdx)
    {
        Scalar T = elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx);
        this->fluidState_.setTemperature(T);
    }

    Vector velocity_;
    Vector gravity_;
    Vector pressureGrad_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
