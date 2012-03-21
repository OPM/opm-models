// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal compositional Stokes box model.
 *
 */
#ifndef DUMUX_STOKES2CNI_LOCAL_RESIDUAL_HH
#define DUMUX_STOKES2CNI_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokes2c/stokes2clocalresidual.hh>

#include "stokes2cnivolumevariables.hh"
#include "stokes2cnifluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal compositional Stokes box model. This class is derived
 *        from the stokes2c box local residual and adds the energy balance equation.
     *
     *  \param result The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
 */
template<class TypeTag>
class Stokes2cniLocalResidual : public Stokes2cLocalResidual<TypeTag>
{
    typedef Stokes2cLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cniIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dim = GridView::dimension };
    enum { energyIdx = Indices::energyIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };

public:
    /*!
     * \brief Evaluate the amount the additional quantities to the stokes2c model
     *        (energy equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(EqVector &result,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // compute the storage term for the transport equation
        ParentType::computeStorage(result, elemCtx, scvIdx, timeIdx);

        // compute the storage of energy
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
        result[energyIdx] =
            volVars.fluidState().density(phaseIdx) *
            volVars.fluidState().internalEnergy(phaseIdx);
    }

    /*!
     * \brief Evaluates the convective energy flux
     *        over a face of a sub-control volume and writes the result in
     *        the flux vector. This method is called by computeFlux in the base class.
     *
     * \param flux The vector for the fluxes over the SCV/boundary face
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int faceIdx,
                              int timeIdx,
                              bool onBoundary) const
    {
        // call computation of the advective fluxes of the stokes model
        // (momentum and mass fluxes)
        ParentType::computeAdvectiveFlux(flux, elemCtx, faceIdx, timeIdx, onBoundary);

        const auto &fluxVars = 
            onBoundary
            ? elemCtx.boundaryFluxVars(faceIdx, timeIdx)
            : elemCtx.fluxVars(faceIdx, timeIdx);

        // vertex data of the upstream and the downstream vertices
        const VolumeVariables &up = elemCtx.volVars(fluxVars.upstreamIdx(), timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(fluxVars.downstreamIdx(), timeIdx);

        const auto &fsUp = up.fluidState();
        const auto &fsDn = dn.fluidState();
        
        Scalar tmp = fluxVars.normalVelocityAtIP();

        tmp *=  
            this->massUpwindWeight_
            * fsUp.density(phaseIdx)
            * fsUp.enthalpy(phaseIdx) 
            +
            (1.0 - this->massUpwindWeight_)
            * fsDn.density(phaseIdx)
            * fsDn.enthalpy(phaseIdx);

        flux[energyIdx] += tmp;
        Valgrind::CheckDefined(flux[energyIdx]);
    }

    /*!
     * \brief Evaluates the conductive energy flux
     *        over the face of a sub-control volume.
     *
     * \param flux The vector for the fluxes over the SCV face
     * \param fluxVars The flux variables at the current SCV face
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int faceIdx,
                              int timeIdx,
                              bool onBoundary) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, elemCtx, faceIdx, timeIdx, onBoundary);
        
        const auto &fluxVars = 
            onBoundary
            ? elemCtx.boundaryFluxVars(faceIdx, timeIdx)
            : elemCtx.fluxVars(faceIdx, timeIdx);
        

        // diffusive heat flux
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            flux[energyIdx] -=
                fluxVars.temperatureGradAtIP()[dimIdx] *
                fluxVars.normal()[dimIdx] *
                fluxVars.heatConductivityAtIP();
    }
};

} // end namespace

#endif
