// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Katherina Baber
 *   Copyright (C) 2008-2009 by Onur Dogan                                   *
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
 * \brief This file contains the data which is required to calculate
 *        the flux of the fluid over a face of a finite volume for the one-phase model.
 *
 *        This means pressure and temperature gradients, phase densities at
 *           the integration point, etc.
 */
#ifndef DUMUX_1P_FLUX_VARIABLES_HH
#define DUMUX_1P_FLUX_VARIABLES_HH

#include "1pproperties.hh"

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the flux of the fluid over a face of a
 *        finite volume for the one-phase model.
 */
template <class TypeTag>
class OnePFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
            

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;
    

    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief Caclulates the quantities required on a sub-control
     *        volume face for the 2p box model.
     */
    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        insideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].j;

        extrusionFactor_ =
            (elemCtx.volVars(insideScvIdx_).extrusionFactor() 
             + elemCtx.volVars(outsideScvIdx_).extrusionFactor()) / 2;
        
        ///////////////
        // calculate the pressure gradient
        ///////////////
        potentialGrad_  = 0.0;
        
        typedef typename FVElementGeometry::SubControlVolumeFace Scvf;
        const Scvf &scvf = elemCtx.fvElemGeom().subContVolFace[scvfIdx];
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx)
        {
            // FE gradient at vertex idx
            const Vector &feGrad = scvf.grad[scvIdx];
            const auto &fs = elemCtx.volVars(scvIdx, /*historyIdx=*/0).fluidState();

            Vector tmp(feGrad);
            tmp *= fs.pressure(/*phaseIdx=*/0);
            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(elemCtx.problem().gravity(elemCtx, insideScvIdx_));
            g += elemCtx.problem().gravity(elemCtx, outsideScvIdx_);
            g /= 2;
            
            const auto &fsI = elemCtx.volVars(insideScvIdx_, /*historyIdx=*/0).fluidState();
            const auto &fsJ = elemCtx.volVars(outsideScvIdx_, /*historyIdx=*/0).fluidState();

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar rhoI = fsI.density(phaseIdx);
                Scalar rhoJ = fsJ.density(phaseIdx);
                Scalar density = (rhoI + rhoJ)/2.0;

                // make gravity acceleration a force
                Vector f(g);
                f *= density;
        
                // calculate the final potential gradient
                potentialGrad_ -= f;
            }
        }

        ///////////////
        // Calculate the velocity
        ///////////////
        const SpatialParameters &spatialParams = elemCtx.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                insideScvIdx_),
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                outsideScvIdx_));

        const Vector &normal = elemCtx.fvElemGeom().subContVolFace[scvfIdx].normal;

        // calculate the flux in the normal direction of the
        // current sub control volume face:
        //
        // v = - (K grad p) * n
        //
        // (the minus comes from the Darcy law which states that
        // the flux is from high to low pressure potentials.)                           
        K.mv(potentialGrad_, filterVelocity_);
        // velocities go along negative pressure gradients
        filterVelocity_ *= -1;

        // scalar product with the face normal
        filterVelocityNormal_ = 0.0;
        for (int i = 0; i < Vector::size; ++i) 
            filterVelocityNormal_ += filterVelocity_[i]*normal[i];
        
        // multiply both with the upstream mobility
        const auto &up = elemCtx.volVars(upstreamIdx(/*phaseIdx=*/0), /*historyIdx=*/0);
        filterVelocityNormal_ *= up.mobility(/*phaseIdx=*/0);
        filterVelocity_ *= up.mobility(/*phaseIdx=*/0);
    }


    /*!
     * \brief Return the extrusion factor of the SCVF.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return the fluid's pressure potential gradient [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return potentialGrad_;
    }

    /*!
     * \brief Return the fluid's filter velocity.
     *
     * For low Reynolds numbers that's the Darcy velocity, for higher
     * ones advanced velocity models like the one suggested by
     * Forchheimer kick in.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &filterVelocity(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return filterVelocity_;
    }

    /*!
     * \brief Return a phase's pressure potential gradient times
     *        intrinsic permeability times times the mobility timesthe
     *        normal of the sub control volume face times the area of
     *        the SCVF.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar filterVelocityNormal(int phaseIdx) const
    { 
        assert(phaseIdx == 0); // this is a single phase model!
        return filterVelocityNormal_;
    }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "inside" of the sub-control volume face.
     */
    short insideIdx() const
    { return insideScvIdx_; }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "outside" of the sub-control volume face.
     */
    short outsideIdx() const
    { return outsideScvIdx_; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIdx(int phaseIdx) const
    { 
        assert(phaseIdx == 0); // this is a single phase model!
        return (filterVelocityNormal_ > 0)?outsideScvIdx_:insideScvIdx_;
    }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIdx(int phaseIdx) const
    {
        assert(phaseIdx == 0); // this is a single phase model!
        return (filterVelocityNormal_ > 0)?insideScvIdx_:outsideScvIdx_;
    }

    /*!
     * \brief Return the weight of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar upstreamWeight(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 0.0; }

private:
    // local indices of the inside and the outside sub-control volumes
    short insideScvIdx_;
    short outsideScvIdx_;

    // extrusion factor
    Scalar extrusionFactor_;

    // potential gradient
    Vector potentialGrad_;

    // filter velocity
    Vector filterVelocity_;

    // normal velocity
    Scalar filterVelocityNormal_;
};

} // end namepace

#endif
