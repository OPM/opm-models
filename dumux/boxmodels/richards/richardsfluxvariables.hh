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
 * \brief Data which is required to calculate the flux of fluid over a
 *        face of a finite volume
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include "richardsproperties.hh"

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup BoxFluxVariables
 * \brief Calculates and stores the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVariables) ElementVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    
    enum { wPhaseIdx = Indices::wPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*
     * \brief The constructor
     */
    RichardsFluxVariables()
    {}

    /*!
     * \brief Caclulates the quantities required on a sub-control
     *        volume face for the 2p box model.
     */
    void update(const ElementVariables &elemVars, int scvfIdx)
    {
        insideScvIdx_ = elemVars.fvElemGeom().subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemVars.fvElemGeom().subContVolFace[scvfIdx].j;

        extrusionFactor_ =
            (elemVars.volVars(insideScvIdx_).extrusionFactor() 
             + elemVars.volVars(outsideScvIdx_).extrusionFactor()) / 2;

        calculateGradients_(elemVars, scvfIdx);
        calculateNormalFluxes_(elemVars, scvfIdx);
    };

    /*!
     * \brief Return the extrusion factor of the SCVF.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return a phase's pressure potential gradient.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx = wPhaseIdx) const
    { 
        assert(phaseIdx == wPhaseIdx);
        return potentialGrad_;
    }

    /*!
     * \brief Return a phase's pressure potential gradient times
     *        intrinsic permeability times the normal of the sub
     *        control volume face times the area of the SCVF.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar normalFlux(int phaseIdx = wPhaseIdx) const
    { 
        assert(phaseIdx == wPhaseIdx);
        return normalFlux_;
    }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    int downstreamIdx(int phaseIdx = wPhaseIdx) const
    { 
        assert(phaseIdx == wPhaseIdx);
        return (normalFlux_ > 0)?outsideScvIdx_:insideScvIdx_;
    }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    int upstreamIdx(int phaseIdx = wPhaseIdx) const
    { 
        assert(phaseIdx == wPhaseIdx);
        return (normalFlux_ > 0)?insideScvIdx_:outsideScvIdx_;
    }

protected:
    void calculateGradients_(const ElementVariables &elemVars,
                             int scvfIdx)
    {
        // reset all wetting phase potential gradient
        potentialGrad_ = Scalar(0);
        
        const auto &scvf = elemVars.fvElemGeom().subContVolFace[scvfIdx];

        // calculate gradients
        for (int scvIdx = 0;
             scvIdx < elemVars.numScv();
             scvIdx ++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = scvf.grad[scvIdx];

            // compute pressure gradient for the wetting phase
            Vector tmp(feGrad);
            tmp *= elemVars.volVars(scvIdx, /*historyIdx=*/0).pressure(wPhaseIdx);
            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(elemVars.problem().gravity(elemVars, insideScvIdx_));
            g += elemVars.problem().gravity(elemVars, outsideScvIdx_);
            g /= 2;
            
            
            // calculate the phase density at the integration point. we
            // only do this if the wetting phase is present in both cells
            Scalar SI = elemVars.volVars(insideScvIdx_, /*historyIdx=*/0).saturation(wPhaseIdx);
            Scalar SJ = elemVars.volVars(outsideScvIdx_, /*historyIdx=*/0).saturation(wPhaseIdx);
            Scalar rhoI = elemVars.volVars(insideScvIdx_, /*historyIdx=*/0).density(wPhaseIdx);
            Scalar rhoJ = elemVars.volVars(outsideScvIdx_, /*historyIdx=*/0).density(wPhaseIdx);
            Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
            Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
            if (fI + fJ == 0)
                // doesn't matter because no wetting phase is present in
                // both cells!
                fI = fJ = 0.5;
            Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);
            
            // make gravity acceleration a force
            Vector f(g);
            f *= density;

            // calculate the final potential gradient
            potentialGrad_ -= f;
        }
    }

    void calculateNormalFluxes_(const ElementVariables &elemVars, 
                                int scvfIdx)
    {
        const SpatialParameters &spatialParams = elemVars.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemVars,
                                                                insideScvIdx_),
                            spatialParams.intrinsicPermeability(elemVars,
                                                                outsideScvIdx_));

        const Vector &normal = elemVars.fvElemGeom().subContVolFace[scvfIdx].normal;

        // calculate the flux in the normal direction of the
        // current sub control volume face:
        //
        // v = - (K grad p) * n
        //
        // (the minus comes from the Darcy law which states that
        // the flux is from high to low pressure potentials.)
        Vector tmpVec;
        K.mv(potentialGrad(wPhaseIdx), tmpVec);
        
        // scalar product with the face normal
        normalFlux_ = 0.0;
        for (int i = 0; i < Vector::size; ++i) 
            normalFlux_ += tmpVec[i]*normal[i];

        // flux is along negative pressure gradients
        normalFlux_ *= -1;
    }

    // local indices of the inside and the outside sub-control volumes
    int insideScvIdx_;
    int outsideScvIdx_;

    // extrusion factor
    Scalar extrusionFactor_;

    // gradients
    Vector potentialGrad_;

    // normal fluxes
    Scalar normalFlux_;
};

} // end namepace

#endif
