// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief This class calculates the pressure potential gradients and
 *        the filter velocities for multi-phase flow in porous media
 */
#ifndef DUMUX_BOX_MULTIPHASE_FLUX_VARIABLES_HH
#define DUMUX_BOX_MULTIPHASE_FLUX_VARIABLES_HH

#include <dumux/common/spline.hh>
#include <dumux/common/propertysystem.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(VolumeVariables);
NEW_PROP_TAG(EnableSmoothUpwinding);


NEW_PROP_TAG(UseTwoPointGradients);
NEW_PROP_TAG(NumPhases);
}

/*!
 * \ingroup BoxFluxVariables
 * \brief This class calculates the pressure potential gradients and
 *        the filter velocities for multi-phase flow in porous media
 */
template <class TypeTag>
class BoxMultiPhaseFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    enum {
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        useTwoPointGradients = GET_PROP_VALUE(TypeTag, UseTwoPointGradients)
    };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
public:
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        insideScvIdx_ = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].j;

        extrusionFactor_ =
            (elemCtx.volVars(insideScvIdx_, timeIdx).extrusionFactor()
             + elemCtx.volVars(outsideScvIdx_, timeIdx).extrusionFactor()) / 2;

        // update the base module (i.e. advection)
        calculateGradients_(elemCtx, scvfIdx, timeIdx);
        calculateVelocities_(elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Return true iff a fluid phase ought is used by the model
     */
    bool usePhase(int phaseIdx)
    { return true; }

    /*!
     * \brief Returns th extrusion factor for the sub-control volume face
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase
     *        at the face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &filterVelocity(int phaseIdx) const
    { return filterVelocity_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point times the phase normal times
     *        the face area [1/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar filterVelocityNormal(int phaseIdx) const
    { return filterVelocityNormal_[phaseIdx]; }

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
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIdx(int phaseIdx) const
    { return (filterVelocityNormal_[phaseIdx] > 0)?insideScvIdx_:outsideScvIdx_; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIdx(int phaseIdx) const
    { return (filterVelocityNormal_[phaseIdx] > 0)?outsideScvIdx_:insideScvIdx_; }

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
    { return 1.0 - upstreamWeight(phaseIdx); }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx,
                             int timeIdx)
    {
        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        
        if (useTwoPointGradients) {
            const auto &fsI = elemCtx.volVars(insideScvIdx_, timeIdx).fluidState();
            const auto &fsJ = elemCtx.volVars(outsideScvIdx_, timeIdx).fluidState();
            const auto &scvI = elemCtx.fvElemGeom(timeIdx).subContVol[insideScvIdx_];
            const auto &scvJ = elemCtx.fvElemGeom(timeIdx).subContVol[outsideScvIdx_];
            
            // distance between the centers of the two SCVs
            Scalar dist = 0;
            for (int i = 0; i < dimWorld; ++ i) {
                Scalar tmp = scvI.global[i] - scvJ.global[i];
                dist += tmp*tmp;
            }
            dist = std::sqrt(dist);
            
            // the "normalized normal" of the scvf divided by the
            // distance of the centers of the two adjacent SCVs
            Vector n = scvf.normal;
            n /= scvf.normal.two_norm()*dist;
            
            // calculate the pressure gradient
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!asImp_().usePhase(phaseIdx)) {
                    potentialGrad_[phaseIdx] = 0;
                    continue;
                }

                potentialGrad_[phaseIdx] = n;
                potentialGrad_[phaseIdx] *= (fsJ.pressure(phaseIdx) - fsI.pressure(phaseIdx));
            }
        }
        else {
            // reset all gradients to 0
            for (int phase = 0; phase < numPhases; ++phase) {
                potentialGrad_[phase] = Scalar(0);
            }

            // calculate gradients
            for (int scvIdx = 0;
                 scvIdx < elemCtx.numScv();
                 scvIdx ++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const Vector &feGrad = scvf.grad[scvIdx];
                const auto &fluidState = elemCtx.volVars(scvIdx, timeIdx).fluidState();

                // compute sum of pressure gradients for each phase
                for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
                {
                    if (!asImp_().usePhase(phaseIdx))
                        continue;

                    // the pressure gradient
                    Vector tmp(feGrad);
                    tmp *= fluidState.pressure(phaseIdx);
                    potentialGrad_[phaseIdx] += tmp;
                }
            }
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(elemCtx.problem().gravity(elemCtx, insideScvIdx_, timeIdx));
            g += elemCtx.problem().gravity(elemCtx, outsideScvIdx_, timeIdx);
            g /= 2;

            const auto &fsIn = elemCtx.volVars(insideScvIdx_, timeIdx).fluidState();
            const auto &fsOut = elemCtx.volVars(outsideScvIdx_, timeIdx).fluidState();
            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                if (!asImp_().usePhase(phaseIdx))
                    continue;

                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = fsIn.saturation(phaseIdx);
                Scalar SJ = fsOut.saturation(phaseIdx);
                Scalar rhoI = fsIn.density(phaseIdx);
                Scalar rhoJ = fsOut.density(phaseIdx);
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
                potentialGrad_[phaseIdx] -= f;
            }
        }
    }

    void calculateVelocities_(const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx)
    {
        const auto &problem = elemCtx.problem();

        // calculate the intrinsic permeability
        Tensor K;
        problem.meanK(K,
                      problem.intrinsicPermeability(elemCtx,
                                                    insideScvIdx_,
                                                    timeIdx),
                      problem.intrinsicPermeability(elemCtx,
                                                    outsideScvIdx_,
                                                    timeIdx));
        
        Vector normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;
        Scalar scvfArea = normal.two_norm();
        normal /= scvfArea;

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            if (!asImp_().usePhase(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0;
                filterVelocityNormal_[phaseIdx] = 0;
                upstreamScvIdx_[phaseIdx] = insideScvIdx_;
                downstreamScvIdx_[phaseIdx] = outsideScvIdx_;
                continue;
            }

            // calculate the "prelimanary" filter velocity not
            // taking the mobility into account
            K.mv(potentialGrad_[phaseIdx], filterVelocity_[phaseIdx]);
            filterVelocity_[phaseIdx] *= -1;
            filterVelocityNormal_[phaseIdx] = (filterVelocity_[phaseIdx] * normal)*scvfArea;

            // determine the upstream index. since this is a
            // semi-smooth non-linear solver, make upstream only look
            // at the evaluation point for the upstream decision
            const auto &evalFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);
            if (&evalFluxVars == this || !GET_PARAM(TypeTag, bool, EnableSmoothUpwinding)) {
                upstreamScvIdx_[phaseIdx] = (filterVelocityNormal_[phaseIdx]>0)?insideScvIdx_:outsideScvIdx_;
                downstreamScvIdx_[phaseIdx] = (filterVelocityNormal_[phaseIdx]>0)?outsideScvIdx_:insideScvIdx_;
            }
            else {
                upstreamScvIdx_[phaseIdx] = evalFluxVars.upstreamIdx(phaseIdx);
                downstreamScvIdx_[phaseIdx] = evalFluxVars.downstreamIdx(phaseIdx);
            }

            // calculate the actual darcy velocities by multiplying
            // the current "filter velocity" with the upstream mobility
            if (!GET_PARAM(TypeTag, bool, EnableSmoothUpwinding)) {
                const VolumeVariables &up = elemCtx.volVars(upstreamIdx(phaseIdx), timeIdx);
                filterVelocity_[phaseIdx] *= up.mobility(phaseIdx);
                filterVelocityNormal_[phaseIdx] *= up.mobility(phaseIdx);
            }
            else {
                handleSmoothUpwinding_(elemCtx, scvfIdx, timeIdx, phaseIdx, normal);
                filterVelocityNormal_[phaseIdx] = (filterVelocity_[phaseIdx]*normal)*scvfArea;
            }
        }
    }

    void handleSmoothUpwinding_(const ElementContext &elemCtx,
                                int scvfIdx,
                                int timeIdx,
                                int phaseIdx,
                                const Vector &normal)
    {
        const VolumeVariables &up = elemCtx.volVars(upstreamIdx(phaseIdx), timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(downstreamIdx(phaseIdx), timeIdx);
        
        // first, calculate the component of the "prelimary velocity"
        // which is parallel to the normal of the sub-control volume
        // face
        Vector parV(normal);
        parV *= normal * filterVelocity_[phaseIdx];
        Scalar x = parV.two_norm();
        assert(x >= 0);

        if (x == 0.0) {
            filterVelocity_[phaseIdx] = 0.0;
            return;
        }

        // slopes of the velocity in upwind, downwind direction and at
        // x = 0
        Scalar mUp = up.mobility(phaseIdx);
        Scalar mDn = dn.mobility(phaseIdx);
        Scalar m0 = Dumux::harmonicMean(mUp, mDn);
        
        Scalar maxMob = std::max(mUp, mDn);
        if (maxMob < 1e-8) {
            filterVelocity_[phaseIdx] = 0.0;
            return;
        }

        // put the mean viscosity and permeanbility in
        // relation to the viscosity of water at
        // approximatly 20 degrees Celsius.
        const Scalar pGradRef = 100; // [Pa/m]
        const Scalar mobRef = 1.0/1e-3; // [1 / (Pa s)]
        const Scalar KRef = 1e-12; // [m^2] = approx 1 Darcy
        const Scalar vRef = mobRef * KRef * pGradRef; // [m/s]

        // calculate the velocity below which smooth upwinding should
        // kick in. For this, we assume that the reference medium is
        // fully saturated with liquid water.
        Scalar eps = vRef / maxMob;
        if (x > eps) {
            // we only do tricks if x is below the epsilon
            // value. Here, this is not the case, so we use
            // the full upwinding scheme
            filterVelocity_[phaseIdx] *= mUp;
        }
        else {
            // interpolate between zero and epsilon using a cubic
            // spline
            Spline<Scalar> sp(/*x0=*/0.0,
                              /*x1=*/eps,
                              /*y0=*/0.0,
                              /*y1=*/mUp*eps,
                              /*m0=*/m0,
                              /*m1=*/mUp);
            
            // set the length of the velocity component which is
            // parallel to the face normal to the one from the spline,
            // and do not modify the perpendicular part
            Vector perpV = filterVelocity_[phaseIdx];
            perpV -= parV;
            perpV *= mUp;

            parV *= sp.eval(x)/parV.two_norm();

            filterVelocity_[phaseIdx] = perpV;
            filterVelocity_[phaseIdx] += parV;
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    // extrusion factor for the sub-control volume face
    Scalar extrusionFactor_;

    // local indices of the inside and the outside sub-control volumes
    short insideScvIdx_;
    short outsideScvIdx_;

    short upstreamScvIdx_[numPhases];
    short downstreamScvIdx_[numPhases];

    // pressure potential gradients of all phases
    Vector potentialGrad_[numPhases];

    // filter velocities of all phases
    Vector filterVelocity_[numPhases];

    // normal velocities, i.e. filter velocity times face normal times
    // face area
    Scalar filterVelocityNormal_[numPhases];
};

} // end namepace

#endif
