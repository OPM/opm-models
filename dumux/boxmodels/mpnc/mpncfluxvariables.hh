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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_MPNC_FLUX_VARIABLES_HH
#define DUMUX_MPNC_FLUX_VARIABLES_HH

#include <dumux/common/spline.hh>

#include "diffusion/fluxvariables.hh"
#include "energy/mpncfluxvariablesenergy.hh"

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, three-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class MPNCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    enum {
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion),
        enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy),
        enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic),
        enableKineticEnergy = GET_PROP_VALUE(TypeTag, EnableKineticEnergy),
        enableGravity = GET_PROP_VALUE(TypeTag, EnableGravity)
    };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;


    typedef MPNCFluxVariablesDiffusion<TypeTag, enableDiffusion> FluxVariablesDiffusion;
    typedef MPNCFluxVariablesEnergy<TypeTag, enableEnergy, enableKineticEnergy> FluxVariablesEnergy;

public:
    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        insideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].j;

        extrusionFactor_ =
            (elemCtx.volVars(insideScvIdx_).extrusionFactor() 
             + elemCtx.volVars(outsideScvIdx_).extrusionFactor()) / 2;
        
        // update the base module (i.e. advection)
        calculateGradients_(elemCtx, scvfIdx);
        calculateVelocities_(elemCtx, scvfIdx);

        // update the flux data of the energy module (i.e. isothermal
        // or non-isothermal)
        energyVars_.update(elemCtx, scvfIdx);

        // update the flux data of the diffusion module (i.e. with or
        // without diffusion)
        diffusionVars_.update(elemCtx, scvfIdx);
    }

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
    { return upstreamWeight_[phaseIdx]; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 1.0 - upstreamWeight_[phaseIdx]; }

    ////////////////////////////////////////////////
    // forward calls to the diffusion module
    Scalar porousDiffCoeffL(int compIdx) const
    { return diffusionVars_.porousDiffCoeffL(compIdx); };

    Scalar porousDiffCoeffG(int compIIdx, int compJIdx) const
    { return diffusionVars_.porousDiffCoeffG(compIIdx, compJIdx); };

    const Scalar moleFraction(int phaseIdx, int compIdx) const
    { return diffusionVars_.moleFraction(phaseIdx, compIdx); };

    const Vector &moleFracGrad(int phaseIdx,
                               int compIdx) const
    { return diffusionVars_.moleFracGrad(phaseIdx, compIdx); };
    // end of forward calls to the diffusion module
    ////////////////////////////////////////////////

    /*!
     * \brief Returns the variables relevant for the energy module
     */
    const FluxVariablesEnergy &energyVars() const
    { return energyVars_; }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx)
    {
        // reset all gradients to 0
        for (int phase = 0; phase < numPhases; ++phase) {
            potentialGrad_[phase] = Scalar(0);
        }
        
        typedef typename FVElementGeometry::SubControlVolumeFace Scvf;
        const Scvf &scvf = elemCtx.fvElemGeom().subContVolFace[scvfIdx];

        // calculate gradients
        for (int scvIdx = 0;
             scvIdx < elemCtx.numScv();
             scvIdx ++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = scvf.grad[scvIdx];
            const auto &fluidState = elemCtx.volVars(scvIdx, /*historyIdx=*/0).fluidState();

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                Vector tmp(feGrad);
                tmp *= fluidState.pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }
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
            
            const auto &fsIn = elemCtx.volVars(insideScvIdx_, /*historyIdx=*/0).fluidState();
            const auto &fsOut = elemCtx.volVars(outsideScvIdx_, /*historyIdx=*/0).fluidState();
            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
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
                              int scvfIdx)
    {
        const SpatialParameters &spatialParams = 
            elemCtx.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                insideScvIdx_),
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                outsideScvIdx_));
        
        const Vector &normal = elemCtx.fvElemGeom().subContVolFace[scvfIdx].normal;

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the "prelimanary" filter velocity not
            // taking the mobility into account
            K.mv(potentialGrad_[phaseIdx], filterVelocity_[phaseIdx]);
            filterVelocity_[phaseIdx] *= -1;
            
            // determine upstream and downstream indices
            int upstreamIdx = outsideScvIdx_;
            int downstreamIdx = insideScvIdx_;
            if (filterVelocity_[phaseIdx] * normal > 0)
                std::swap(upstreamIdx, downstreamIdx);
            
            // calculate the actual darcy velocities and the upstream
            // and downstream weights.
            if (!GET_PARAM(TypeTag, bool, EnableSmoothUpwinding)) {
                const VolumeVariables &up = elemCtx.volVars(upstreamIdx);
                filterVelocity_[phaseIdx] *= up.mobility(phaseIdx);
                upstreamWeight_[phaseIdx] = 1.0;
            }
            else {
                const VolumeVariables &up = elemCtx.volVars(upstreamIdx);
                const VolumeVariables &dn = elemCtx.volVars(downstreamIdx);
                
                Scalar x = filterVelocity_[phaseIdx].two_norm();
                
                Scalar mUp = up.mobility(phaseIdx);
                Scalar mDn = dn.mobility(phaseIdx);
                Scalar m0 = Dumux::harmonicMean(mUp, mDn);
                
                // approximate the mean viscosity at the face
                Scalar meanVisc =
                    (up.fluidState().viscosity(phaseIdx) 
                     + dn.fluidState().viscosity(phaseIdx))
                    / 2;
                
                // put the mean viscosity and permeanbility in
                // relation to the viscosity of water at
                // approximatly 20 degrees Celsius.
                const Scalar pGradRef = 10; // [Pa/m]
                const Scalar muRef = 1e-3; // [Ns/m^2]
                const Scalar Kref = 1e-12; // [m^2] = approx 1 Darcy
                
                Scalar eps = pGradRef * Kref * meanVisc/muRef;
                if (0 <= x || eps < x) {
                    // we only do tricks if x is below the epsilon
                    // value. Here, this is not the case, so we use
                    // the full upwinding scheme
                    filterVelocity_[phaseIdx] *= mUp;
                    
                    upstreamWeight_[phaseIdx] = 1.0;
                }
                else {
                    // interpolate between zero and epsilon
                    Scalar xPos[] = { 0, eps };
                    Scalar yPos[] = { 0.0, mUp*eps };
                    Spline<Scalar> sp2(xPos, yPos, m0, mUp);
                    Scalar absV = sp2.eval(x);
                    Scalar vUp = x * mUp;
                    Scalar vDn = x * mDn;
                    
                    // the velocity keeps the direction and has the
                    // magnitude of 'absV'. To normalize the velocity,
                    // we use fact, that 'x' is the magnitude of the
                    // 'non-smooth' velocity (also, x != 0, because we
                    // already checked that in the if statement.)
                    filterVelocity_[phaseIdx] *= absV/x;

                    upstreamWeight_[phaseIdx] = (absV - vDn)/(vUp - vDn);;
                }
            }

            // normal velocity is the scalar product of the filter
            // velocity with the face normal
            filterVelocityNormal_[phaseIdx] = 0.0;
            for (int i = 0; i < Vector::size; ++i) 
                filterVelocityNormal_[phaseIdx] += filterVelocity_[phaseIdx][i]*normal[i];
        }
    }

    // local indices of the inside and the outside sub-control volumes
    short insideScvIdx_;
    short outsideScvIdx_;

    // extrusion factor for the sub-control volume face
    Scalar extrusionFactor_;

    // pressure potential gradients of all phases
    Vector potentialGrad_[numPhases];

    // filter velocities of all phases
    Vector filterVelocity_[numPhases];

    // normal velocities, i.e. filter velocity times face normal times
    // face area
    Scalar filterVelocityNormal_[numPhases];

    // the relative weight of the values of the upstream sub-control
    // volume compared to the ones downstream
    Scalar upstreamWeight_[numPhases];

    // data for the diffusion and the energy modules
    FluxVariablesDiffusion diffusionVars_;
    FluxVariablesEnergy energyVars_;
};

} // end namepace

#endif
