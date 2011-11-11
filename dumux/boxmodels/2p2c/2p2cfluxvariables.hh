// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief   This file contains the data which is required to calculate
 *          all fluxes of components over a face of a finite volume for
 *          the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_FLUX_VARIABLES_HH
#define DUMUX_2P2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "2p2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents =  GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> Vector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

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
     * \brief Returns th extrusion factor for the sub-control volume face
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase at the
     *        face's integration point [Pa/m]
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
    { return 0.0; }

    Scalar porousDiffCoeff(int phaseIdx, int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return porousDiffCoeff_[phaseIdx];
    };

    Scalar molarDensity(int phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return molarDensity_[phaseIdx];
    };

    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return moleFrac_[phaseIdx][compIdx];
    };

    const Vector &moleFracGrad(int phaseIdx, int compIdx) const
    {
        return moleFracGrad_[phaseIdx][compIdx];
    };

    const Vector &temperatureGrad() const
    { return temperatureGrad_; };

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx,
                             int timeIdx)
    {
        // reset all quantities to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            molarDensity_[phaseIdx] = Scalar(0);
            porousDiffCoeff_[phaseIdx] = Scalar(0);
            potentialGrad_[phaseIdx] = Scalar(0);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFracGrad_[phaseIdx][compIdx] = Scalar(0);
                moleFrac_[phaseIdx][compIdx] = Scalar(0);
            }
        }
        temperatureGrad_ = Scalar(0);

        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];

        // calculate gradients
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx)
        {
            // FE gradient at vertex
            const Vector &feGrad = scvf.grad[scvIdx];
            Scalar shapeValue = scvf.shapeValue[scvIdx];
            const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
            const auto &fluidState = volVars.fluidState();
            Vector tmp;

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= fluidState.pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;

                molarDensity_[phaseIdx] += shapeValue * fluidState.molarDensity(phaseIdx);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    tmp = feGrad;
                    tmp *= fluidState.moleFraction(phaseIdx, compIdx);
                    moleFracGrad_[phaseIdx][compIdx] += tmp;

                    moleFrac_[phaseIdx][compIdx] +=
                        shapeValue * fluidState.moleFraction(phaseIdx, compIdx);
                }
            }

            tmp = feGrad;
            tmp *= fluidState.temperature(/*phaseIdx=*/0);
            temperatureGrad_ += tmp;
        }


        const auto &volVarsI = elemCtx.volVars(insideScvIdx_, timeIdx);
        const auto &volVarsJ = elemCtx.volVars(outsideScvIdx_, timeIdx);
        const auto &fsI = volVarsI.fluidState();
        const auto &fsJ = volVarsJ.fluidState();

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient
            Scalar tauI =
                1.0/(volVarsI.porosity() * volVarsI.porosity())
                *
                pow(volVarsI.porosity() * fsI.saturation(phaseIdx),
                    7.0/3);
            Scalar tauJ =
                1.0/(volVarsJ.porosity() * volVarsJ.porosity())
                *
                pow(volVarsJ.porosity() * fsJ.saturation(phaseIdx),
                    7.0/3);

            // Diffusion coefficient in the porous medium
            // -> harmonic mean
            porousDiffCoeff_[phaseIdx] =
                harmonicMean(volVarsI.porosity()
                             * fsI.saturation(phaseIdx)
                             * tauI
                             * volVarsI.diffCoeff(phaseIdx)
                             ,
                             volVarsJ.porosity()
                             * fsJ.saturation(phaseIdx)
                             * tauJ
                             * volVarsJ.diffCoeff(phaseIdx));
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

            const auto &fsIn = volVarsI.fluidState();
            const auto &fsOut = volVarsJ.fluidState();
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
                              int scvfIdx,
                              int timeIdx)
    {
        const SpatialParameters &spatialParams =
            elemCtx.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                insideScvIdx_,
                                                                timeIdx),
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                outsideScvIdx_,
                                                                timeIdx));

        const Vector &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the "prelimanary" filter velocity not
            // taking the mobility into account
            K.mv(potentialGrad_[phaseIdx], filterVelocity_[phaseIdx]);
            filterVelocity_[phaseIdx] *= -1;

            // normal velocity is the scalar product of the filter
            // velocity with the face normal
            filterVelocityNormal_[phaseIdx] = 0.0;
            for (int i = 0; i < Vector::size; ++i)
                filterVelocityNormal_[phaseIdx] += filterVelocity_[phaseIdx][i]*normal[i];

            // multiply both with the upstream mobility
            const auto &up = elemCtx.volVars(upstreamIdx(phaseIdx), timeIdx);
            filterVelocityNormal_[phaseIdx] *= up.mobility(phaseIdx);
            filterVelocity_[phaseIdx] *= up.mobility(phaseIdx);
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

    Scalar porousDiffCoeff_[numPhases];
    Scalar molarDensity_[numPhases];
    Scalar moleFrac_[numPhases][numComponents];
    Vector moleFracGrad_[numPhases][numComponents];
    Vector temperatureGrad_;
};

} // end namepace

#endif
