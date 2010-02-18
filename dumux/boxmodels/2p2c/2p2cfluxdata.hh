// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \ingroup TwoPTwoCBoxModel
 * \brief This file contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities, etc. at the integration points of the control volume
 */
#ifndef DUMUX_2P2C_FLUX_DATA_HH
#define DUMUX_2P2C_FLUX_DATA_HH

#include <dumux/auxiliary/math.hh>

namespace Dune
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoPTwoCFluxData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))    Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef std::vector<VertexData>                      VertexDataArray;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume             SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace         SCVFace;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
    };

public:
    TwoPTwoCFluxData(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int faceIdx,
                     const VertexDataArray &elemDat)
        : fvElemGeom(elemGeom)
    {
        scvfIdx = faceIdx;
        
        densityAtIP = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            potentialGrad[phase] = Scalar(0);
            concentrationGrad[phase] = Scalar(0);
        }

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
    };

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const VertexDataArray &elemDat)
    {
        // calculate gradients
        GlobalPosition tmp(0.0);
        for (int idx = 0;
             idx < fvElemGeom.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const LocalPosition &feGrad = face().grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure(phase);
                potentialGrad[phase] += tmp;

                // phase density
                densityAtIP[phase]
                    +=
                    elemDat[idx].density(phase) *
                    face().shapeValue[idx];
            }

            // the concentration gradient of the non-wetting
            // component in the wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].phaseState().massFrac(wPhaseIdx, nCompIdx);
            concentrationGrad[wPhaseIdx] += tmp;

            // the concentration gradient of the wetting component
            // in the non-wetting phase
            tmp = feGrad;
            tmp *= elemDat[idx].phaseState().massFrac(nPhaseIdx, wCompIdx);
            concentrationGrad[nPhaseIdx] += tmp;
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
        {
            tmp = problem.gravity();
            tmp *= densityAtIP[phase];

            potentialGrad[phase] -= tmp;
        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const VertexDataArray &elemDat)
    {
        // multiply the pressure potential with the intrinsic
        // permeability
        GlobalPosition Kmvp;
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            problem.spatialParameters().Kmv(Kmvp,
                                            potentialGrad[phaseIdx],
                                            element,
                                            fvElemGeom,
                                            scvfIdx);
            KmvpNormal[phaseIdx] = Kmvp * face().normal;
        }

        // set the upstream and downstream vertices
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            upstreamIdx[phaseIdx] = face().i;
            downstreamIdx[phaseIdx] = face().j;

            if (KmvpNormal[phaseIdx] < 0) {
                std::swap(upstreamIdx[phaseIdx],
                          downstreamIdx[phaseIdx]);
            }
        }
    }

    void calculateDiffCoeffPM_(const Problem &problem,
                               const Element &element,
                               const VertexDataArray &elemDat)
    {
        const VertexData &vDat_i = elemDat[face().i];
        const VertexData &vDat_j = elemDat[face().j];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to only calculate diffusion coefficents
            // for phases which exist in both finite volumes
            if (vDat_i.saturation(phaseIdx) <= 0 ||
                vDat_j.saturation(phaseIdx) <= 0)
            {
                diffCoeffPM[phaseIdx] = 0.0;
                continue;
            }

            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient

            Scalar tau_i =
                1.0/(vDat_i.porosity() * vDat_i.porosity()) *
                pow(vDat_i.porosity() * vDat_i.saturation(phaseIdx), 7.0/3);
            Scalar tau_j =
                1.0/(vDat_j.porosity() * vDat_j.porosity()) *
                pow(vDat_j.porosity() * vDat_j.saturation(phaseIdx), 7.0/3);
            // Diffusion coefficient in the porous medium

            // -> arithmetic mean
            diffCoeffPM[phaseIdx]
                = 1./2*(vDat_i.porosity() * vDat_i.saturation(phaseIdx) * tau_i * vDat_i.diffCoeff(phaseIdx) +
                        vDat_j.porosity() * vDat_j.saturation(phaseIdx) * tau_j * vDat_j.diffCoeff(phaseIdx));
            // -> harmonic mean
            // = harmonicMean_(vDat_i.porosity * vDat_i.saturation[phaseIdx] * tau_i * vDat_i.diffCoeff[phaseIdx],
            //                 vDat_j.porosity * vDat_j.saturation[phaseIdx] * tau_j * vDat_j.diffCoeff[phaseIdx]);

        }
    }

public:
    const SCVFace &face() const
    { return fvElemGeom.subContVolFace[scvfIdx]; }

    const FVElementGeometry &fvElemGeom;
    int                      scvfIdx;

    // gradients
    GlobalPosition potentialGrad[numPhases];
    GlobalPosition concentrationGrad[numPhases];

    // density of each face at the integration point
    PhasesVector densityAtIP;

    // darcy velocities of each phase (without the mobility)
    GlobalPosition vDarcy[numPhases];

    // intrinsic permeability times pressure potential gradient
    // projected on the face normal
    Scalar KmvpNormal[numPhases];

    // local index of the upwind vertex for each phase
    int upstreamIdx[numPhases];
    // local index of the downwind vertex for each phase
    int downstreamIdx[numPhases];

    // the diffusion coefficient for the porous medium
    PhasesVector diffCoeffPM;
};

} // end namepace

#endif
