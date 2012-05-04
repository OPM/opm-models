// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief Implements a boundary vector for the fully implicit Stokes model.
 */
#ifndef DUMUX_BOX_STOKES_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_BOX_STOKES_BOUNDARY_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>

#include "stokesvolumevariables.hh"

namespace Dumux
{
/*!
 * \ingroup StokesModel
 *
 * \brief Implements a boundary vector for the fully implicit Stokes model.
 */
template <class TypeTag>
class StokesBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numComponents = FluidSystem::numComponents };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*!
     * \brief Default constructor
     */
    StokesBoundaryRateVector()
        : ParentType()
    { };

    /*!
     * \brief Constructor with assignment from scalar
     */
    StokesBoundaryRateVector(Scalar value)
        : ParentType(value)
    { };

    /*!
     * \brief Copy constructor
     */
    StokesBoundaryRateVector(const StokesBoundaryRateVector &value)
        : ParentType(value)
    { };

    /*!
     * \brief Specify a free-flow boundary
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context &context,
                     int spaceIdx,
                     int timeIdx,
                     const Vector &velocity, 
                     const FluidState &fs) 
    {
        const auto &fvElemGeom = context.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.boundaryFace[spaceIdx];

        int insideScvIdx = context.insideScvIndex(spaceIdx, timeIdx);
        const auto &insideScv = fvElemGeom.subContVol[insideScvIdx];
        const auto &insideVolVars = context.volVars(spaceIdx, timeIdx);

        // the outer unit normal
        auto normal = context.fvElemGeom(timeIdx).boundaryFace[spaceIdx].normal;
        normal /= normal.two_norm();
      
        // distance between the center of the SCV and center of the boundary face
        Vector distVec = scvf.ipGlobal;
        const auto &scvCenter = context.element().geometry().global(insideScv.localCenter);
        for (int i = 0; i < dimWorld; ++ i)
            distVec[i] -= scvCenter[i];
        Scalar dist = distVec.two_norm();

        Vector gradv[dimWorld];
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // Approximation of the pressure gradient at the boundary
            // segment's integration point.
            gradv[axisIdx] = normal;
            gradv[axisIdx] *= (velocity[axisIdx] - insideVolVars.velocity()[axisIdx])/dist;
        }

        // specify the mass fluxes over the boundary
        Scalar vTimesN = velocity*normal;

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fs, phaseIdx);
        Scalar molarDensity = 
            FluidSystem::density(fs, paramCache, phaseIdx)
            / fs.averageMolarMass(phaseIdx);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            (*this)[conti0EqIdx + compIdx] = 
                vTimesN
                * molarDensity
                * fs.moleFraction(phaseIdx, compIdx);
        }
        
        // calculate the momentum flux over the boundary
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // calculate a row of grad v + (grad v)^T
            Vector tmp(0.0);
            for (int j = 0; j < dimWorld; ++j) {
                tmp[j] += gradv[axisIdx][j];
                tmp[j] += gradv[j][axisIdx];
            }

            // the momentum flux due to viscous forces
            (*this)[momentum0EqIdx + axisIdx] = 
                insideVolVars.fluidState().viscosity(phaseIdx) 
                * (tmp * normal);
        }
    }

    /*!
     * \brief Specify an inflow boundary
     */
    template <class Context, class FluidState>
    void setInFlow(const Context &context,
                   int spaceIdx,
                   int timeIdx,
                   const Vector &velocity, 
                   const FluidState &fs) 
    {
        const auto &volVars = context.volVars(spaceIdx, timeIdx);

        setFreeFlow(context, spaceIdx, timeIdx, velocity, fs);

        // don't let mass flow out
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] = std::min(0.0, (*this)[conti0EqIdx + compIdx]);

        // don't let momentum flow out
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            (*this)[momentum0EqIdx + axisIdx] = std::min(0.0, (*this)[momentum0EqIdx + axisIdx]);
    }

    /*!
     * \brief Specify an outflow boundary
     */
    template <class Context>
    void setOutFlow(const Context &context,
                    int spaceIdx,
                    int timeIdx)
    {
        const auto &volVars = context.volVars(spaceIdx, timeIdx);

        Vector velocity = volVars.velocity();
        const auto &fs = volVars.fluidState();

        setFreeFlow(context, spaceIdx, timeIdx, velocity, fs);

        // don't let mass flow in
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] = std::max(0.0, (*this)[conti0EqIdx + compIdx]);
        
        // don't let momentum flow in
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            (*this)[momentum0EqIdx + axisIdx] = std::max(0.0, (*this)[momentum0EqIdx + axisIdx]);
    }
    
    /*!
     * \brief Specify a no-flow boundary.
     */
    template <class Context>
    void setNoFlow(const Context &context,
                   int spaceIdx,
                   int timeIdx)
    { 
        static Vector v0(0.0);

        const auto &volVars = context.volVars(spaceIdx, timeIdx);
        const auto &fs = volVars.fluidState(); // don't care

        // no flow of mass and no slip for the momentum
        setFreeFlow(context, spaceIdx, timeIdx, 
                    /*velocity = */v0,
                    fs);
    }

protected:
    Implementation &asImp_() 
    { return *static_cast<Implementation *>(this); }
/*
    template <class FluidState>
    void enthalpyFlux_(const FluxVariables &fluxVars,
                       const VolumeVariables &insideVolVars,
                       const FluidState &fs,
                       const typename FluidSystem::ParameterCache &paramCache,
                       int phaseIdx,
                       Scalar density)
    { }
*/
};

} // end namepace

#endif
