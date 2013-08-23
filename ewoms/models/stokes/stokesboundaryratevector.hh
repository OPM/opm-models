// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Ewoms::StokesBoundaryRateVector
 */
#ifndef EWOMS_STOKES_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_STOKES_BOUNDARY_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/common/valgrind.hh>
#include <opm/material/constraintsolvers/ncpflash.hh>

#include "stokesvolumevariables.hh"

namespace Ewoms {

/*!
 * \ingroup StokesModel
 *
 * \brief Implements a boundary vector for the fully implicit (Navier-)Stokes model.
 */
template <class TypeTag>
class StokesBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numComponents = FluidSystem::numComponents };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef VcfvEnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    StokesBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    StokesBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector &)
     */
    StokesBoundaryRateVector(const StokesBoundaryRateVector &value)
        : ParentType(value)
    { }

    /*!
     * \param velocity The velocity vector [m/s] at the boundary.
     *
     * \param context The execution context for which the boundary rate should be specified.
     * \param bfIdx The local index of the boundary segment (-> local space index).
     * \param timeIdx The index used by the time discretization.
     * \param velocity The velocity vector [m/s] at the boundary.
     * \param fluidState The repesentation of the thermodynamic state of the system on the integration point of the boundary segment.
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context &context,
                     int bfIdx,
                     int timeIdx,
                     const DimVector &velocity,
                     const FluidState &fluidState)
    {
        const auto &fvElemGeom = context.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.boundaryFace[bfIdx];

        int insideScvIdx = context.insideScvIndex(bfIdx, timeIdx);
        //const auto &insideScv = fvElemGeom.subContVol[insideScvIdx];
        const auto &insideVolVars = context.volVars(bfIdx, timeIdx);

        // the outer unit normal
        auto normal = scvf.normal;
        normal /= normal.two_norm();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec
            = context.element().geometry().global(fvElemGeom.subContVol[insideScvIdx].localGeometry->center());
        const auto &scvPos = context.element().geometry().corner(insideScvIdx);
        distVec.axpy(-1, scvPos);
        Scalar dist = std::abs(distVec * normal);

        DimVector gradv[dimWorld];
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // Approximation of the pressure gradient at the boundary
            // segment's integration point.
            gradv[axisIdx] = normal;
            gradv[axisIdx] *= (velocity[axisIdx] - insideVolVars.velocity()[axisIdx])/dist;
            Valgrind::CheckDefined(gradv[axisIdx]);
        }

        // specify the mass fluxes over the boundary
        Scalar volumeFlux = velocity*normal;

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);
        Scalar density = FluidSystem::density(fluidState, paramCache, phaseIdx);
        Scalar molarDensity = density / fluidState.averageMolarMass(phaseIdx);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            (*this)[conti0EqIdx + compIdx] =
                volumeFlux
                * molarDensity
                * fluidState.moleFraction(phaseIdx, compIdx);
        }

        // calculate the momentum flux over the boundary
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // calculate a row of grad v + (grad v)^T
            DimVector tmp(0.0);
            for (int j = 0; j < dimWorld; ++j) {
                tmp[j] = gradv[axisIdx][j] + gradv[j][axisIdx];
            }


            // the momentum flux due to viscous forces
            (*this)[momentum0EqIdx + axisIdx] =
                - insideVolVars.fluidState().viscosity(phaseIdx)
                * (tmp * normal);
        }

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volumeFlux);
    }

    /*!
     * \brief Set a in-flow boundary in the (Navier-)Stoke model
     *
     * \param context The execution context for which the boundary rate should be specified.
     * \param bfIdx The local index of the boundary segment (-> local space index).
     * \param timeIdx The index used by the time discretization.
     * \param velocity The velocity vector [m/s] at the boundary.
     * \param fluidState The repesentation of the thermodynamic state of the system on the integration point of the boundary segment.
     */
    template <class Context, class FluidState>
    void setInFlow(const Context &context,
                   int bfIdx,
                   int timeIdx,
                   const DimVector &velocity,
                   const FluidState &fluidState)
    {
        const auto &volVars = context.volVars(bfIdx, timeIdx);

        setFreeFlow(context, bfIdx, timeIdx, velocity, fluidState);

        // don't let mass flow out
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] = std::min(0.0, (*this)[conti0EqIdx + compIdx]);

        // don't let momentum flow out
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            (*this)[momentum0EqIdx + axisIdx] = std::min(0.0, (*this)[momentum0EqIdx + axisIdx]);
    }

    /*!
     * \brief Set a out-flow boundary in the (Navier-)Stoke model
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    void setOutFlow(const Context &context,
                    int spaceIdx,
                    int timeIdx)
    {
        const auto &volVars = context.volVars(spaceIdx, timeIdx);

        DimVector velocity = volVars.velocity();
        const auto &fluidState = volVars.fluidState();

        setFreeFlow(context, spaceIdx, timeIdx, velocity, fluidState);

        // don't let mass flow in
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            (*this)[conti0EqIdx + compIdx] = std::max(0.0, (*this)[conti0EqIdx + compIdx]);

        // don't let momentum flow in
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            (*this)[momentum0EqIdx + axisIdx] = std::max(0.0, (*this)[momentum0EqIdx + axisIdx]);
    }

    /*!
     * \brief Set a no-flow boundary in the (Navier-)Stoke model
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    void setNoFlow(const Context &context,
                   int spaceIdx,
                   int timeIdx)
    {
        static DimVector v0(0.0);

        const auto &volVars = context.volVars(spaceIdx, timeIdx);
        const auto &fluidState = volVars.fluidState(); // don't care

        // no flow of mass and no slip for the momentum
        setFreeFlow(context, spaceIdx, timeIdx,
                    /*velocity = */v0,
                    fluidState);
    }
};

} // end namepace

#endif
