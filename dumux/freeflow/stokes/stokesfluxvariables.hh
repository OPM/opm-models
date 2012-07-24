// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Katherina Baber                              *
 *   Copyright (C) 2010-2011 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 *        the fluxes of the Stokes model over a face of a finite volume.
 *
 * This means pressure gradients, phase densities at the integration point, etc.
 */
#ifndef DUMUX_STOKES_FLUX_VARIABLES_HH
#define DUMUX_STOKES_FLUX_VARIABLES_HH

#include "stokesproperties.hh"

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux
{

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the mass and momentum fluxes over the face of a
 *        sub-control volume for the Stokes box model.
 *
 * This means pressure gradients, phase densities, viscosities, etc.
 * at the integration point of the sub-control-volume face
 */
template <class TypeTag>
class StokesFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx, bool isBoundaryFace = false)
    {
        const auto &fvGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = 
            isBoundaryFace?
            fvGeom.boundaryFace[scvfIdx]:
            fvGeom.subContVolFace[scvfIdx];

        onBoundary_ = isBoundaryFace;
        normal_ = scvf.normal;
        Valgrind::CheckDefined(normal_);

        // calculate gradients and secondary variables at IPs
        DimVector tmp(0.0);
        density_ = Scalar(0);
        molarDensity_ = Scalar(0);
        viscosity_ = Scalar(0);
        pressure_ = Scalar(0);
        normalVelocity_ = Scalar(0);
        velocity_ = Scalar(0);
        pressureGrad_ = Scalar(0);
        
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            velocityGrad_[dimIdx] = 0.0;

        for (int idx = 0; idx < elemCtx.numScv(); idx++)
        {
            const auto &volVars = elemCtx.volVars(idx, timeIdx);
            const auto &fs = volVars.fluidState();

            // phase density and viscosity at IP
            density_ += 
                fs.density(phaseIdx)
                * scvf.shapeValue[idx];
            molarDensity_ += 
                fs.molarDensity(phaseIdx)
                * scvf.shapeValue[idx];
            viscosity_ +=
                fs.viscosity(phaseIdx) 
                * scvf.shapeValue[idx];
            pressure_ +=
                fs.pressure(phaseIdx)
                * scvf.shapeValue[idx];

            // velocity at the IP (fluxes)
            DimVector velocityTimesShapeValue = volVars.velocityCenter();
            velocityTimesShapeValue *= scvf.shapeValue[idx];
            Valgrind::CheckDefined(scvf.shapeValue[idx]);
            velocity_ += velocityTimesShapeValue;

            // the pressure gradient
            tmp = scvf.grad[idx];
            tmp *= fs.pressure(phaseIdx);
            pressureGrad_ += tmp;
            // take gravity into account
            tmp = elemCtx.problem().gravity(elemCtx, idx, timeIdx);
            tmp *= density_;
            // pressure gradient including influence of gravity
            pressureGrad_ -= tmp;

            // the velocity gradients
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            {
                tmp = scvf.grad[idx];
                tmp *= volVars.velocityCenter()[dimIdx];
                velocityGrad_[dimIdx] += tmp;
            }
        }
        Valgrind::CheckDefined(velocity_);

        normalVelocity_ = velocity_ * normal_;
        Valgrind::CheckDefined(normalVelocity_);

        // set the upstream and downstream vertices
        upstreamIdx_ = scvf.i;
        downstreamIdx_ = scvf.j;
        if (normalVelocity_ < 0)
            std::swap(upstreamIdx_, downstreamIdx_);

        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(viscosity_);
        Valgrind::CheckDefined(velocity_);
        Valgrind::CheckDefined(pressureGrad_);
        Valgrind::CheckDefined(velocityGrad_);
    }

public:
    /*!
     * \brief Return the pressure \f$\mathrm{[Pa]}\f$ at the integration
     *        point.
     */
    Scalar pressure() const
    { return pressure_; }

    /*!
     * \brief Return the mass density \f$ \mathrm{[kg/m^3]} \f$ at the integration
     *        point.
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Return the molar density \f$ \mathrm{[mol/m^3]} \f$ at the integration point.
     */
    const Scalar molarDensity() const
    { return molarDensity_; }

    /*!
     * \brief Return the viscosity \f$ \mathrm{[m^2/s]} \f$ at the integration
     *        point.
     */
    Scalar viscosity() const
    { return viscosity_; }

    /*!
     * \brief Return the velocity \f$ \mathrm{[m/s]} \f$ at the integration
     *        point multiplied by the normal and the area.
     */
    Scalar normalVelocity() const
    { return normalVelocity_; }

    /*!
     * \brief Return the pressure gradient at the integration point.
     */
    const DimVector &pressureGrad() const
    { return pressureGrad_; }

    /*!
     * \brief Return the velocity vector at the integration point.
     */
    const DimVector &velocity() const
    { return velocity_; }

    /*!
     * \brief Return the velocity gradient at the integration
     *        point of a face.
     */
    const DimVector &velocityGrad(int axisIdx) const
    { return velocityGrad_[axisIdx]; }

    /*!
     * \brief Return the eddy viscosity (if implemented).
     */
    Scalar eddyViscosity() const
    { return 0; }

     /*!
     * \brief Return the eddy diffusivity (if implemented).
     */
    Scalar eddyDiffusivity() const
    { return 0; }

    /*!
     * \brief Return the local index of the upstream sub-control volume.
     */
    int upstreamIdx() const
    { return upstreamIdx_; }

    /*!
     * \brief Return the local index of the downstream sub-control volume.
     */
    int downstreamIdx() const
    { return downstreamIdx_; }

    /*!
     * \brief Indicates if a face is on a boundary. Used for in the
     *        face() method (e.g. for outflow boundary conditions).
     */
    bool onBoundary() const
    { return onBoundary_; }
    
    /*!
     * \brief Returns the extrusionFactor of the face.
     */
    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \brief Returns normal vector of the face of the flux variables.
     */
    const DimVector &normal() const
    { return normal_; }

protected:
    bool onBoundary_;

    // values at the integration point
    Scalar density_;
    Scalar molarDensity_;
    Scalar viscosity_;
    Scalar pressure_;
    Scalar normalVelocity_;
    DimVector velocity_;
    DimVector normal_;

    // gradients at the IPs
    DimVector pressureGrad_;
    DimVector velocityGrad_[dimWorld];

    // local index of the upwind vertex
    int upstreamIdx_;
    // local index of the downwind vertex
    int downstreamIdx_;
};

} // end namepace

#endif
