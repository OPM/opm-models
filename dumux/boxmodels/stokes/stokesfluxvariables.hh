// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
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
 * \copydoc Dumux::StokesFluxVariables
 */
#ifndef DUMUX_STOKES_FLUX_VARIABLES_HH
#define DUMUX_STOKES_FLUX_VARIABLES_HH

#include "stokesproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/math.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxFluxVariables
 *
 * \brief Contains the data which is required to calculate the mass
 *        and momentum fluxes over the face of a sub-control volume
 *        for the Stokes box model.
 *
 * This means pressure gradients, phase densities, viscosities, etc.
 * at the integration point of the sub-control-volume face
 */
template <class TypeTag>
class StokesFluxVariables
    : public BoxMultiPhaseEnergyFluxVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dimWorld = GridView::dimensionworld };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef BoxMultiPhaseEnergyFluxVariables<TypeTag, enableEnergy> EnergyFluxVariables;

public:
    /*!
     * \brief Register all run-time parameters for the flux variables.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Update all quantities which are required on an
     *        intersection between two finite volumes.
     *
     * \param elemCtx The current execution context.
     * \param scvfIdx The local index of the sub-control volume face.
     * \param timeIdx The index relevant for the time discretization.
     * \param isBoundaryFace Specifies whether the sub-control-volume
     *                       face is on the domain boundary or not.
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx, bool isBoundaryFace = false)
    {
        const auto &fvGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = 
            isBoundaryFace?
            fvGeom.boundaryFace[scvfIdx]:
            fvGeom.subContVolFace[scvfIdx];

        insideIdx_ = scvf.i;
        outsideIdx_ = scvf.j;

        onBoundary_ = isBoundaryFace;
        normal_ = scvf.normal;
        Valgrind::CheckDefined(normal_);

        // calculate gradients and secondary variables at IPs
        DimVector tmp(0.0);
        density_ = Scalar(0);
        molarDensity_ = Scalar(0);
        viscosity_ = Scalar(0);
        pressure_ = Scalar(0);
        volumeFlux_ = Scalar(0);
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

        volumeFlux_ = velocity_ * normal_;
        Valgrind::CheckDefined(volumeFlux_);

        // set the upstream and downstream vertices
        upstreamIdx_ = scvf.i;
        downstreamIdx_ = scvf.j;
        if (volumeFlux_ < 0)
            std::swap(upstreamIdx_, downstreamIdx_);

        EnergyFluxVariables::update_(elemCtx, scvfIdx, timeIdx);

        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(viscosity_);
        Valgrind::CheckDefined(velocity_);
        Valgrind::CheckDefined(pressureGrad_);
        Valgrind::CheckDefined(velocityGrad_);
    }

    /*!
     * \copydoc BoxMultiPhaseFluxVariables::updateBoundary
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context, 
                        int bfIdx, 
                        int timeIdx, 
                        const FluidState &fs, 
                        typename FluidSystem::ParameterCache &paramCache)
    {      
        update(context, bfIdx, timeIdx, fs, paramCache, /*isOnBoundary=*/true);
    }

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
     * \brief Return the volume flux of mass
     */
    Scalar volumeFlux(int phaseIdx) const
    { return volumeFlux_; }

    /*!
     * \brief Return the weight of the upstream index
     */
    Scalar upstreamWeight(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream index
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 0.0; }

    /*!
     * \brief Return the local index of the upstream sub-control volume.
     */
    int upstreamIndex(int phaseIdx) const
    { return upstreamIdx_; }

    /*!
     * \brief Return the local index of the downstream sub-control volume.
     */
    int downstreamIndex(int phaseIdx) const
    { return downstreamIdx_; }

    /*!
     * \brief Return the local index of the sub-control volume which is located in negative normal direction.
     */
    int insideIndex() const
    { return insideIdx_; }

    /*!
     * \brief Return the local index of the sub-control volume which is located in negative normal direction.
     */
    int outsideIndex() const
    { return outsideIdx_; }

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

private:
    bool onBoundary_;

    // values at the integration point
    Scalar density_;
    Scalar molarDensity_;
    Scalar viscosity_;
    Scalar pressure_;
    Scalar volumeFlux_;
    DimVector velocity_;
    DimVector normal_;

    // gradients at the IPs
    DimVector pressureGrad_;
    DimVector velocityGrad_[dimWorld];

    // local index of the upwind vertex
    int upstreamIdx_;
    // local index of the downwind vertex
    int downstreamIdx_;

    int insideIdx_;
    int outsideIdx_;
};

} // end namespace

#endif
