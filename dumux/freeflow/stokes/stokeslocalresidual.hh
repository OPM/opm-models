// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2012 by Katherina Baber                              *
 *   Copyright (C) 2009-2012 Klaus Mosthaf                                   *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief The local residual function for problems using the
 *        Stokes box model.
 */

#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_HH

#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"
#include "stokesproperties.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dune/grid/common/grid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxLocalResidual
 * \brief The local residual function for problems using the
 *        Stokes box model.
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        dimWorld = GridView::dimensionworld,
        phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex),
        numComponents = FluidSystem::numComponents
    };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

 public:
    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (mass and momentum) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();

        // mass storage
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] = 
                fs.molarity(phaseIdx, compIdx);

        // momentum balance
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            storage[momentum0EqIdx + axisIdx] =
                fs.density(phaseIdx) * volVars.velocity()[axisIdx];
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume. The face may be within
     *        an element (SCV face) or on the boundary. The advective and
     *        the diffusive fluxes are computed.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param faceIdx The index of the SCV face (may also be a boundary face)
     */
    void computeFlux(RateVector &flux, 
                     const ElementContext &elemCtx,
                     int faceIdx,
                     int timeIdx) const
    {
        flux = 0.0;
        asImp_().computeAdvectiveFlux(flux, elemCtx, faceIdx, timeIdx);
        Valgrind::CheckDefined(flux);
        asImp_().computeDiffusiveFlux(flux, elemCtx, faceIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective fluxes over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int faceIdx,
                              int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(faceIdx, timeIdx);

        // data attached to upstream vertex
        const VolumeVariables &up = elemCtx.volVars(fluxVars.upstreamIdx(), timeIdx);

        auto normal = fluxVars.normal();
        Scalar faceArea = normal.two_norm();
        normal /= faceArea;

        // mass fluxes
        Scalar vTimesN = fluxVars.velocity() * normal;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            flux[conti0EqIdx + compIdx] = 
                up.fluidState().molarity(phaseIdx, compIdx)
                * vTimesN;
        
        // momentum flux
        Scalar mu = up.fluidState().viscosity(phaseIdx);
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
        {
            // deal with the surface forces, i.e. the $\div[ \mu
            // (\grad[v] + \grad[v^T])]$ term on the right hand side
            // of the equation
            DimVector tmp;
            for (int j = 0; j < dimWorld; ++j) {
                tmp[j] = fluxVars.velocityGrad(/*velocityComp=*/axisIdx)[j];
                tmp[j] += fluxVars.velocityGrad(/*velocityComp=*/j)[axisIdx];
            }
            
            flux[momentum0EqIdx + axisIdx] = - mu * (tmp * normal);
            
            // this adds the convective momentum flux term $rho v
            // div[v]$ to the Stokes equation, transforming it to
            // Navier-Stokes.
            if (enableNavierTerm_()) {
                flux[momentum0EqIdx + axisIdx] +=
                    up.velocity()[axisIdx]
                    * (up.velocity() * normal);
            }
        }
        
        flux *= faceArea;
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        a SCV face or a boundary face.
     *
     * It doesn't do anything in the Stokes model but is used by the
     * transport and non-isothermal models to calculate diffusive and
     * conductive fluxes.
     *
     * \param flux The diffusive flux over the SCV face or boundary face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    { }

    /*!
     * \brief Calculate the source term of all equations.
     *        The pressure gradient at the center of a SCV is computed
     *        and the gravity term evaluated.
     *
     * \param q The source/sink in the sub control volume for each component
     * \param localVertexIdx The local index of the sub-control volume
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        assert(timeIdx == 0);
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);

        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);

        const auto &gravity = volVars.gravity();
        const auto &gradp = volVars.pressureGradient();
        Scalar density = volVars.fluidState().density(phaseIdx);

        Valgrind::CheckDefined(gravity);
        Valgrind::CheckDefined(gradp);
        Valgrind::CheckDefined(density);

        // deal with the pressure and volume terms
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            source[momentum0EqIdx + axisIdx] += gradp[axisIdx] - density*gravity[axisIdx];
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    static bool enableNavierTerm_()
    { return GET_PARAM(TypeTag, bool, EnableNavierTerm); }
};

}

#endif
