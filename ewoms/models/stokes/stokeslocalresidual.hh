// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser
  Copyright (C) 2012 by Christoph Grueninger

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::StokesLocalResidual
 */
#ifndef EWOMS_STOKES_LOCAL_RESIDUAL_HH
#define EWOMS_STOKES_LOCAL_RESIDUAL_HH

#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"
#include "stokesproperties.hh"

#include <opm/material/Valgrind.hpp>

#include <dune/grid/common/grid.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup StokesModel
 * \ingroup LocalResidual
 *
 * \brief The local residual function for problems using the
 *        Stokes model.
 */
template<class TypeTag>
class StokesLocalResidual
    : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dimWorld = GridView::dimensionworld };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { numComponents = FluidSystem::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \brief Register all run-time parameters for the local residual.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableNavierTerm,
                             "Enable the Navier term (convective flux term).");
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        const auto &volVars = elemCtx.volVars(dofIdx, timeIdx);
        const auto &fs = volVars.fluidState();

        // mass storage
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            storage[conti0EqIdx + compIdx] = fs.molarity(phaseIdx, compIdx);
        }
        Valgrind::CheckDefined(storage);

        // momentum balance
        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx) {
            storage[momentum0EqIdx + axisIdx] =
                fs.density(phaseIdx) * volVars.velocity()[axisIdx];
        }
        Valgrind::CheckDefined(storage);

        EnergyModule::addPhaseStorage(storage, elemCtx.volVars(dofIdx, timeIdx), phaseIdx);
        Valgrind::CheckDefined(storage);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux, const ElementContext &elemCtx,
                     int scvfIdx, int timeIdx) const
    {
        flux = 0.0;
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        // data attached to upstream DOF
        const VolumeVariables &up
            = elemCtx.volVars(fluxVars.upstreamIndex(phaseIdx), timeIdx);

        auto normal = fluxVars.normal();

        // mass fluxes
        Scalar vTimesN = fluxVars.velocity() * normal;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            flux[conti0EqIdx + compIdx]
                = up.fluidState().molarity(phaseIdx, compIdx) * vTimesN;

        // momentum flux
        Scalar mu = up.fluidState().viscosity(phaseIdx)
                    + fluxVars.eddyViscosity();
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // deal with the surface forces, i.e. the $\div[ \mu
            // (\grad[v] + \grad[v^T])]$ term on the right hand side
            // of the equation
            DimVector tmp;
            for (int j = 0; j < dimWorld; ++j) {
                tmp[j] = fluxVars.velocityGrad(/*velocityComp=*/axisIdx)[j];
                tmp[j] += fluxVars.velocityGrad(/*velocityComp=*/j)[axisIdx];
            }

            flux[momentum0EqIdx + axisIdx] = -mu * (tmp * normal);

            // this adds the convective momentum flux term $rho v
            // div[v]$ to the Stokes equation, transforming it to
            // Navier-Stokes.
            if (enableNavierTerm_()) {
                flux[momentum0EqIdx + axisIdx] += up.velocity()[axisIdx]
                                                  * (up.velocity() * normal);
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector &flux, const ElementContext &elemCtx,
                          int scvfIdx, int timeIdx) const
    {
        // heat conduction
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        assert(timeIdx == 0);
        const auto &volVars = elemCtx.volVars(dofIdx, timeIdx);

        // retrieve the source term intrinsic to the problem
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(source);

        const auto &gravity = volVars.gravity();
        const auto &gradp = volVars.pressureGradient();
        Scalar density = volVars.fluidState().density(phaseIdx);

        assert(std::isfinite(gradp.two_norm()));
        assert(std::isfinite(density));
        assert(std::isfinite(source.two_norm()));

        Valgrind::CheckDefined(gravity);
        Valgrind::CheckDefined(gradp);
        Valgrind::CheckDefined(density);

        // deal with the pressure and volume terms
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
            source[momentum0EqIdx + axisIdx] += gradp[axisIdx]
                                                - density * gravity[axisIdx];
    }

private:
    static bool enableNavierTerm_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableNavierTerm); }
};

} // namespace Ewoms

#endif
