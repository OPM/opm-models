// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::StokesLocalResidual
 */
#ifndef EWOMS_STOKES_LOCAL_RESIDUAL_HH
#define EWOMS_STOKES_LOCAL_RESIDUAL_HH

#include "stokesintensivequantities.hh"
#include "stokesextensivequantities.hh"
#include "stokesproperties.hh"

#include <opm/material/common/MathToolbox.hpp>

#include <opm/common/Valgrind.hpp>

#include <dune/grid/common/grid.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup StokesModel
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
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dimWorld = GridView::dimensionworld };
    enum { stokesPhaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numComponents = FluidSystem::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Opm::MathToolbox<Evaluation> Toolbox;

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
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        if (elemCtx.focusDofIndex() == dofIdx)
            computeStorage_<LhsEval>(storage, elemCtx, dofIdx, timeIdx);
        else
            computeStorage_<Scalar>(storage, elemCtx, dofIdx, timeIdx);
    }

    template <class DensityEval, class LhsEval>
    void computeStorage_(Dune::FieldVector<LhsEval, numEq>& storage,
                         const ElementContext& elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx) const
    {
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        storage = 0.0;

        // mass
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] = Opm::decay<DensityEval>(fs.molarity(stokesPhaseIdx, compIdx));
        Opm::Valgrind::CheckDefined(storage);

        // momentum
        for (unsigned axisIdx = 0; axisIdx < dimWorld; ++ axisIdx) {
            storage[momentum0EqIdx + axisIdx] =
                Opm::decay<DensityEval>(fs.density(stokesPhaseIdx)) *
                Opm::decay<LhsEval>(intQuants.velocity()[axisIdx]);
        }
        Opm::Valgrind::CheckDefined(storage);

        // energy
        EnergyModule::addPhaseStorage(storage, elemCtx.intensiveQuantities(dofIdx, timeIdx), stokesPhaseIdx);
        Opm::Valgrind::CheckDefined(storage);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        flux = 0.0;
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Opm::Valgrind::CheckDefined(flux);
        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Opm::Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        auto upIdx = extQuants.upstreamIndex(stokesPhaseIdx);
        auto focusDofIdx = elemCtx.focusDofIndex();

        // data attached to upstream DOF
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        auto normal = extQuants.normal();

        // mass fluxes
        Evaluation vTimesN = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            vTimesN += extQuants.velocity()[dimIdx]*normal[dimIdx];

        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            flux[conti0EqIdx + compIdx] = up.fluidState().molarity(stokesPhaseIdx, compIdx) * vTimesN;

        // momentum flux
        Evaluation mu;
        if (upIdx == focusDofIdx)
            mu = up.fluidState().viscosity(stokesPhaseIdx) + extQuants.eddyViscosity();
        else
            mu =
                Opm::scalarValue(up.fluidState().viscosity(stokesPhaseIdx)) +
                Opm::scalarValue(extQuants.eddyViscosity());
        for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx) {
            // deal with the surface forces, i.e. the $\div[ \mu
            // (\grad[v] + \grad[v^T])]$ term on the right hand side
            // of the equation
            EvalDimVector tmp;
            for (unsigned j = 0; j < dimWorld; ++j) {
                tmp[j] = extQuants.velocityGrad(/*velocityComp=*/axisIdx)[j];
                tmp[j] += extQuants.velocityGrad(/*velocityComp=*/j)[axisIdx];
            }

            Evaluation alpha = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                alpha += tmp[dimIdx] * normal[dimIdx];
            flux[momentum0EqIdx + axisIdx] = -mu*alpha;

            // this adds the convective momentum flux term $rho v
            // div[v]$ to the Stokes equation, transforming it to
            // Navier-Stokes.
            if (enableNavierTerm_()) {
                for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    flux[momentum0EqIdx + axisIdx] +=
                        up.velocity()[axisIdx] * up.velocity()[dimIdx]*normal[dimIdx];
            }
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector& flux,
                          const ElementContext& elemCtx,
                          unsigned scvfIdx,
                          unsigned timeIdx) const
    {
        // heat conduction
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        assert(timeIdx == 0);
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);

        // retrieve the source term intrinsic to the problem
        Opm::Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(source);

        // if the current DOF is not the one we currently focus on, we need
        // to throw away the derivatives!
        if (!std::is_same<Scalar, Evaluation>::value
            && dofIdx != elemCtx.focusDofIndex())
        {
            for (int i = 0; i < numEq; ++i)
                source[i] = Toolbox::value(source[i]);
        }

        const auto& gravity = intQuants.gravity();
        const auto& gradp = intQuants.pressureGradient();
        const auto& fs = intQuants.fluidState();

#ifndef NDEBUG
        assert(Opm::isfinite(fs.density(stokesPhaseIdx)));
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++ dimIdx) {
            assert(Opm::isfinite(gradp[dimIdx]));
            assert(Opm::isfinite(source[dimIdx]));
        }

        Opm::Valgrind::CheckDefined(gravity);
        Opm::Valgrind::CheckDefined(gradp);
        Opm::Valgrind::CheckDefined(fs.density(stokesPhaseIdx));
#endif

        auto focusDofIdx = elemCtx.focusDofIndex();

        // deal with the pressure and volumetric terms
        if (dofIdx == focusDofIdx) {
            for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                source[momentum0EqIdx + axisIdx] +=
                    gradp[axisIdx] -
                    fs.density(stokesPhaseIdx)*gravity[axisIdx];
        }
        else {
            for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                source[momentum0EqIdx + axisIdx] +=
                    gradp[axisIdx] -
                    Opm::scalarValue(fs.density(stokesPhaseIdx))*gravity[axisIdx];
        }
    }

private:
    static bool enableNavierTerm_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableNavierTerm); }
};

} // namespace Ewoms

#endif
