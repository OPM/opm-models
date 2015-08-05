// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser
  Copyright (C) 2012 by Klaus Mosthaf

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
 * \copydoc Ewoms::StokesIntensiveQuantities
 */
#ifndef EWOMS_STOKES_INTENSIVE_QUANTITIES_HH
#define EWOMS_STOKES_INTENSIVE_QUANTITIES_HH

#include "stokesproperties.hh"

#include <ewoms/models/common/energymodule.hh>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <dune/geometry/quadraturerules.hh>

#include <dune/common/fvector.hh>

#include <vector>

namespace Ewoms {

/*!
 * \ingroup StokesModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the intensive quantities of the Stokes model.
 */
template <class TypeTag>
class StokesIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public EnergyIntensiveQuantities<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { numComponents = FluidSystem::numComponents };
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { pressureIdx = Indices::pressureIdx };
    enum { moleFrac1Idx = Indices::moleFrac1Idx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Ewoms::EnergyIntensiveQuantities<TypeTag, enableEnergy> EnergyIntensiveQuantities;

public:
    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        fluidState_.setPressure(phaseIdx, priVars[pressureIdx]);
        Valgrind::CheckDefined(fluidState_.pressure(phaseIdx));

        // set the saturation of the phase to 1. for the stokes model,
        // saturation is not a meaningful quanity, but it allows us to
        // reuse infrastructure written for the porous media models
        // more easily (e.g. the energy module)
        fluidState_.setSaturation(phaseIdx, 1.0);

        // set the phase composition
        Scalar sumx = 0;
        for (int compIdx = 1; compIdx < numComponents; ++compIdx) {
            fluidState_.setMoleFraction(phaseIdx, compIdx,
                                        priVars[moleFrac1Idx + compIdx - 1]);
            sumx += priVars[moleFrac1Idx + compIdx - 1];
        }
        fluidState_.setMoleFraction(phaseIdx, 0, 1 - sumx);

        // create NullParameterCache and do dummy update
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        fluidState_.setDensity(phaseIdx, FluidSystem::density(fluidState_, paramCache, phaseIdx));
        fluidState_.setViscosity(phaseIdx, FluidSystem::viscosity(fluidState_, paramCache, phaseIdx));

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // the effective velocity of the control volume
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            velocityCenter_[dimIdx] = priVars[Indices::velocity0Idx + dimIdx];

        // the gravitational acceleration applying to the material
        // inside the volume
        gravity_ = elemCtx.problem().gravity();
    }

    /*!
     * \copydoc IntensiveQuantities::updateScvGradients
     */
    void updateScvGradients(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    {
        // calculate the pressure gradient at the SCV using finite
        // element gradients
        pressureGrad_ = 0.0;
        for (int i = 0; i < elemCtx.numDof(/*timeIdx=*/0); ++i) {
            const auto &feGrad = elemCtx.stencil(timeIdx).subControlVolume(dofIdx).gradCenter[i];
            Valgrind::CheckDefined(feGrad);
            DimVector tmp(feGrad);
            tmp *= elemCtx.intensiveQuantities(i, timeIdx).fluidState().pressure(phaseIdx);
            Valgrind::CheckDefined(tmp);

            pressureGrad_ += tmp;
        }

        // integrate the velocity over the sub-control volume
        // const auto &elemGeom = elemCtx.element().geometry();
        const auto &stencil = elemCtx.stencil(timeIdx);
        const auto &scvLocalGeom = stencil.subControlVolume(dofIdx).localGeometry();

        Dune::GeometryType geomType = scvLocalGeom.type();
        static const int quadratureOrder = 2;
        const auto &rule = Dune::QuadratureRules<Scalar, dimWorld>::rule(geomType, quadratureOrder);

        // integrate the veloc over the sub-control volume
        velocity_ = 0.0;
        for (auto it = rule.begin(); it != rule.end(); ++it) {
            const auto &posScvLocal = it->position();
            const auto &posElemLocal = scvLocalGeom.global(posScvLocal);

            DimVector velocityAtPos = velocityAtPos_(elemCtx, timeIdx, posElemLocal);
            Scalar weight = it->weight();
            Scalar detjac = 1.0;
            // scvLocalGeom.integrationElement(posScvLocal) *
            // elemGeom.integrationElement(posElemLocal);
            velocity_.axpy(weight * detjac, velocityAtPos);
        }

        // since we want the average velocity, we have to divide the
        // integrated value by the volume of the SCV
        //velocity_ /= stencil.subControlVolume(dofIdx).volume;
    }

    /*!
     * \brief Returns the thermodynamic state of the fluid for the
     * control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the porosity of the medium
     *
     * For the Navier-Stokes model this quantity does not make sense
     * because there is no porous medium. The method is included here
     * to allow the Navier-Stokes model share the energy module with
     * the porous-media models.
     */
    Scalar porosity() const
    { return 1.0; }

    /*!
     * \brief Returns the average velocity in the sub-control volume.
     */
    const DimVector &velocity() const
    { return velocity_; }

    /*!
     * \brief Returns the velocity at the center in the sub-control volume.
     */
    const DimVector &velocityCenter() const
    { return velocityCenter_; }

    /*!
     * \brief Returns the pressure gradient in the sub-control volume.
     */
    const DimVector &pressureGradient() const
    { return pressureGrad_; }

    /*!
     * \brief Returns the gravitational acceleration vector in the
     *        sub-control volume.
     */
    const DimVector &gravity() const
    { return gravity_; }

private:
    DimVector velocityAtPos_(const ElementContext &elemCtx,
                             int timeIdx,
                             const LocalPosition &localPos) const
    {
        auto &feCache =
            elemCtx.gradientCalculator().localFiniteElementCache();
        const auto &localFiniteElement =
            feCache.get(elemCtx.element().type());

        typedef Dune::FieldVector<Scalar, 1> ShapeValue;
        std::vector<ShapeValue> shapeValue;

        localFiniteElement.localBasis().evaluateFunction(localPos, shapeValue);

        DimVector result(0.0);
        for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); dofIdx++) {
            result.axpy(shapeValue[dofIdx][0], elemCtx.intensiveQuantities(dofIdx, timeIdx).velocityCenter());
        }

        return result;
    }

    DimVector velocity_;
    DimVector velocityCenter_;
    DimVector gravity_;
    DimVector pressureGrad_;
    FluidState fluidState_;
};

} // namespace Ewoms

#endif
