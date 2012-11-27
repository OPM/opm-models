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
 * \copydoc Ewoms::StokesVolumeVariables
 */
#ifndef EWOMS_STOKES_VOLUME_VARIABLES_HH
#define EWOMS_STOKES_VOLUME_VARIABLES_HH

#include "stokesproperties.hh"

#include <ewoms/disc/vcfv/vcfvvolumevariables.hh>
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup VCFVStokesModel
 * \ingroup VcfvVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the Stokes VCVF discretization.
 */
template <class TypeTag>
class StokesVolumeVariables
    : public VcfvVolumeVariables<TypeTag>
    , public VcfvEnergyVolumeVariables<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy) >
{
    typedef VcfvVolumeVariables<TypeTag> ParentType;
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
    typedef VcfvEnergyVolumeVariables<TypeTag, enableEnergy> EnergyVolumeVariables;

public:
    /*!
     * \copydoc VcfvVolumeVariables::update
     */
    void update(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvIdx, timeIdx);

        EnergyVolumeVariables::updateTemperatures_(fluidState_, elemCtx, scvIdx, timeIdx);

        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
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
            fluidState_.setMoleFraction(phaseIdx, compIdx, priVars[moleFrac1Idx + compIdx - 1]);
            sumx += priVars[moleFrac1Idx + compIdx - 1];
        }
        fluidState_.setMoleFraction(phaseIdx, 0, 1 - sumx);

        // create NullParameterCache and do dummy update
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        fluidState_.setDensity(phaseIdx,
                               FluidSystem::density(fluidState_,
                                                    paramCache,
                                                    phaseIdx));
        fluidState_.setViscosity(phaseIdx,
                                 FluidSystem::viscosity(fluidState_,
                                                        paramCache,
                                                        phaseIdx));

        // energy related quantities
        EnergyVolumeVariables::update_(fluidState_, paramCache, elemCtx, scvIdx, timeIdx);

        // the effective velocity of the control volume
        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            velocityCenter_[dimIdx] = priVars[Indices::velocity0Idx + dimIdx];

        // the gravitational acceleration applying to the material
        // inside the volume
        gravity_ = elemCtx.problem().gravity();
    }

    /*!
     * \copydoc VcfvVolumeVariables::updateScvGradients
     */
    void updateScvGradients(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        // calculate the pressure gradient at the SCV using finite
        // element gradients
        pressureGrad_ = 0.0;
        for (int i = 0; i < elemCtx.numScv(); ++i) {
            const auto &feGrad = elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx].gradCenter[i];
            DimVector tmp(feGrad);
            tmp *= elemCtx.volVars(i, timeIdx).fluidState().pressure(phaseIdx);

            pressureGrad_ += tmp;
        }

        // integrate the velocity over the sub-control volume
        //const auto &elemGeom = elemCtx.element().geometry();
        const auto &fvElemGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvLocalGeom = *fvElemGeom.subContVol[scvIdx].localGeometry;

        Dune::GeometryType geomType = scvLocalGeom.type();
        static const int quadratureOrder = 2;
        const auto &rule = Dune::QuadratureRules<Scalar,dimWorld>::rule(geomType, quadratureOrder);

        // integrate the veloc over the sub-control volume
        velocity_ = 0.0;
        for (auto it = rule.begin(); it != rule.end(); ++ it)
        {
            const auto &posScvLocal = it->position();
            const auto &posElemLocal = scvLocalGeom.global(posScvLocal);

            DimVector velocityAtPos = velocityAtPos_(elemCtx, timeIdx, posElemLocal);
            Scalar weight = it->weight();
            Scalar detjac = 1.0;
            //scvLocalGeom.integrationElement(posScvLocal) *
            //elemGeom.integrationElement(posElemLocal);
            velocity_.axpy(weight * detjac,  velocityAtPos);
        }

        // since we want the average velocity, we have to divide the
        // integrated value by the volume of the SCV
        //velocity_ /= fvElemGeom.subContVol[scvIdx].volume;
    }


    /*!
     * \brief Returns the thermodynamic state of the fluid for the control-volume.
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
    DimVector velocityAtPos_(const ElementContext elemCtx,
                          int timeIdx,
                          const LocalPosition &localPos) const
    {
        const auto &localFiniteElement =
            elemCtx.fvElemGeom(timeIdx).localFiniteElement();

        typedef Dune::FieldVector<Scalar, 1> ShapeValue;
        std::vector<ShapeValue> shapeValue;

        localFiniteElement.localBasis().evaluateFunction(localPos, shapeValue);

        DimVector result(0.0);
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++) {
            result.axpy(shapeValue[scvIdx][0], elemCtx.volVars(scvIdx, timeIdx).velocityCenter());
        }

        return result;
    }

    DimVector velocity_;
    DimVector velocityCenter_;
    DimVector gravity_;
    DimVector pressureGrad_;
    FluidState fluidState_;
};

}

#endif
