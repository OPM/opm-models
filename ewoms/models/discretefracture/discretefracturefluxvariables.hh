// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2010 by Melanie Darcis                               *
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
 * \copydoc Ewoms::DiscreteFractureFluxVariables
 */
#ifndef EWOMS_DISCRETE_FRACTURE_FLUX_VARIABLES_HH
#define EWOMS_DISCRETE_FRACTURE_FLUX_VARIABLES_HH

#include <ewoms/models/immiscible/immisciblefluxvariables.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup DiscreteFractureVcfvModel
 * \ingroup VCFVFluxVariables
 *
 * \brief This class provides the data all quantities that are required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the discrete-fracture immiscible multi-phase model.
 */
template <class TypeTag>
class DiscreteFractureFluxVariables
    : public ImmiscibleFluxVariables<TypeTag>
{
    typedef ImmiscibleFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = FluidSystem::numPhases };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \copydoc VcfvMultiPhaseFluxVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);

        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &fvElemGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.subContVolFace[scvfIdx];
        int insideScvIdx = scvf.i;
        int outsideScvIdx = scvf.j;

        int globalI = elemCtx.globalSpaceIndex(insideScvIdx, timeIdx);
        int globalJ = elemCtx.globalSpaceIndex(outsideScvIdx, timeIdx);
        const auto &fractureMapper = elemCtx.problem().fractureMapper();
        if (!fractureMapper.isFractureEdge(globalI, globalJ))
            // do nothing if no fracture goes though the current edge
            return;

        // average the intrinsic permeability of the fracture
        const auto &volVarsInside = elemCtx.volVars(insideScvIdx, timeIdx);
        const auto &volVarsOutside = elemCtx.volVars(outsideScvIdx, timeIdx);
        elemCtx.problem().meanK(fractureIntrinsicPermeability_,
                                volVarsInside.fractureIntrinsicPermeability(),
                                volVarsOutside.fractureIntrinsicPermeability());

        auto distDirection = elemCtx.pos(outsideScvIdx, timeIdx);
        distDirection -= elemCtx.pos(insideScvIdx, timeIdx);
        distDirection /= distDirection.two_norm();

        const auto &problem = elemCtx.problem();
        fractureWidth_ = problem.fractureWidth(elemCtx, insideScvIdx, outsideScvIdx, timeIdx);
        assert(fractureWidth_ < scvf.normal.two_norm());

        const auto &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto &pGrad = fluxVars.potentialGrad(phaseIdx);

            int upstreamIdx = evalPointFluxVars.upstreamIndex(phaseIdx);
            const auto &up = elemCtx.volVars(upstreamIdx, timeIdx);

            // multiply with the fracture mobility of the upstream vertex
            fractureIntrinsicPermeability_.mv(pGrad, fractureFilterVelocity_[phaseIdx]);
            fractureFilterVelocity_[phaseIdx] *= -up.fractureMobility(phaseIdx);

            // divide the volume flux by two. This is required because
            // a fracture is always shared by two sub-control-volume
            // faces.
            fractureVolumeFlux_[phaseIdx] =
                (fractureFilterVelocity_[phaseIdx]*distDirection)
                * fractureWidth_/2.0;
        }
    }

public:
    const DimMatrix &fractureIntrinsicPermeability() const
    { return fractureIntrinsicPermeability_; }

    Scalar fractureVolumeFlux(int phaseIdx) const
    { return fractureVolumeFlux_[phaseIdx]; }

    Scalar fractureWidth() const
    { return fractureWidth_; }

    const DimVector &fractureFilterVelocity(int phaseIdx) const
    { return fractureFilterVelocity_[phaseIdx]; }

private:
    DimMatrix fractureIntrinsicPermeability_;
    DimVector fractureFilterVelocity_[numPhases];
    Scalar fractureVolumeFlux_[numPhases];
    Scalar fractureWidth_;
};

} // end namepace

#endif
