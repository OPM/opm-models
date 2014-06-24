/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Ewoms::DiscreteFractureExtensiveQuantities
 */
#ifndef EWOMS_DISCRETE_FRACTURE_EXTENSIVE_QUANTITIES_HH
#define EWOMS_DISCRETE_FRACTURE_EXTENSIVE_QUANTITIES_HH

#include <ewoms/models/immiscible/immiscibleextensivequantities.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup DiscreteFractureModel
 * \ingroup ExtensiveQuantities
 *
 * \brief This class expresses all intensive quantities of the discrete fracture model.
 */
template <class TypeTag>
class DiscreteFractureExtensiveQuantities : public ImmiscibleExtensiveQuantities<TypeTag>
{
    typedef ImmiscibleExtensiveQuantities<TypeTag> ParentType;

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
     * \copydoc MultiPhaseBaseExtensiveQuantities::update()
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);

        const auto &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const auto &stencil = elemCtx.stencil(timeIdx);
        const auto &scvf = stencil.interiorFace(scvfIdx);
        int insideScvIdx = scvf.interiorIndex();
        int outsideScvIdx = scvf.exteriorIndex();

        int globalI = elemCtx.globalSpaceIndex(insideScvIdx, timeIdx);
        int globalJ = elemCtx.globalSpaceIndex(outsideScvIdx, timeIdx);
        const auto &fractureMapper = elemCtx.problem().fractureMapper();
        if (!fractureMapper.isFractureEdge(globalI, globalJ))
            // do nothing if no fracture goes though the current edge
            return;

        // average the intrinsic permeability of the fracture
        const auto &intQuantsInside = elemCtx.intensiveQuantities(insideScvIdx, timeIdx);
        const auto &intQuantsOutside = elemCtx.intensiveQuantities(outsideScvIdx, timeIdx);
        elemCtx.problem().meanK(fractureIntrinsicPermeability_,
                                intQuantsInside.fractureIntrinsicPermeability(),
                                intQuantsOutside.fractureIntrinsicPermeability());

        auto distDirection = elemCtx.pos(outsideScvIdx, timeIdx);
        distDirection -= elemCtx.pos(insideScvIdx, timeIdx);
        distDirection /= distDirection.two_norm();

        const auto &problem = elemCtx.problem();
        fractureWidth_ = problem.fractureWidth(elemCtx, insideScvIdx,
                                               outsideScvIdx, timeIdx);
        assert(fractureWidth_ < scvf.area());

        const auto &evalPointExtQuants = elemCtx.evalPointExtensiveQuantities(scvfIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto &pGrad = extQuants.potentialGrad(phaseIdx);

            int upstreamIdx = evalPointExtQuants.upstreamIndex(phaseIdx);
            const auto &up = elemCtx.intensiveQuantities(upstreamIdx, timeIdx);

            // multiply with the fracture mobility of the upstream vertex
            fractureIntrinsicPermeability_.mv(pGrad,
                                              fractureFilterVelocity_[phaseIdx]);
            fractureFilterVelocity_[phaseIdx] *= -up.fractureMobility(phaseIdx);

            // divide the volume flux by two. This is required because
            // a fracture is always shared by two sub-control-volume
            // faces.
            fractureVolumeFlux_[phaseIdx] = 0;
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                fractureVolumeFlux_[phaseIdx] +=
                    (fractureFilterVelocity_[phaseIdx][dimIdx] * distDirection[dimIdx])
                    * (fractureWidth_ / 2.0) / scvf.area();
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

} // namespace Ewoms

#endif
