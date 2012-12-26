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
 * \copydoc Ewoms::DiscreteFractureLocalResidual
 */
#ifndef EWOMS_DISCRETE_FRACTURE_LOCAL_RESIDUAL_BASE_HH
#define EWOMS_DISCRETE_FRACTURE_LOCAL_RESIDUAL_BASE_HH

#include <ewoms/models/immiscible/immisciblelocalresidual.hh>

namespace Ewoms {

/*!
 * \ingroup DiscreteFractureVcfvModel
 *
 * \brief Calculates the local residual of the discrete fracture
 *        immiscible multi-phase model.
 */
template<class TypeTag>
class DiscreteFractureLocalResidual
    : public ImmiscibleLocalResidual<TypeTag>
{
    typedef ImmiscibleLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef VcfvEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Adds the amount all conservation quantities (e.g. phase
     *        mass) within a single fluid phase
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::vcfvScvCtxParams
     * \copydetails Doxygen::phaseIdxParam
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int scvIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        EqVector phaseStorage(0.0);
        ParentType::addPhaseStorage(phaseStorage, elemCtx, scvIdx, timeIdx, phaseIdx);

        const auto &problem = elemCtx.problem();
        const auto &fractureMapper = problem.fractureMapper();
        int globalIdx = elemCtx.globalSpaceIndex(scvIdx, timeIdx);

        if (true || !fractureMapper.isFractureVertex(globalIdx)) {
            // don't do anything in addition to the immiscible model
            // for finite volumes that do not feature fractures
            storage += phaseStorage;
            return;
        }

        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &scv = elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx];

        // reduce the matrix storage by the fracture volume
        phaseStorage *= 1 - volVars.fractureVolume()/scv.volume;

        // add the storage term inside the fractures
        const auto &fsFracture = volVars.fractureFluidState();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            phaseStorage[conti0EqIdx + phaseIdx] =
                volVars.fracturePorosity()*
                fsFracture.saturation(phaseIdx) *
                fsFracture.density(phaseIdx) *
                volVars.fractureVolume()/scv.volume;
        }

        EnergyModule::addFracturePhaseStorage(phaseStorage, volVars, scv, phaseIdx);

        // add the result to the overall storage term
        storage += phaseStorage;
    }

    /*!
     * \copydoc VcfvLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        ParentType::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        int i = fluxVars.insideIndex();
        int j = fluxVars.outsideIndex();
        int I = elemCtx.globalSpaceIndex(i, timeIdx);
        int J = elemCtx.globalSpaceIndex(j, timeIdx);
        const auto &fractureMapper = elemCtx.problem().fractureMapper();
        if (!fractureMapper.isFractureEdge(I, J))
            // do nothing if the edge from i to j is not part of a
            // fracture
            return;

        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        Scalar scvfArea = scvf.normal.two_norm();

        // advective mass fluxes of all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            // reduce the matrix mass flux by the width of the scv
            // face that is occupied by the fracture. As usual, the
            // fracture is shared between two SCVs, so the its width
            // needs to be divided by two.
            flux[conti0EqIdx + phaseIdx] *= 1 - fluxVars.fractureWidth()/(2*scvfArea);

            // vertex data of the upstream and the downstream vertices
            int upIdx = evalPointFluxVars.upstreamIndex(phaseIdx);
            const auto &up = elemCtx.volVars(upIdx, timeIdx);
            flux[conti0EqIdx + phaseIdx] +=
                fluxVars.fractureVolumeFlux(phaseIdx)
                * up.fractureFluidState().density(phaseIdx);
        }

        EnergyModule::handleFractureFlux(flux, elemCtx, scvfIdx, timeIdx);
    }
};

}

#endif
