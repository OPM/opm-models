// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief This file contains parts to calculate the diffusive flux in
 *        the compositional NCP model
 */
#ifndef DUMUX_NCP_DIFFUSION_HH
#define DUMUX_NCP_DIFFUSION_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/boxmodels/ncp/ncpproperties.hh>

namespace Dumux {

template <class TypeTag, bool enableDiffusion>
class NcpDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;


    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> DiffMatrix;
    typedef Dune::FieldVector<Scalar, numComponents> DiffVector;
    typedef Dune::FieldVector<Scalar, numComponents> CompVector;

public:
    static void flux(CompVector &fluxes,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx,
                     int phaseIdx,
                     Scalar molarDensity)
    {
        if (phaseIdx == gPhaseIdx)
            gasFlux_(fluxes, elemCtx, scvfIdx, timeIdx, molarDensity);
        else if (phaseIdx == lPhaseIdx)
            liquidFlux_(fluxes, elemCtx, scvfIdx, timeIdx, molarDensity);
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Invalid phase index: " << phaseIdx);
    }

protected:
    static void liquidFlux_(CompVector &fluxes,
                            const ElementContext &elemCtx,
                            int scvfIdx,
                            int timeIdx,
                            Scalar molarDensity)
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            // TODO: tensorial diffusion coefficients
            Scalar xGrad = fluxVars.moleFracGrad(lPhaseIdx, compIdx)*normal;
            fluxes[compIdx] =
                - xGrad *
                molarDensity *
                normal.two_norm() * // because we want a mole flux and not an area specific flux
                fluxVars.porousDiffCoeffL(compIdx) ;
        }
    }

    static void gasFlux_(CompVector &fluxes,
                         const ElementContext &elemCtx,
                         int scvfIdx,
                         int timeIdx,
                         Scalar molarDensity)
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;

        // Stefan-Maxwell equation
        //
        // See: R. Reid, et al.: "The Properties of Liquids and
        // Gases", 4th edition, 1987, McGraw-Hill, p 596

        // TODO: tensorial diffusion coefficients
        DiffMatrix M(0);

        for (int i = 0; i < numComponents - 1; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                Scalar Dij = fluxVars.porousDiffCoeffG(i, j);
                if (Dij) {
                    M[i][j] += fluxVars.moleFraction(gPhaseIdx, i) / Dij;
                    M[i][i] -= fluxVars.moleFraction(gPhaseIdx, j) / Dij;
                }
            }
        };

        for (int i = 0; i < numComponents; ++i) {
            M[numComponents - 1][i] = 1.0;
        }

        DiffVector rightHandSide ; // see source cited above
        for (int i = 0; i < numComponents - 1; ++i) {
            rightHandSide[i] = molarDensity*(fluxVars.moleFracGrad(gPhaseIdx, i)*normal);
        }
        rightHandSide[numComponents - 1] = 0.0;

        M.solve(fluxes, rightHandSide);
    }

    // return whether a concentration can be assumed to be a trace
    // component in the context of diffusion
    static Scalar isTraceComp_(Scalar x)
    { return x < 0.5/numComponents; }
};

/*!
 * \brief Specialization of the diffusion module for the case where
 *        diffusion is disabled.
 *
 * This class just does nothing.
 */
template <class TypeTag>
class NcpDiffusion<TypeTag, false>
{
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef Dune::FieldVector<Scalar, numComponents> CompVector;

public:
    static void flux(CompVector &fluxes,
                     const ElementContext &fluxVars,
                     int scvfIdx,
                     int timeIdx,
                     int phaseIdx,
                     Scalar totalConcentration)
    { fluxes = 0; }
};

}

#endif // DUMUX_NCP_DIFFUSION_HH
