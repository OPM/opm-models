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
 * \copydoc Ewoms::StokesModel
 */
#ifndef EWOMS_STOKES_MODEL_HH
#define EWOMS_STOKES_MODEL_HH

#include "stokeslocalresidual.hh"
#include "stokesproblem.hh"
#include "stokesproperties.hh"

#include <ewoms/disc/vcfv/vcfvmodel.hh>

#include <dune/common/fvector.hh>

#include <ewoms/disc/vcfv/vcfvmodel.hh>

namespace Ewoms {

/*!
 * \ingroup VCFVStokesModel
 *
 * \brief Fully-implicit discretization of the Navier-Stokes equations
 *        using the vertex-centered finite volume scheme.
 *
 * This model implements Navier-Stokes flow of a single fluid. By
 * default, it solves the momentum balance of the time-dependent Stokes
 * equations, i.e.
 * \f[
 * \frac{\partial \left(\varrho\,\mathbf{v}\right)}{\partial t}
 * + \boldsymbol{\nabla} p
 * - \nabla \cdot
 * \left(
 * \mu \left(\boldsymbol{\nabla} \mathbf{v} + \boldsymbol{\nabla} \mathbf{v}^T\right)
 * \right)
 * - \varrho\,\mathbf{g}
 * = 0\;,
 * \f]
 * and the mass balance equation
 * \f[
 * \frac{\partial \varrho}{\partial t}
 * + \nabla \cdot\left(\varrho\,\mathbf{v}\right)
 * - q
 * = 0 \;.
 * \f]
 *
 * If the property \c EnableNavierStokes is set to \c true, an
 * additional convective momentum flux term (Navier term) gets
 * included into the momentum conservation equations which allows to
 * capture turbolent flow regimes. This additional term is given by
 *  \f[
 * \varrho \left(\mathbf{v} \cdot \boldsymbol{\nabla} \right) \mathbf{v} \;.
 * \f]
 *
 * These equations are discretized by a fully-coupled vertex-centered
 * finite volume scheme in space and using the implicit Euler method
 * in time. Be aware, that this discretization scheme is quite
 * unstable for the Navier-Stokes equations and quickly leads to
 * unphysical oscillations in the calculated solution. We intend to
 * use a more appropriate discretization scheme in the future, though.
 */
template<class TypeTag>
class StokesModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dimWorld = GridView::dimensionworld };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { numComponents = FluidSystem::numComponents };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Returns true iff a fluid phase is used by the model.
     *
     * \param phaseIdxQueried The index of the fluid phase in question
     */
    bool phaseIsConsidered(int phaseIdxQueried) const
    { return phaseIdxQueried == phaseIdx; }

    /*!
     * \copydoc VcfvModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;
        if (pvIdx == Indices::pressureIdx)
            oss << "pressure";
        else if (Indices::moleFrac1Idx <= pvIdx && pvIdx < Indices::moleFrac1Idx + numComponents - 1)
            oss << "moleFraction^" << FluidSystem::componentName(pvIdx - Indices::moleFrac1Idx + 1);
        else if (Indices::velocity0Idx <= pvIdx && pvIdx < Indices::velocity0Idx + dimWorld)
            oss << "velocity_" << pvIdx - Indices::velocity0Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc VcfvModel::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(eqIdx - Indices::conti0EqIdx);
        else if (Indices::momentum0EqIdx <= eqIdx && eqIdx < Indices::momentum0EqIdx + dimWorld)
            oss << "momentum_" << eqIdx - Indices::momentum0EqIdx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc VcfvModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        // for stokes flow the pressure gradients are often quite
        // small, so we need higher precision for pressure. TODO: find
        // a good weight for the pressure.
        if (Indices::pressureIdx == pvIdx)
            return 1.0/this->solution(/*timeIdx=*/0)[globalVertexIdx][Indices::pressureIdx];

        return 1;
    }

    /*!
     * \copydoc VcfvModel::addOutputVtkFields
     */
    template <class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer) const
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dimWorld);
        ScalarField &pressure = *writer.allocateManagedBuffer(numVertices);
        ScalarField &density = *writer.allocateManagedBuffer(numVertices);
        ScalarField &temperature = *writer.allocateManagedBuffer(numVertices);
        ScalarField &viscosity = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<double, dimWorld> (numVertices);
        ScalarField *moleFraction[numComponents];
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            moleFraction[compIdx] = writer.allocateManagedBuffer(numVertices);

        // iterate over grid
        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updateAll(*elemIt);

            int numScv = elemCtx.numScv();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                int globalIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/scvIdx, /*timeIdx=*/0);
                const auto &volVars = elemCtx.volVars(/*spaceIdx=*/scvIdx, /*timeIdx=*/0);
                const auto &fluidState = volVars.fluidState();

                pressure[globalIdx] = fluidState.pressure(phaseIdx);
                density[globalIdx] = fluidState.density(phaseIdx);
                temperature[globalIdx] = fluidState.temperature(phaseIdx);
                viscosity[globalIdx] = fluidState.viscosity(phaseIdx);
                velocity[globalIdx] = volVars.velocityCenter();

                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    (*moleFraction[compIdx])[globalIdx] = fluidState.moleFraction(phaseIdx, compIdx);
            };
        }

        std::ostringstream tmp;

        tmp.str(""); tmp << "pressure_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(pressure, tmp.str());

        tmp.str(""); tmp << "density_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(density, tmp.str());

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            tmp.str(""); tmp << "moleFraction_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
            writer.attachVertexData(*moleFraction[compIdx], tmp.str());
        }

        tmp.str(""); tmp << "temperature_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(temperature, tmp.str());

        tmp.str(""); tmp << "viscosity_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(viscosity, tmp.str());

        tmp.str(""); tmp << "velocity_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(velocity, tmp.str(), dimWorld);
    }
};
}

#include "stokespropertydefaults.hh"

#endif
