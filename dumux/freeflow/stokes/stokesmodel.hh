// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2010-2012 by Katherina Baber, Klaus Mosthaf               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
#ifndef DUMUX_STOKES_MODEL_HH
#define DUMUX_STOKES_MODEL_HH

#include "stokeslocalresidual.hh"
#include "stokesnewtoncontroller.hh"
#include "stokesproblem.hh"
#include "stokesproperties.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dune/common/fvector.hh>

#include <dumux/boxmodels/common/boxmodel.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \brief Adaption of the box scheme to the Stokes model.
 *
 * This model implements laminar Stokes flow of a single fluid, solving a momentum balance:
 * \f[
\frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \left(p_g {\bf {I}}
- \mu_g \left(\boldsymbol{\nabla} \boldsymbol{v}_g
+ \boldsymbol{\nabla} \boldsymbol{v}_g^T\right)\right)
- \varrho_g {\bf g} = 0,
 * \f]
 *
 * and the mass balance equation:
 * \f[
\frac{\partial \varrho_g}{\partial t} + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class StokesModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        if (pvIdx == Indices::pressureIdx)
            oss << "pressure";
        else if (Indices::velocity0Idx <= pvIdx && pvIdx < Indices::velocity0Idx + dim)
            oss << "velocity_" << pvIdx - Indices::velocity0Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        if (eqIdx == Indices::massBalanceIdx) {
            oss << "continuity";
        }
        else if (Indices::momentum0Idx <= eqIdx && eqIdx < Indices::momentum0Idx + dim)
            oss << "momentum_" << eqIdx - Indices::momentum0Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return 1.0/1e5;
        return 1;
    }

    //! \copydoc BoxModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer) const
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        // iterate over grid
        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updateAll(*elemIt);

            int numLocalVerts = elemIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/i, /*timeIdx=*/0);
                const auto &volVars = elemCtx.volVars(/*spaceIdx=*/i, /*timeIdx=*/0);
                const auto &fs = volVars.fluidState();
                
                pN[globalIdx] = fs.pressure(phaseIdx);
                delP[globalIdx] = fs.pressure(phaseIdx) - 1e5;
                rho[globalIdx] = fs.density(phaseIdx);
                mu[globalIdx] = fs.viscosity(phaseIdx);
                velocity[globalIdx] = volVars.velocity();
            };
        }

        std::ostringstream tmp;

        tmp.str(""); tmp << "pressure_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(pN, tmp.str());

        tmp.str(""); tmp << "delta pressure_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(delP, tmp.str());

        tmp.str(""); tmp << "density_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(rho, tmp.str());

        tmp.str(""); tmp << "viscosity_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(mu, tmp.str());

        tmp.str(""); tmp << "velocity_" << FluidSystem::phaseName(phaseIdx);
        writer.attachVertexData(velocity, tmp.str(), dim);
    }
};
}

#include "stokespropertydefaults.hh"

#endif
