// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
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
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate the
 *        component fluxes over a face of a finite volume.
 *
 * This means concentration gradients, diffusion coefficients, mass fractions, etc.
 * at the integration point.
 */
#ifndef DUMUX_STOKES2C_FLUX_VARIABLES_HH
#define DUMUX_STOKES2C_FLUX_VARIABLES_HH

#include <dumux/freeflow/stokes/stokesfluxvariables.hh>
#include <dumux/common/math.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Stokes2cIndices); //!< Enumerations for the compositional stokes models
}

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for the compositional Stokes model.
 *
 * This means concentration gradients, diffusion coefficient, mass fractions, etc.
 * at the integration point of a SCV or boundary face.
 */
template <class TypeTag>
class Stokes2cFluxVariables : public StokesFluxVariables<TypeTag>
{
    typedef StokesFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;

    enum { dim = GridView::dimension };
    enum { lCompIdx = Indices::lCompIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIndex) };

    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;

public:
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx, bool isBoundaryFace = false)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx, isBoundaryFace);

        diffusionCoeffAtIP_ = Scalar(0);
        moleFractionGradAtIP_ = Scalar(0);

        const auto &fvElemGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.subContVolFace[scvfIdx];

        // calculate gradients and secondary variables at IPs
        for (int idx = 0;
             idx < elemCtx.numScv();
             idx++) // loop over vertices of the element
        {
            const auto &volVars = elemCtx.volVars(idx, timeIdx);

            diffusionCoeffAtIP_ += 
                volVars.diffusionCoeff()
                * scvf.shapeValue[idx];

            // the gradient of the mass fraction at the IP
            for (int dimIdx=0; dimIdx<dim; ++dimIdx)
            {
                moleFractionGradAtIP_ +=
                    scvf.grad[idx][dimIdx] *
                    volVars.fluidState().moleFraction(phaseIdx, lCompIdx);
            }
        };

        Valgrind::CheckDefined(diffusionCoeffAtIP_);
        Valgrind::CheckDefined(moleFractionGradAtIP_);
    }

    /*!
     * \brief Return the gradient of the mole fraction at the integration point.
     */
    const ScalarGradient &moleFractionGradAtIP() const
    { return moleFractionGradAtIP_; }
    
    /*!
     * \brief Return the diffusion coefficient at the integration point.
     */
    Scalar diffusionCoeffAtIP() const
    { return diffusionCoeffAtIP_; }

protected:
    Scalar diffusionCoeffAtIP_;
    ScalarGradient moleFractionGradAtIP_;
};

} // end namespace

#endif
