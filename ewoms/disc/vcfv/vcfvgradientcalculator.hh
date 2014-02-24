/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \copydoc Ewoms::VcfvGradientCalculator
 */
#ifndef EWOMS_VCFV_GRADIENT_CALCULATOR_HH
#define EWOMS_VCFV_GRADIENT_CALCULATOR_HH

#include "vcfvproperties.hh"

#include <ewoms/disc/common/fvbasegradientcalculator.hh>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/common/fvector.hh>

#include <vector>

namespace Ewoms {
/*!
 * \ingroup VcfvDiscretization
 *
 * \brief This class calculates gradients of arbitrary quantities at
 *        flux integration points for the vertex centered finite
 *        volume (VCFV) discretization
 */
template<class TypeTag>
class VcfvGradientCalculator : public FvBaseGradientCalculator<TypeTag>
{
    typedef FvBaseGradientCalculator<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dim = GridView::dimension };

    // set the maximum number of degrees of freedom and the maximum
    // number of flux approximation points per elements. For this, we
    // assume cubes as the type of element with the most vertices...
    enum { maxDof = (2 << dim) };
    enum { maxFap = maxDof };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;

public:
    /*!
     * \brief Register all run-time parameters for the gradient calculator
     *        of the VCVF discretization.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, UseTwoPointGradients,
                             "Use two-point gradients instead of P1-finite element ones");
    }

    /*!
     * \brief Precomputes the common values to calculate gradients and
     *        values of quantities at any flux approximation point.
     *
     * \param elemCtx The current execution context
     */
    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext &elemCtx, int timeIdx)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, UseTwoPointGradients)) {
            ParentType::template prepare<prepareValues, prepareGradients>(elemCtx, timeIdx);
            return;
        }

        const auto &stencil = elemCtx.stencil(timeIdx);

        const LocalFiniteElement& localFE = feCache_.get(elemCtx.element().type());
        localFiniteElement_ = &localFE;

        // loop over all face centeres
        for (int faceIdx = 0; faceIdx < stencil.numInteriorFaces(); ++faceIdx) {
            const auto &localFacePos = stencil.interiorFace(faceIdx).localPos();

            // Evaluate the P1 shape functions and their gradients at all
            // flux approximation points.
            if (prepareValues)
                localFE.localBasis().evaluateFunction(localFacePos, p1Value_[faceIdx]);

            if (prepareGradients) {
                // first, get the shape function's gradient in local coordinates
                std::vector<ShapeJacobian> localGradient;
                localFE.localBasis().evaluateJacobian(localFacePos, localGradient);

                // convert to a gradient in global space by
                // multiplying with the inverse transposed jacobian of
                // the position
                const auto &geom = elemCtx.element().geometry();
                const auto &jacInvT =
                    geom.jacobianInverseTransposed(localFacePos);

                int numVertices = elemCtx.numDof(timeIdx);
                for (int vertIdx = 0; vertIdx < numVertices; vertIdx++) {
                    jacInvT.mv(/*xVector=*/localGradient[vertIdx][0],
                               /*destVector=*/p1Gradient_[faceIdx][vertIdx]);
                }
            }
        }
    }


    /*!
     * \brief Calculates the value of an arbitrary quantity at any
     *        interior flux approximation point.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback, class QuantityType = Scalar>
    QuantityType calculateValue(const ElementContext &elemCtx,
                                int fapIdx,
                                const QuantityCallback &quantityCallback) const
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, UseTwoPointGradients)) {
            return ParentType::template calculateValue<QuantityCallback,
                                                       QuantityType>(elemCtx,
                                                                     fapIdx,
                                                                     quantityCallback);
        }

        // If the user does not want to use two-point gradients, we
        // use P1 finite element gradients..
        QuantityType value(0.0);
        QuantityType tmp;
        for (int vertIdx = 0; vertIdx < elemCtx.numDof(/*timeIdx=*/0); ++vertIdx) {
            tmp = quantityCallback(vertIdx);
            tmp *= p1Value_[fapIdx][vertIdx];
            value += tmp;
        }
        return value;
    }

    /*!
     * \brief Calculates the gradient of an arbitrary quantity at any
     *        flux approximation point.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateGradient(DimVector &quantityGrad,
                           const ElementContext &elemCtx,
                           int fapIdx,
                           const QuantityCallback &quantityCallback) const
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, UseTwoPointGradients)) {
            ParentType::calculateGradient(quantityGrad, elemCtx, fapIdx, quantityCallback);
            return;
        }

        // If the user does not want two-point gradients, we use P1
        // finite element gradients...
        quantityGrad = 0.0;
        for (int vertIdx = 0; vertIdx < elemCtx.numDof(/*timeIdx=*/0); ++vertIdx) {
            Scalar dofVal = quantityCallback(vertIdx);

            DimVector tmp(p1Gradient_[fapIdx][vertIdx]);
            tmp *= dofVal;
            quantityGrad += tmp;
        }
    }

    /*!
     * \brief Calculates the value of an arbitrary quantity at any
     *        flux approximation point on the grid boundary.
     *
     * Boundary values are always calculated using the two-point
     * approximation.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    Scalar calculateBoundaryValue(const ElementContext &elemCtx,
                                  int fapIdx,
                                  const QuantityCallback &quantityCallback)
    { return ParentType::calculateBoundaryValue(elemCtx, fapIdx, quantityCallback); }

    /*!
     * \brief Calculates the gradient of an arbitrary quantity at any
     *        flux approximation point on the boundary.
     *
     * Boundary gradients are always calculated using the two-point
     * approximation.
     *
     * \param elemCtx The current execution context
     * \param fapIdx The local index of the flux approximation point
     *               in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateBoundaryGradient(DimVector &quantityGrad,
                                   const ElementContext &elemCtx,
                                   int fapIdx,
                                   const QuantityCallback &quantityCallback) const
    { ParentType::calculateBoundaryGradient(quantityGrad, elemCtx, fapIdx, quantityCallback); }

    static LocalFiniteElementCache& localFiniteElementCache()
    { return feCache_; }

private:
    static LocalFiniteElementCache feCache_;

    const LocalFiniteElement* localFiniteElement_;
    std::vector<Dune::FieldVector<Scalar, 1>> p1Value_[maxFap];
    DimVector p1Gradient_[maxFap][maxDof];
};

template<class TypeTag>
typename VcfvGradientCalculator<TypeTag>::LocalFiniteElementCache
VcfvGradientCalculator<TypeTag>::feCache_;
} // namespace Ewoms

#endif
