// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
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
#include <dune/common/fvector.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#endif

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

#if HAVE_DUNE_LOCALFUNCTIONS
    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;
    typedef typename LocalFiniteElement::Traits::LocalBasisType::Traits LocalBasisTraits;
    typedef typename LocalBasisTraits::JacobianType ShapeJacobian;
#else
    typedef int LocalFiniteElementCache;
    typedef int LocalFiniteElement;
#endif


public:
    /*!
     * \brief Precomputes the common values to calculate gradients and
     *        values of quantities at any flux approximation point.
     *
     * \param elemCtx The current execution context
     */
    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext &elemCtx, unsigned timeIdx)
    {
        if (GET_PROP_VALUE(TypeTag, UseTwoPointGradients)) {
            ParentType::template prepare<prepareValues, prepareGradients>(elemCtx, timeIdx);
            return;
        }

        const auto &stencil = elemCtx.stencil(timeIdx);

#if HAVE_DUNE_LOCALFUNCTIONS
        const LocalFiniteElement& localFE = feCache_.get(elemCtx.element().type());
        localFiniteElement_ = &localFE;

        // loop over all face centeres
        for (unsigned faceIdx = 0; faceIdx < stencil.numInteriorFaces(); ++faceIdx) {
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

                unsigned numVertices = elemCtx.numDof(timeIdx);
                for (unsigned vertIdx = 0; vertIdx < numVertices; vertIdx++) {
                    jacInvT.mv(/*xVector=*/localGradient[vertIdx][0],
                               /*destVector=*/p1Gradient_[faceIdx][vertIdx]);
                }
            }
        }
#else
        OPM_THROW(std::runtime_error,"VcfvGradientCalculator needs the dune-localfunctions package which has not been found!");
#endif
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
    template <class QuantityCallback, class Dummy=int>
    auto calculateValue(const ElementContext &elemCtx,
                        typename std::enable_if<!GET_PROP_VALUE(TypeTag, UseTwoPointGradients),
                                                Dummy>::type fapIdx,
                        const QuantityCallback& quantityCallback) const
        ->  typename std::remove_reference<typename QuantityCallback::ResultType>::type
    {
        static_assert(std::is_integral<Dummy>::value,
                      "The 'Dummy' template parameter must _not_ be specified explicitly."
                      "It is only required to conditionally disable this method!");

        typedef typename std::remove_reference<typename QuantityCallback::ResultType>::type QuantityConstType;
        typedef typename std::remove_const<QuantityConstType>::type QuantityType;

        // If the user does not want to use two-point gradients, we
        // use P1 finite element gradients..
        QuantityType value(0.0);
        QuantityType tmp;
        for (unsigned vertIdx = 0; vertIdx < elemCtx.numDof(/*timeIdx=*/0); ++vertIdx) {
            tmp = quantityCallback(vertIdx);
            tmp *= p1Value_[fapIdx][vertIdx];
            value += tmp;
        }
        return value;
    }

    template <class QuantityCallback, class Dummy = int>
    auto calculateValue(const ElementContext &elemCtx,
                        typename std::enable_if<GET_PROP_VALUE(TypeTag, UseTwoPointGradients),
                                                Dummy>::type fapIdx,
                        const QuantityCallback& quantityCallback) const
        ->  decltype(ParentType::calculateValue(elemCtx, fapIdx, quantityCallback))
    {
        static_assert(std::is_integral<Dummy>::value,
                      "The 'Dummy' template parameter must _not_ be specified explicitly."
                      "It is only required to conditionally disable this method!");

        return ParentType::calculateValue(elemCtx, fapIdx, quantityCallback);
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
    template <class QuantityCallback, class EvalDimVector, class Dummy = int>
    void calculateGradient(EvalDimVector &quantityGrad,
                           const ElementContext &elemCtx,
                           typename std::enable_if<!GET_PROP_VALUE(TypeTag, UseTwoPointGradients),
                                                   Dummy>::type fapIdx,
                           const QuantityCallback& quantityCallback) const
    {
        static_assert(std::is_integral<Dummy>::value,
                      "The 'Dummy' template parameter must _not_ be specified explicitly."
                      "It is only required to conditionally disable this method!");

        // If the user does not want two-point gradients, we use P1
        // finite element gradients...
        quantityGrad = 0.0;
        for (unsigned vertIdx = 0; vertIdx < elemCtx.numDof(/*timeIdx=*/0); ++vertIdx) {
            Scalar dofVal = quantityCallback(vertIdx);

            auto tmp = p1Gradient_[fapIdx][vertIdx];
            tmp *= dofVal;
            quantityGrad += tmp;
        }
    }

    template <class QuantityCallback, class EvalDimVector, class Dummy=int>
    void calculateGradient(EvalDimVector &quantityGrad,
                           const ElementContext &elemCtx,
                           typename std::enable_if<GET_PROP_VALUE(TypeTag, UseTwoPointGradients),
                                                   Dummy>::type fapIdx,
                           const QuantityCallback& quantityCallback) const
    {
        static_assert(std::is_integral<Dummy>::value,
                      "The 'Dummy' template parameter must _not_ be specified explicitly."
                      "It is only required to conditionally disable this method!");

        return ParentType::calculateGradient(quantityGrad, elemCtx, fapIdx, quantityCallback);
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
    auto calculateBoundaryValue(const ElementContext &elemCtx,
                                unsigned fapIdx,
                                const QuantityCallback &quantityCallback)
        ->  decltype(ParentType::calculateBoundaryValue(elemCtx, fapIdx, quantityCallback))
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
    template <class QuantityCallback, class EvalDimVector>
    void calculateBoundaryGradient(EvalDimVector &quantityGrad,
                                   const ElementContext &elemCtx,
                                   unsigned fapIdx,
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
