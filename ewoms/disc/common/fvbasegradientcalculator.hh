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
 * \copydoc Ewoms::FvBaseGradientCalculator
 */
#ifndef EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH
#define EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH

#include "fvbaseproperties.hh"

#include <opm/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Ewoms {
template<class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief This class calculates gradients of arbitrary quantities at
 *        flux integration points using the two-point approximation scheme
 */
template<class TypeTag>
class FvBaseGradientCalculator
{

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // maximum number of flux approximation points. to calculate this,
    // we assume that the geometry with the most pointsq is a cube.
    enum { maxFap = 2 << dim };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Evaluation, dim> EvalDimVector;

    static_assert(std::is_same<Evaluation, Scalar>::value ||
                  std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value,
                  "So far, automatic differentiation is only available for the "
                  "element-centered finite volume discretization!");

public:
    /*!
     * \brief Register all run-time parameters for the gradient calculator
     *        of the base class of the discretization.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Precomputes the common values to calculate gradients and values of
     *        quantities at every interior flux approximation point.
     *
     * \param elemCtx The current execution context
     * \param timeIdx The index used by the time discretization.
     */
    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext& elemCtx OPM_UNUSED, unsigned timeIdx OPM_UNUSED)
    { /* noting to do */ }

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
    template <class QuantityCallback>
    auto calculateValue(const ElementContext& elemCtx,
                        unsigned fapIdx,
                        const QuantityCallback& quantityCallback) const
        -> typename std::remove_reference<decltype(quantityCallback.operator()(0))>::type
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.interiorFace(fapIdx);

        // calculate the distances of the position of the interior and of the exterior
        // finite volume to the position of the integration point.
        const auto& normal = face.normal();
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).globalPos();
        const auto& exteriorPos = stencil.subControlVolume(face.exteriorIndex()).globalPos();
        const auto& integrationPos = face.integrationPos();

        Scalar interiorDistance = 0.0;
        Scalar exteriorDistance = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            interiorDistance +=
                (interiorPos[dimIdx] - integrationPos[dimIdx])
                * normal[dimIdx];

            exteriorDistance +=
                (exteriorPos[dimIdx] - integrationPos[dimIdx])
                * normal[dimIdx];
        }

        interiorDistance = std::sqrt(std::abs(interiorDistance));
        exteriorDistance = std::sqrt(std::abs(exteriorDistance));

        // use the average weighted by distance...
        auto value(quantityCallback(face.interiorIndex()));
        value *= interiorDistance;
        auto tmp(quantityCallback(face.exteriorIndex()));
        tmp *= exteriorDistance;
        value += tmp;
        value /= interiorDistance + exteriorDistance;

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
     *               of the quantity given the index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateGradient(EvalDimVector& quantityGrad,
                           const ElementContext& elemCtx,
                           unsigned fapIdx,
                           const QuantityCallback& quantityCallback) const
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.interiorFace(fapIdx);

        const auto& exteriorPos = stencil.subControlVolume(face.exteriorIndex()).globalPos();
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).globalPos();

        // this is slightly hacky because the derivatives of the quantity for the
        // exterior DOF are thrown away and this code thus assumes that the exterior DOF
        // is not a primary degree of freedom. Basically this means that two-point flux
        // approximation scheme is used. A way to fix this in a conceptionally elegant
        // way would be to introduce the concept of extensive evaluations, but
        // unfortunately this makes things quite a bit slower. (and also quite a bit
        // harder to comprehend :/ )
        Evaluation deltay =
            Toolbox::value(quantityCallback(face.exteriorIndex()))
            - quantityCallback(face.interiorIndex());

        Scalar distSquared = 0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            quantityGrad[dimIdx] = deltay;

            Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
            distSquared += tmp*tmp;
        }

        // divide the gradient by the squared distance between the positions of the
        // control volumes: the gradient is the normalized directional vector between the
        // two centers times the ratio of the difference of the values and their
        // distance, i.e., d/abs(d) * delta y / abs(d) = d*delta y / abs(d)^2.
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = exteriorPos[dimIdx] - interiorPos[dimIdx];
            quantityGrad[dimIdx] *= tmp/distSquared;
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
     *               of the quantity given the index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    auto calculateBoundaryValue(const ElementContext& elemCtx OPM_UNUSED,
                                unsigned fapIdx OPM_UNUSED,
                                const QuantityCallback& quantityCallback)
        -> decltype(quantityCallback.boundaryValue())
    { return quantityCallback.boundaryValue(); }

    /*!
     * \brief Calculates the gradient of an arbitrary quantity at any
     *        flux approximation point on the boundary.
     *
     * Boundary gradients are always calculated using the two-point
     * approximation.
     *
     * \param elemCtx The current execution context
     * \param faceIdx The local index of the flux approximation point
     *                in the current element's stencil.
     * \param quantityCallback A callable object returning the value
     *               of the quantity at an index of a degree of
     *               freedom
     */
    template <class QuantityCallback>
    void calculateBoundaryGradient(EvalDimVector& quantityGrad,
                                   const ElementContext& elemCtx,
                                   unsigned faceIdx,
                                   const QuantityCallback& quantityCallback) const
    {
        const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto& face = stencil.boundaryFace(faceIdx);

        auto deltay = quantityCallback.boundaryValue() - quantityCallback(face.interiorIndex());

        const auto& boundaryFacePos = face.integrationPos();
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).center();

        Scalar distSquared = 0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = boundaryFacePos[dimIdx] - interiorPos[dimIdx];
            distSquared += tmp*tmp;

            quantityGrad[dimIdx] = deltay;
        }

        // divide the gradient by the squared distance between the center of the
        // sub-control and the center of the boundary face: the gradient is the
        // normalized directional vector between the two centers times the ratio of the
        // difference of the values and their distance, i.e., d/abs(d) * deltay / abs(d)
        // = d*deltay / abs(d)^2.
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            Scalar tmp = boundaryFacePos[dimIdx] - interiorPos[dimIdx];
            quantityGrad[dimIdx] *= tmp/distSquared;
        }
    }
};
} // namespace Ewoms

#endif
