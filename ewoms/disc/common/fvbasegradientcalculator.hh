/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::FvBaseGradientCalculator
 */
#ifndef EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH
#define EWOMS_FV_BASE_GRADIENT_CALCULATOR_HH

#include "fvbaseproperties.hh"

#include <dune/common/fvector.hh>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief This class calculates gradients of arbitrary quantities at
 *        flux integration points using the two-point approximation scheme
 */
template<class TypeTag>
class FvBaseGradientCalculator
{

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // maximum number of flux approximation points. to calculate this,
    // we assume that the geometry with the most pointsq is a cube.
    enum { maxFap = 2 << dim };

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \brief Register all run-time parameters for the gradient calculator
     *        of the base class of the discretization.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Precomputes the common values to calculate gradients and
     *        values of quantities at every interior flux
     *        approximation point.
     *
     * \param elemCtx The current execution context
     * \param timeIdx The index used by the time discretization.
     */
    template <bool prepareValues = true, bool prepareGradients = true>
    void prepare(const ElementContext &elemCtx, int timeIdx)
    {
        const auto &stencil = elemCtx.stencil(timeIdx);
        for (int fapIdx = 0; fapIdx < stencil.numInteriorFaces(); ++ fapIdx) {
            const auto &scvf = stencil.interiorFace(fapIdx);
            const auto &normal = scvf.normal();
            const auto &interiorPos = stencil.subControlVolume(scvf.interiorIndex()).globalPos();
            const auto &exteriorPos = stencil.subControlVolume(scvf.exteriorIndex()).globalPos();

            interiorDistance_[fapIdx] = 0;
            exteriorDistance_[fapIdx] = 0;
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                interiorDistance_[fapIdx] +=
                    (interiorPos[dimIdx] - scvf.integrationPos()[dimIdx])
                    * normal[dimIdx];

                exteriorDistance_[fapIdx] +=
                    (exteriorPos[dimIdx] - scvf.integrationPos()[dimIdx])
                    * normal[dimIdx];
            }

            interiorDistance_[fapIdx] = std::abs(interiorDistance_[fapIdx]);
            exteriorDistance_[fapIdx] = std::abs(exteriorDistance_[fapIdx]);
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
        const auto &face = elemCtx.stencil(/*timeIdx=*/0).interiorFace(fapIdx);

        // average weighted by distance to DOF coordinate...
        QuantityType value(quantityCallback(face.interiorIndex()));
        value *= interiorDistance_[fapIdx];
        QuantityType tmp(quantityCallback(face.exteriorIndex()));
        tmp *= exteriorDistance_[fapIdx];
        value += tmp;
        value /= interiorDistance_[fapIdx] + exteriorDistance_[fapIdx];

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
    void calculateGradient(DimVector &quantityGrad,
                           const ElementContext &elemCtx,
                           int fapIdx,
                           const QuantityCallback &quantityCallback) const
    {
        const auto &stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto &face = stencil.interiorFace(fapIdx);
        const auto &normal = face.normal();

        Scalar dy = quantityCallback(face.exteriorIndex()) - quantityCallback(face.interiorIndex());
        Scalar dx = exteriorDistance_[fapIdx] + interiorDistance_[fapIdx];

        for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            quantityGrad[dimIdx] = normal[dimIdx]*dy/dx;
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
    Scalar calculateBoundaryValue(const ElementContext &elemCtx,
                                  int fapIdx,
                                  const QuantityCallback &quantityCallback)
    { return quantityCallback.boundaryValue(); }

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
    {
        const auto &stencil = elemCtx.stencil(/*timeIdx=*/0);
        const auto &face = stencil.boundaryFace(fapIdx);
        const auto &normal = face.normal();

        Scalar dy = quantityCallback.boundaryValue()
            - quantityCallback(face.interiorIndex());

        DimVector distVec(face.integrationPos());
        distVec -= stencil.subControlVolume(face.interiorIndex()).center();

        // we assume that the value on the boundary belongs equivalent
        // to a finite volume which is equivalent to the finite volume
        // on the interior of the face, but mirrored on the boundary
        // face. This means the distance is twice the one between the
        // integration point and the center of the FV on the interior
        // of the domain.
        Scalar dx = 2*distVec.two_norm();

        quantityGrad = normal;
        quantityGrad *= dy/dx;
    }

private:
    // distance [m] of the the flux approximation point to the center
    // of the control volume to the inside of the face
    Scalar interiorDistance_[maxFap];

    // distance [m] of the the flux approximation point to the center
    // of the control volume to the outside of the face
    Scalar exteriorDistance_[maxFap];
};
} // namespace Ewoms

#endif
