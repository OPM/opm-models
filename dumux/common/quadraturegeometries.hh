// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
/**
 * \file
 *
 * \brief Geometries for quadrature which allow to specify all corners.
 */
#ifndef DUMUX_QUADRATURE_GEOMETRIES_HH
#define DUMUX_QUADRATURE_GEOMETRIES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {
/*!
 * \brief Quadrilateral quadrature geometry.
 */
template <class Scalar, int dim>
class QuadrialteralQuadratureGeometry
{
public:
    enum { numCorners = (1 << dim) }; 

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

    Dune::GeometryType type() const
    { return Dune::GeometryType(Dune::GeometryType::cube, dim); }
    
    template <class CornerContainer>
    void setCorners(const CornerContainer &corners)
    {
        int cornerIdx = 0;
        for (const auto &corner : corners) {
            for (int j = 0; j < dim; ++ j)
                corners_[cornerIdx][j] = corner[j];
            ++ cornerIdx;
            if (cornerIdx == numCorners)
                break;
        }
        assert(cornerIdx == numCorners);

        center_ = 0;
        for (int cornerIdx = 0; cornerIdx < numCorners; ++ cornerIdx)
            center_ += corners_[cornerIdx];
        center_ /= numCorners;
    }

    /*!
     * \brief Returns the center of weight of the polyhedron.
     */
    const GlobalPosition &center() const
    { return center_; }
    
    /*!
     * \brief Convert a local coordinate into a global one.
     */
    GlobalPosition global(const LocalPosition &localPos) const
    {
        GlobalPosition globalPos(0.0);

        for (int cornerIdx = 0; cornerIdx < numCorners; ++ cornerIdx)
            globalPos.axpy(cornerWeight(localPos, cornerIdx), corners_[cornerIdx]);

        return globalPos;
    }

    /*!
     * \brief Returns the Jacobian matrix of the local to global
     *        mapping at a given local position.
     */
    void jacobian(Dune::FieldMatrix<Scalar, dim, dim> &jac, const LocalPosition &localPos) const
    {
        jac = 0.0;
        for (int cornerIdx = 0; cornerIdx < numCorners; ++ cornerIdx) {
            for (int k = 0; k < dim; ++k) {
                Scalar dWeight_dk = (cornerIdx & (1 << k)) ? 1 : -1;
                for (int j = 0; j < dim; ++j) {
                    if (k != j) {
                        if (cornerIdx & (1 << j))
                            dWeight_dk *= localPos[j];
                        else
                            dWeight_dk *= 1 - localPos[j];;
                    }
                }

                jac[k].axpy(dWeight_dk, corners_[cornerIdx]); 
            }
        }
    }

    /*!
     * \brief Return the determinant of the Jacobian of the mapping
     *        from local to global coordinates at a given local
     *        position.
     */
    Scalar integrationElement(const LocalPosition &localPos) const
    {
        Dune::FieldMatrix<Scalar, dim, dim> jac;
        jacobian(jac, localPos);
        return jac.determinant();
    }

    /*!
     * \brief Return the weight of an individual corner for the local
     *        to global mapping.
     */
    Scalar cornerWeight(const LocalPosition &localPos, int cornerIdx) const
    {
        GlobalPosition globalPos(0.0);

        // this code is based on the Q1 finite element code from
        // dune-localfunctions
        Scalar weight = 1.0;
        for (int j=0; j < dim; ++ j)
            weight *= (cornerIdx & (1<<j)) ? localPos[j] : (1 - localPos[j]);

        return weight;
    }

private:
    GlobalPosition corners_[numCorners];
    GlobalPosition center_;
};

}

#endif // DUMUX_QUADRATURE_GEOMETRY_HH
