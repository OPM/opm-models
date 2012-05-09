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
#include "config.h"

#include<dune/common/exceptions.hh>
#include<dune/grid/common/quadraturerules.hh>

#include <dumux/common/quadraturegeometries.hh>

template <class GlobalCoord>
typename GlobalCoord::field_type f(const GlobalCoord &pos)
{
    typename GlobalCoord::field_type result = 1;
    for (int i = 0; i < GlobalCoord::dimension; ++i)
        result *= pos[i];
    return result;
}

int main(int argc, char** argv)
{
    const int dim = 3;
    typedef double Scalar;
    typedef Dumux::QuadrialteralQuadratureGeometry<Scalar, dim> QuadratureGeom;
    typedef QuadratureGeom::LocalPosition LocalPosition;
    typedef QuadratureGeom::GlobalPosition GlobalPosition;

    QuadratureGeom foo;
    
    Scalar corners[][3] = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 1, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 0, 1, 1 },
        { 1, 1, 1 }
    };
    foo.setCorners(corners);

    std::cout << "testing local to global mapping...\n";
    int n = 100;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                LocalPosition localPos;
                
                localPos[0] = Scalar(i)/(n - 1);
                localPos[1] = Scalar(j)/(n - 1);
                localPos[2] = Scalar(k)/(n - 1);
                
                GlobalPosition globalPos = foo.global(localPos);
                
                GlobalPosition diff(localPos);
                diff -= globalPos;
                assert(diff.two_norm() < 1e-10);
            }
        }
    }

    std::cout << "testing quadrature rule...\n";
    Dune::GeometryType geomType = foo.type();
    static const int quadratureOrder = 2;
    const auto &rule = Dune::QuadratureRules<Scalar,dim>::rule(geomType, quadratureOrder);
    
    // compute approximate integral
    Scalar result=0;
    for (auto it = rule.begin(); it != rule.end(); ++ it)
    {
          Scalar fval = f(foo.global(it->position()));
          Scalar weight = it->weight();
          Scalar detjac = foo.integrationElement(it->position());
          assert(std::abs(1 - detjac) < 1e-10);
          result += fval * weight * detjac;
    }
    
    std::cout << "result: " << result << " (expected value: " << 1.0/(1 << dim) << ")\n";

    return 0;
}
