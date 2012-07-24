// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief A preconditioner which solves the local system of equations.
 */
#ifndef DUMUX_SOLVER_PRECONDITIONER_HH
#define DUMUX_SOLVER_PRECONDITIONER_HH

namespace Dumux {
namespace Linear {

template <class Matrix, class DomainVector, class RangeVector>
class SolverPreconditioner 
    : public Dune::Preconditioner<DomainVector, RangeVector>
{
    typedef Dune::MatrixAdapter<Matrix, DomainVector, RangeVector> InnerOperator;
    typedef Dune::SeqScalarProduct<DomainVector> InnerScalarProduct;
    typedef Dune::SeqILU0<Matrix, DomainVector, RangeVector> InnerPreConditioner;
    typedef Dune::BiCGSTABSolver<DomainVector> InnerSolver;
    typedef typename DomainVector::field_type Scalar;

public:
    typedef DomainVector domain_type;
    typedef RangeVector range_type;

    enum { category = Dune::SolverCategory::overlapping };

    SolverPreconditioner(const Matrix &matrix, int order, Scalar relaxationFactor)
    {
        innerOperator_ = new InnerOperator(matrix);
        innerScalarProduct_ = new InnerScalarProduct;
        innerPreCond_ = new InnerPreConditioner(matrix, 
                                                /*relaxation=*/1.0);

        Scalar tolerance = 1e-6;
        int maxIter = 10;
        int verbosity = 0;
        innerSolver_ = new InnerSolver(*innerOperator_,
                                       *innerScalarProduct_,
                                       *innerPreCond_,
                                       tolerance,
                                       maxIter,
                                       verbosity);
    }
    
    ~SolverPreconditioner() 
    {
        delete innerSolver_;
        delete innerOperator_;
        delete innerScalarProduct_;
        delete innerPreCond_;
    }

    void pre(domain_type &x, range_type &y)
    {}

    void apply(domain_type &x, const range_type &d)
    {
        range_type dd(d);
        Dune::InverseOperatorResult result;
        innerSolver_->apply(x, dd, result);
    }

    void post(domain_type &x)
    {}

private:
    InnerOperator *innerOperator_;
    InnerScalarProduct *innerScalarProduct_;
    InnerPreConditioner *innerPreCond_;
    InnerSolver *innerSolver_;
};

} // namespace Linear
} // namespace Dumux

#endif
