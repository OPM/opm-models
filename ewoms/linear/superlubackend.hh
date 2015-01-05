/*
  Copyright (C) 2012-2013 by Andreas Lauser
  Copyright (C) 2011 by Markus Wolff

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
 * \copydoc Ewoms::Linear::SuperLUBackend
 */
#ifndef EWOMS_SUPER_LU_BACKEND_HH
#define EWOMS_SUPER_LU_BACKEND_HH

#if HAVE_SUPERLU

#include <ewoms/common/parametersystem.hh>

#include <ewoms/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {
namespace Properties {
// forward declaration of the required property tags
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverBackend);

NEW_TYPE_TAG(SuperLULinearSolver);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
namespace Linear {
template <class Scalar, class TypeTag, class Problem, class Matrix, class Vector>
class SuperLUSolve_;

/*!
 * \ingroup Linear
 * \brief A linear solver backend for the SuperLU sparse matrix library.
 */
template <class TypeTag>
class SuperLUBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    SuperLUBackend(const Problem &problem) : problem_(problem)
    {}

    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");
    }

    template <class Matrix, class Vector>
    bool solve(const Matrix &A, Vector &x, const Vector &b)
    { return SuperLUSolve_<Scalar, TypeTag, Problem, Matrix, Vector>::solve_(problem_, A, x, b); }

private:
    const Problem &problem_;
};

template <class Scalar, class TypeTag, class Problem, class Matrix, class Vector>
class SuperLUSolve_
{
public:
    static bool solve_(const Problem &problem, const Matrix &A, Vector &x, const Vector &b)
    {
        Vector bTmp(b);

        int verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        Dune::InverseOperatorResult result;
        Dune::SuperLU<Matrix> solver(A, verbosity > 0);
        solver.apply(x, bTmp, result);

        if (result.converged) {
            // make sure that the result only contains finite values.
            Scalar tmp = 0;
            for (unsigned i = 0; i < x.size(); ++i) {
                const auto &xi = x[i];
                for (int j = 0; j < Vector::block_type::dimension; ++j)
                    tmp += xi[j];
            }
            result.converged = std::isfinite(tmp);
        }

        return result.converged;
    }
};

// the following is required to make the SuperLU adapter of dune-istl happy with
// quadruple precision math on Dune 2.4. this is because the most which SuperLU can
// handle is double precision (i.e., the linear systems of equations are always solved
// with at most double precision if chosing SuperLU as the linear solver...)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4) && HAVE_QUAD
template <class TypeTag, class Problem, class Matrix, class Vector>
class SuperLUSolve_<__float128, TypeTag, Problem, Matrix, Vector>
{
public:
    static bool solve_(const Problem &problem,
                       const Matrix &A,
                       Vector &x,
                       const Vector &b)
    {
        static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);
        typedef Dune::FieldVector<double, numEq> DoubleEqVector;
        typedef Dune::FieldMatrix<double, numEq, numEq> DoubleEqMatrix;
        typedef Dune::BlockVector<DoubleEqVector> DoubleVector;
        typedef Dune::BCRSMatrix<DoubleEqMatrix> DoubleMatrix;

        // copy the inputs into the double precision data structures
        DoubleVector bDouble(b);
        DoubleVector xDouble(x);
        DoubleMatrix ADouble(A);

        bool res =
            SuperLUSolve_<double, TypeTag, Problem, Matrix, Vector>::solve_(problem,
                                                                            ADouble,
                                                                            xDouble,
                                                                            bDouble);

        // copy the result back into the quadruple precision vector.
        x = xDouble;

        return res;
    }
};
#endif

} // namespace Linear
} // namespace Ewoms

namespace Opm {
namespace Properties {
SET_INT_PROP(SuperLULinearSolver, LinearSolverVerbosity, 0);
SET_TYPE_PROP(SuperLULinearSolver, LinearSolverBackend,
              Ewoms::Linear::SuperLUBackend<TypeTag>);
} // namespace Properties
} // namespace Opm

#endif // HAVE_SUPERLU

#endif
