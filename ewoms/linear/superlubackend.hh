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

private:
    const Problem &problem_;
};

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
