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
 * \copydoc Ewoms::Linear::ParallelIstlSolverBackend
 */
#ifndef EWOMS_PARALLEL_ISTL_BACKEND_HH
#define EWOMS_PARALLEL_ISTL_BACKEND_HH

#include "parallelbasebackend.hh"
#include "istlsolverwrappers.hh"

#include <dune/common/version.hh>

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(ParallelIstlLinearSolver, INHERITS_FROM(ParallelBaseLinearSolver));

NEW_PROP_TAG(LinearSolverWrapper);

//! number of iterations between solver restarts for the GMRES solver
NEW_PROP_TAG(GMResRestart);
} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides all unmodified linear solvers from dune-istl
 *
 * To set the linear solver, use
 * \code
 * SET_TYPE_PROP(YourTypeTag, LinearSolverWrapper,
 *               Ewoms::Linear::SolverWrapper$SOLVER<TypeTag>);
 * \endcode
 *
 * The possible choices for '\c $SOLVER' are:
 * - \c Richardson: A fixpoint solver using the Richardson iteration
 * - \c SteepestDescent: The steepest descent solver
 * - \c ConjugatedGradients: A conjugated gradients solver
 * - \c BiCGStab: A stabilized bi-conjugated gradients solver
 * - \c MinRes: A solver based on the  minimized residual algorithm
 * - \c RestartedGMRes: A restarted GMRES solver
 *
 * Chosing the preconditioner works in an analogous way:
 * \code
 * SET_TYPE_PROP(YourTypeTag, PreconditionerWrapper,
 *               Ewoms::Linear::PreconditionerWrapper$PRECONDITIONER<TypeTag>);
 * \endcode
 *
 * Where the choices possible for '\c $PRECONDITIONER' are:
 * - \c Jacobi: A Jacobi preconditioner
 * - \c GaussSeidel: A Gauss-Seidel preconditioner
 * - \c SSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - \c SOR: A successive overrelaxation (SOR) preconditioner
 * - \c ILUn: An ILU(n) preconditioner
 * - \c ILU0: A specialized (and optimized) ILU(0) preconditioner
 */
template <class TypeTag>
class ParallelIstlSolverBackend : public ParallelBaseBackend<TypeTag>
{
    typedef ParallelBaseBackend<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverWrapper) LinearSolverWrapper;

    typedef typename ParentType::ParallelOperator ParallelOperator;
    typedef typename ParentType::OverlappingVector OverlappingVector;
    typedef typename ParentType::ParallelPreconditioner ParallelPreconditioner;
    typedef typename ParentType::ParallelScalarProduct ParallelScalarProduct;

    typedef typename LinearSolverWrapper::RawSolver RawLinearSolver;

    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;

public:
    ParallelIstlSolverBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        LinearSolverWrapper::registerParameters();
    }

protected:
    friend ParentType;

    std::shared_ptr<RawLinearSolver> prepareSolver_(ParallelOperator& parOperator,
                                                    ParallelScalarProduct& parScalarProduct,
                                                    ParallelPreconditioner& parPreCond)
    {
        return solverWrapper_.get(parOperator,
                                  parScalarProduct,
                                  parPreCond);
    }

    void cleanupSolver_()
    {
        solverWrapper_.cleanup();
    }

    bool runSolver_(std::shared_ptr<RawLinearSolver> solver)
    {
        Dune::InverseOperatorResult result;
        solver->apply(*this->overlappingx_, *this->overlappingb_, result);
        return result.converged;
    }

    LinearSolverWrapper solverWrapper_;
};

}} // namespace Linear, Ewoms

namespace Ewoms {
namespace Properties {
SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverBackend,
              Ewoms::Linear::ParallelIstlSolverBackend<TypeTag>);

SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverWrapper,
              Ewoms::Linear::SolverWrapperBiCGStab<TypeTag>);

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
SET_TYPE_PROP(ParallelIstlLinearSolver,
              PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILU<TypeTag>);
#else
SET_TYPE_PROP(ParallelIstlLinearSolver,
              PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILU0<TypeTag>);
#endif

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(ParallelIstlLinearSolver, GMResRestart, 10);
}} // namespace Properties, Ewoms

#endif
