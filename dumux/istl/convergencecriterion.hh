// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::ConvergenceCriterion
 */
#ifndef DUMUX_ISTL_CONVERGENCE_CRITERION_HH
#define DUMUX_ISTL_CONVERGENCE_CRITERION_HH

#include <cmath>

namespace Dumux {
/*! \addtogroup ISTL_Solvers
 * \{
 */

/*!
 * \file
 * \brief Define some base class for the convergence criteria of the linear solvers of DUNE-ISTL.
 */

/*!
 * \brief Base class for all convergence criteria which only defines an virtual API.
 */
template <class Vector>
class ConvergenceCriterion
{
  typedef typename Vector::field_type Scalar;
public:
  /*!
   * \brief Destructor.
   *
   * In the ConvergenceCriterion it does not do anything, but it is
   * required to be declared virtual.
   */
  virtual ~ConvergenceCriterion()
  {}

  /*!
   * \brief Set the maximum acceptable accuracy the solution for
   *        which the iterative solution ought to be considered to be
   *        converged.
   *
   * \param tol The maximum acceptable error for convergence of the
   *            iterative solution
   */
  virtual void setTolerance(Scalar tol) = 0;

  /*!
   * \brief Return the maximum acceptable accuracy the solution for
   *        which the iterative solution ought to be considered to be
   *        converged.
   */
  virtual Scalar tolerance() const = 0;

  /*!
   * \brief Return the accuracy of the solution given as an argument
   *        to the last call of the update() method.
   */
  virtual Scalar accuracy() const = 0;

  /*!
   * \brief Set the initial solution of the linear system of equations.
   *
   * This version of the method does NOT take the two-norm of the
   * residual as argument. If the two-norm of the defect is available
   * for the linear solver, the version of the update() method with it
   * should be called.
   *
   * \param curSol The current iterative solution of the linear system
   *               of equations
   * \param curResid The residual vector of the current iterative
   *                 solution of the linear system of equations
   */
  virtual void setInitial(const Vector &curSol,
                          const Vector &curResid) = 0;

  /*!
   * \brief Update the internal members of the convergence criterion
   *        with the current solution.
   *
   * This version of the method does NOT take the two-norm of the
   * residual as argument. If the two-norm of the defect is available
   * for the linear solver, the version of the update() method with it
   * should be called.
   *
   * \param curSol The current iterative solution of the linear system
   *               of equations
   * \param curResid The residual vector of the current iterative
   *                 solution of the linear system of equations
   */
  virtual void update(const Vector &curSol, 
                      const Vector &curResid) = 0;

  /*!
   * \brief Returns true if and only if the convergence criterion is
   *        met.
   */
  virtual bool converged() const = 0;
};

//! \} end documentation

} // end namespace Dumux

#endif
