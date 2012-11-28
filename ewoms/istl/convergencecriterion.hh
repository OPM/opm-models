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
 * \copydoc Ewoms::ConvergenceCriterion
 */
#ifndef EWOMS_ISTL_CONVERGENCE_CRITERION_HH
#define EWOMS_ISTL_CONVERGENCE_CRITERION_HH

#include <cmath>
#include <iostream>
#include <iomanip>

namespace Ewoms {
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

  /*!
   * \brief Prints the initial information about the convergence behaviour.
   *
   * This method is called after setInitial() if the solver thinks
   * it's a good idea to be verbose. In practice, "printing the
   * initial information" means printing column headers and the
   * initial state.
   *
   * \param os The output stream to which the message gets written.
   */
  virtual void printInitial(std::ostream &os=std::cout) const
  { };

  /*!
   * \brief Prints the information about the convergence behaviour for
   *        the current iteration.
   *
   * \param iter The iteration number. The semantics of this parameter
   *             are chosen by the linear solver.
   * \param os The output stream to which the message gets written.
   */
  virtual void print(Scalar iter, std::ostream &os=std::cout) const
  { };
};

//! \} end documentation

} // end namespace Ewoms

#endif
