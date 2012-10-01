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
 * \copydoc Dumux::ResidReductionCriterion
 */
#ifndef DUMUX_ISTL_RESID_REDUCTION_CRITERION_HH
#define DUMUX_ISTL_RESID_REDUCTION_CRITERION_HH

#include "convergencecriterion.hh"

#include <dune/istl/scalarproducts.hh>

#include <dune/common/mpihelper.hh>

namespace Dumux {
/*! \addtogroup ISTL_Solvers
 * \{
 */

/*!
 * \brief Provides a convergence criterion which looks at the
 *        reduction of the two-norm of the residual for the linear
 *        solvers.
 *
 * For the ResidReductionCriterion, the error of the solution is defined 
 * as
 * \f[ e^k = \frac{\left| A x_k - b \right|}{\left| A x_0 - b \right|}\;, \f]
 */
template<class Vector>
class ResidReductionCriterion : public ConvergenceCriterion<Vector>
{
  typedef typename Vector::field_type Scalar;
  
public:
  ResidReductionCriterion(Dune::ScalarProduct<Vector> &scalarProduct)
    : scalarProduct_(scalarProduct)
  {}

  ResidReductionCriterion(Dune::ScalarProduct<Vector> &scalarProduct, Scalar reduction)
    : scalarProduct_(scalarProduct)
    , defectReduction_(reduction)
  {}

  /*!
   * \copydoc ConvergenceCriterion::setTolerance(Scalar)
   */
  void setTolerance(Scalar tol)
  {
    defectReduction_ = tol;
  }

  /*!
   * \copydoc ConvergenceCriterion::tolerance()
   */
  Scalar tolerance() const
  { 
    return defectReduction_;
  }

  /*!
   * \copydoc ConvergenceCriterion::accuracy()
   */
  Scalar accuracy() const
  {
    return curDefect_/initialDefect_;
  }

  /*!
   * \copydoc ConvergenceCriterion::setInitial(const Vector &, const Vector &)
   */
  void setInitial(const Vector &curSol,
                  const Vector &curResid)
  {
    // make sure that we don't allow an initial error of 0 to avoid
    // divisions by zero
    initialDefect_ = std::max(scalarProduct_.norm(curResid), 1e-20);
    curDefect_ = initialDefect_;
  }

  /*!
   * \copydoc ConvergenceCriterion::update(const Vector &, const Vector &)
   */
  void update(const Vector &curSol, 
              const Vector &curResid)
  {
    curDefect_ = scalarProduct_.norm(curResid);
  }

  /*!
   * \copydoc ConvergenceCriterion::converged()
   */
  bool converged() const
  {
    return accuracy() <= tolerance();
  }

private:
  Dune::ScalarProduct<Vector> &scalarProduct_;

  Scalar initialDefect_;
  Scalar curDefect_;
  Scalar defectReduction_;
};

//! \} end documentation

} // end namespace Dumux

#endif

