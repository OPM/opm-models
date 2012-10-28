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
 * \copydoc Dumux::WeightedResidReductionCriterion
 */
#ifndef DUMUX_ISTL_WEIGHTED_RESID_REDUCTION_CRITERION_HH
#define DUMUX_ISTL_WEIGHTED_RESID_REDUCTION_CRITERION_HH

#include "convergencecriterion.hh"

#include <dune/common/mpihelper.hh>

namespace Dumux {
/*! \addtogroup ISTL_Solvers
 * \{
 */

/*!
 * \brief Convergence criterion which looks at the weighted absolute
 *        value of the residual
 *
 * For the WeightedResidReductionCriterion, the error of the solution is defined
 * as
 * \f[ e^k = \max_i\{ \left| w_i r^k_i \right| \}\;, \f]
 *
 * where \f$r^k = \mathbf{A} x^k - b \f$ is the residual for the
 * k-th iterative solution vector \f$x^k\f$ and \f$w_i\f$ is the
 * weight of the \f$i\f$-th linear equation.
 *
 * In addition, to the reduction of the maximum defect, the linear
 * solver is also considered to be converged, if the defect goes below
 * a given absolute limit.
 */
template <class Vector, class CollectiveCommunication>
class WeightedResidReductionCriterion : public ConvergenceCriterion<Vector>
{
  typedef typename Vector::field_type Scalar;
  typedef typename Vector::block_type BlockType;

public:
  WeightedResidReductionCriterion(const CollectiveCommunication &comm)
    : comm_(comm)
  { }

  WeightedResidReductionCriterion(const CollectiveCommunication &comm,
                                  const Vector &weights,
                                  Scalar reduction,
                                  Scalar absTolerance = 0.0)
    : comm_(comm)
    , weightVec_(weights)
    , tolerance_(reduction)
    , absTolerance_(absTolerance)
  { }

  /*!
   * \brief Sets the relative weight of each equation.
   *
   * For the WeightedResidReductionCriterion, the error of the solution is defined
   * as
   * \f[ e^k = \max_i\{ \left| w_i r^k_i \right| \}\;, \f]
   *
   * where \f$r^k = \mathbf{A} x^k - b \f$ is the residual for the
   * k-th iterative solution vector \f$x^k\f$ and \f$w_i\f$ is the
   * weight of the \f$i\f$-th linear equation.
   *
   * This method is not part of the generic ConvergenceCriteria interface.
   *
   * \param weightVec A Dune::BlockVector<Dune::FieldVector<Scalar, n> >
   *                  with the relative weights of the linear equations
   */
  void setWeight(const Vector &weightVec)
  {
    weightVec_ = weightVec;
  }

  /*!
   * \brief Return the relative weight of a primary variable
   *
   * For the FixPointCriterion, the error of the solution is defined
   * as
   * \f[ e^k = \max_i\{ \left| w_i \delta^k_i \right| \}\;, \f]
   *
   * where \f$\delta = x^k - x^{k + 1} \f$ is the difference between
   * two consequtive iterative solution vectors \f$x^k\f$ and \f$x^{k + 1}\f$
   * and \f$w_i\f$ is the weight of the \f$i\f$-th degree of freedom.
   *
   * This method is specific to the FixPointCriterion.
   *
   * \param outerIdx The index of the outer vector (i.e. Dune::BlockVector)
   * \param innerIdx The index of the inner vector (i.e. Dune::FieldVector)
   */
  Scalar weight(int outerIdx, int innerIdx) const
  {
    return (weightVec_.size() == 0)?1.0:weightVec_[outerIdx][innerIdx];
  }

  /*!
   * \copydoc ConvergenceCriterion::setTolerance(Scalar)
   */
  void setTolerance(Scalar tol)
  {
    tolerance_ = tol;
  }

  /*!
   * \brief Returns true if and only if the convergence criterion is
   *        met.
   */
  Scalar tolerance() const
  {
    return tolerance_;
  }

  /*!
   * \copydoc ConvergenceCriterion::accuracy()
   */
  Scalar accuracy() const
  {
    return error_/initialError_;
  }

  /*!
   * \copydoc ConvergenceCriterion::setInitial(const Vector &, const Vector &)
   */
  void setInitial(const Vector &curSol,
                  const Vector &curResid)
  {
    lastError_ = 1e100;
    updateError_(curSol, curResid);

    // make sure that we don't allow an initial error of 0 to avoid
    // divisions by zero
    error_ = std::max(error_, 1e-20);
    initialError_ = error_;
  }

  /*!
   * \copydoc ConvergenceCriterion::update(const Vector &, const Vector &)
   */
  void update(const Vector &curSol,
              const Vector &curResid)
  {
    lastError_ = error_;
    updateError_(curSol, curResid);
  }

  /*!
   * \copydoc ConvergenceCriterion::converged()
   */
  bool converged() const
  {
    // we're converged if difference between two iterations is either
    // smaller than the numerical precision, or if the solution is
    // better than the tolerance. TODO: are we actually converged in
    // the first case?
    return
      accuracy() <= tolerance() ||
      error_ <= absTolerance_;
  }

  /*!
   * \copydoc ConvergenceCriterion::printInitial()
   */
  void printInitial(std::ostream &os=std::cout) const
  {
    os << std::setw(20) << " Iter ";
    os << std::setw(20) << " Accuracy ";
    os << std::setw(20) << " Defect ";
    os << std::setw(20) << " Rate " << std::endl;

    os << std::setw(20) << 0 << " ";
    os << std::setw(20) << accuracy() << " ";
    os << std::setw(20) << error_ << " ";
    os << std::setw(20) << 0 << " ";
    os << std::endl;
  }

  /*!
   * \copydoc ConvergenceCriterion::print()
   */
  void print(Scalar iter, std::ostream &os=std::cout) const
  {
    os << std::setw(20) << iter << " ";
    os << std::setw(20) << accuracy() << " ";
    os << std::setw(20) << error_ << " ";
    os << std::setw(20) << lastError_/std::max(error_, 1e-80) << " ";
    os << std::endl;
  };

private:
  // update the weighted absolute residual
  void updateError_(const Vector &curSol, const Vector &curResid)
  {
    error_ = 0.0;
    for (size_t i = 0; i < curResid.size(); ++i) {
      for (size_t j = 0; j < BlockType::dimension; ++j) {
        error_ = std::max<Scalar>(error_, weight(i, j)*std::abs(curResid[i][j]));
      }
    }

    error_ = comm_.max(error_);
  }

  const CollectiveCommunication &comm_;

  Vector weightVec_; // the weights of the components of the residual

  Scalar error_; // the maximum of the absolute weighted difference of the last two iterations
  Scalar lastError_; // the maximum of the absolute weighted difference of the last two iterations
  Scalar initialError_; // the maximum of the absolute weighted difference of the last two iterations
  Scalar tolerance_; // the maximum allowed relative tolerance for the solution to be considered converged
  Scalar absTolerance_; // the maximum allowed absolute tolerance for the solution to be considered converged
};

//! \} end documentation

} // end namespace Dumux

#endif

