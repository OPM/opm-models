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
 * \copydoc Ewoms::WeightedResidReductionCriterion
 */
#ifndef EWOMS_ISTL_WEIGHTED_RESID_REDUCTION_CRITERION_HH
#define EWOMS_ISTL_WEIGHTED_RESID_REDUCTION_CRITERION_HH

#include "convergencecriterion.hh"

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <iostream>

namespace Ewoms {
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
                                  const Vector &fixPointWeights,
                                  const Vector &residWeights,
                                  Scalar fixPointTolerance,
                                  Scalar residTolerance,
                                  Scalar absResidTolerance = 0.0)
    : comm_(comm)
    , fixPointWeightVec_(fixPointWeights)
    , residWeightVec_(residWeights)
    , fixPointTolerance_(fixPointTolerance)
    , residTolerance_(residTolerance)
    , absResidTolerance_(absResidTolerance)
  { }

  /*!
   * \brief Sets the relative weight of each row of the residual.
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
   * \param residWeightVec A Dune::BlockVector<Dune::FieldVector<Scalar, n> >
   *                  with the relative weights of the linear equations
   */
  void setResidualWeight(const Vector &residWeightVec)
  {
    residWeightVec_ = residWeightVec;
  }

  /*!
   * \brief Sets the relative weight of each row of the solution.
   *
   * \param fixPointWeightVec A Dune::BlockVector<Dune::FieldVector<Scalar, n> >
   *                 with the relative weights of the linear equations
   */
  void setFixPointWeight(const Vector &fixPointWeightVec)
  {
    fixPointWeightVec_ = fixPointWeightVec;
  }

  /*!
   * \brief Return the relative weight of a row of the residual.
   *
   * \param outerIdx The index of the outer vector (i.e. Dune::BlockVector)
   * \param innerIdx The index of the inner vector (i.e. Dune::FieldVector)
   */
  Scalar residualWeight(int outerIdx, int innerIdx) const
  {
    return (residWeightVec_.size() == 0)?1.0:residWeightVec_[outerIdx][innerIdx];
  }

  /*!
   * \brief Return the relative weight of a row of the solution vector.
   *
   * \param outerIdx The index of the outer vector (i.e. Dune::BlockVector)
   * \param innerIdx The index of the inner vector (i.e. Dune::FieldVector)
   */
  Scalar fixPointWeight(int outerIdx, int innerIdx) const
  {
    return (fixPointWeightVec_.size() == 0)?1.0:fixPointWeightVec_[outerIdx][innerIdx];
  }

  /*!
   * \brief Sets the residual reduction tolerance.
   */
  void setResidReductionTolerance(Scalar tol)
  {
    residTolerance_ = tol;
  }

  /*!
   * \brief Returns the tolerance of the residual reduction of the solution.
   */
  Scalar residReductionTolerance() const
  {
    return residTolerance_;
  }

  /*!
   * \brief Sets the maximum absolute tolerated residual.
   */
  void setResidTolerance(Scalar tol)
  {
    absResidTolerance_ = tol;
  }

  /*!
   * \brief Returns the maximum absolute tolerated residual.
   */
  Scalar absResidTolerance() const
  {
    return absResidTolerance_;
  }

  /*!
   * \brief Returns the reduction of the weighted maximum of the
   *        residual compared to the initial solution.
   */
  Scalar residAccuracy() const
  {
    return residError_/initialResidError_;
  }

  /*!
   * \brief Sets the fix-point tolerance.
   */
  void setFixPointTolerance(Scalar tol)
  {
    fixPointTolerance_ = tol;
  }

  /*!
   * \brief Returns the maximum tolerated difference between two
   *        iterations to be met before a solution is considered to be
   *        converged.
   */
  Scalar fixPointTolerance() const
  {
    return fixPointTolerance_;
  }

  /*!
   * \brief Returns the weighted maximum of the difference
   *        between the last two iterative solutions.
   */
  Scalar fixPointAccuracy() const
  {
    return fixPointError_;
  }

  /*!
   * \copydoc ConvergenceCriterion::setInitial(const Vector &, const Vector &)
   */
  void setInitial(const Vector &curSol,
                  const Vector &curResid)
  {
    lastResidError_ = 1e100;

    lastSolVec_ = curSol;
    updateErrors_(curSol, curResid);
    // the fix-point error is not applicable for the initial solution!
    fixPointError_ = 1e100;

    // make sure that we don't allow an initial error of 0 to avoid
    // divisions by zero
    residError_ = std::max(residError_, 1e-20);
    initialResidError_ = residError_;
  }

  /*!
   * \copydoc ConvergenceCriterion::update(const Vector &, const Vector &)
   */
  void update(const Vector &curSol,
              const Vector &curResid)
  {
    lastResidError_ = residError_;
    updateErrors_(curSol, curResid);
  }

  /*!
   * \copydoc ConvergenceCriterion::converged()
   */
  bool converged() const
  {
    // we're converged if the solution is better than the tolerance
    // fix-point and residual tolerance.
    return fixPointAccuracy() < fixPointTolerance() &&
           ( residAccuracy() <= residReductionTolerance() ||
             residError_ <= absResidTolerance_);
  }

  /*!
   * \copydoc ConvergenceCriterion::printInitial()
   */
  void printInitial(std::ostream &os=std::cout) const
  {
    os << std::setw(20) << " Iter ";
    os << std::setw(20) << " Delta ";
    os << std::setw(20) << " ResidRed ";
    os << std::setw(20) << " Defect ";
    os << std::setw(20) << " Rate " << std::endl;

    os << std::setw(20) << 0 << " ";
    os << std::setw(20) << 1e100 << " ";
    os << std::setw(20) << residAccuracy() << " ";
    os << std::setw(20) << residError_ << " ";
    os << std::setw(20) << 0 << " ";
    os << std::endl;
  }

  /*!
   * \copydoc ConvergenceCriterion::print()
   */
  void print(Scalar iter, std::ostream &os=std::cout) const
  {
    os << std::setw(20) << iter << " ";
    os << std::setw(20) << fixPointAccuracy() << " ";
    os << std::setw(20) << residAccuracy() << " ";
    os << std::setw(20) << residError_ << " ";
    os << std::setw(20) << lastResidError_/std::max(residError_, 1e-80) << " ";
    os << std::endl;
  };

private:
  // update the weighted absolute residual
  void updateErrors_(const Vector &curSol, const Vector &curResid)
  {
    residError_ = 0.0;
    fixPointError_ = 0.0;
    for (size_t i = 0; i < curResid.size(); ++i) {
      for (size_t j = 0; j < BlockType::dimension; ++j) {
        residError_ = std::max<Scalar>(residError_, residualWeight(i, j)*std::abs(curResid[i][j]));
        fixPointError_ = std::max<Scalar>(fixPointError_, fixPointWeight(i, j)*std::abs(curSol[i][j] - lastSolVec_[i][j]));
      }
    }
    lastSolVec_ = curSol;

    residError_ = comm_.max(residError_);
    fixPointError_ = comm_.max(fixPointError_);
  }

  const CollectiveCommunication &comm_;

  Vector fixPointWeightVec_; // the weights of the solution vector
  Vector residWeightVec_; // the weights of the components of the residual
  Vector lastSolVec_; // solution vector of the last iteration

  Scalar fixPointError_; // the maximum of the weighted difference between the last two iterations
  Scalar fixPointTolerance_; // the maximum allowed relative tolerance for difference of the solution of two iterations

  Scalar residError_; // the maximum of the absolute weighted residual of the last iteration
  Scalar lastResidError_; // the maximum of the absolute weighted difference of the last iteration
  Scalar initialResidError_; // the maximum of the absolute weighted residual of the initial solution
  Scalar residTolerance_; // the maximum allowed relative tolerance of the residual for the solution to be considered converged
  Scalar absResidTolerance_; // the maximum allowed absolute tolerance of the residual for the solution to be considered converged
};

//! \} end documentation

} // end namespace Ewoms

#endif
