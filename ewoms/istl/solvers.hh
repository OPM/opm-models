/*
  Copyright (C) 2004-2012 by Christian Engwer
  Copyright (C) 2005-2011 by Markus Blatt
  Copyright (C) 2007 by Robert Kloefkorn
  Copyright (C) 2004-2006 by Peter Bastian
  Copyright (C) 2010 by Jorrit Fahlke
  Copyright (C) 2007-2012 by Oliver Sander
  Copyright (C) 2012 by Andreas Lauser

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, version 2.

  As a special exception, you may use the DUNE library without
  restriction.  Specifically, if other files instantiate templates or
  use macros or inline functions from one or more of the DUNE source
  files, or you compile one or more of the DUNE source files and link
  them with other files to produce an executable, this does not by
  itself cause the resulting executable to be covered by the GNU
  General Public License.  This exception does not however invalidate
  any other reasons why the executable file might be covered by the
  GNU General Public License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \brief Copy of dune-istl's linear solvers with added support for
 *        pluggable convergence criteria.
 *
 * For eWoms, pluggable convergence criteria for the linear solvers
 * are an important feature. Unfortunatly, the DUNE developers don't
 * seem to care, so this could not go directly into ISTL. For
 * details, see
 *
 * http://www.dune-project.org/flyspray/index.php?do=details&task_id=1018
 */
#ifndef EWOMS_SOLVERS_HH
#define EWOMS_SOLVERS_HH

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "residreductioncriterion.hh"
#include "weightedresidreductioncriterion.hh"
#include "fixpointcriterion.hh"

#include <dune/istl/istlexception.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/common/version.hh>
#include <dune/common/array.hh>
#include <dune/common/timer.hh>
#include <dune/common/ftraits.hh>

#include <type_traits>
#include <algorithm>
#include <memory>

namespace Ewoms {
/*
  @ingroup Linear
*/
/** @addtogroup Linear
    \{
*/

//=====================================================================
/*!
  \brief Abstract base class for all solvers.

  An InverseOperator computes the solution of \f$ A(x)=b\f$ where
  \f$ A : X \to Y \f$ is an operator.
  Note that the solver "knows" which operator
  to invert and which preconditioner to apply (if any). The
  user is only interested in inverting the operator.
  InverseOperator might be a Newton scheme, a Krylov subspace method,
  or a direct solver or just anything.
*/
template <class X, class Y>
class InverseOperator
{
public:
    //! \brief Type of the domain of the operator to be inverted.
    typedef X domain_type;

    //! \brief Type of the range of the operator to be inverted.
    typedef Y range_type;

    /** \brief The field type of the operator. */
    typedef typename X::field_type field_type;

    /*!
     * \brief Return the criterion to be used to check for convergence of the
     * linear solver.
     */
    virtual const Ewoms::ConvergenceCriterion<X> &convergenceCriterion() const
    { return *convergenceCriterion_; }

    /*!
     * \brief Return the criterion to be used to check for convergence of the
     * linear solver.
     */
    virtual Ewoms::ConvergenceCriterion<X> &convergenceCriterion()
    { return *convergenceCriterion_; }

    /*!
     * \brief Set the criterion to be used to check for convergence of the
     * linear solver.
     */
    virtual void setConvergenceCriterion(std::shared_ptr<Ewoms::ConvergenceCriterion<X> > convCrit)
    { convergenceCriterion_ = convCrit; }

    /*!
      \brief Apply inverse operator,

      \warning Note: right hand side b may be overwritten!

      \param x The left hand side to store the result in.
      \param b The right hand side
      \param res Object to store the statistics about applying the operator.
    */
    virtual void apply(X& x, Y& b, Dune::InverseOperatorResult& res) = 0;

    //! \brief Destructor
    virtual ~InverseOperator()
    {}

private:
    std::shared_ptr<ConvergenceCriterion<X> > convergenceCriterion_;
};

//=====================================================================
// Implementation of this interface
//=====================================================================

/*!
  \brief Preconditioned loop solver.

  Implements a preconditioned loop.
  Using this class every Preconditioner can be turned
  into a solver. The solver will apply one preconditioner
  step in each iteration loop.
*/
template <class X>
class LoopSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator that we do the inverse for.
    typedef X domain_type;
    //! \brief The range type of the operator that we do the inverse for.
    typedef X range_type;
    //! \brief The field type of the operator that we do the inverse for.
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif

    /*!
      \brief Set up Loop solver.

      \param op The operator we solve.
      \param prec The preconditioner to apply in each iteration of the loop.
      Has to inherit from Preconditioner.
      \param reduction The relative defect reduction to achieve when applying
      the operator.
      \param maxit The maximum number of iteration steps allowed when applying
      the operator.
      \param verbose The verbosity level.

      Verbose levels are:
      <ul>
      <li> 0 : print nothing </li>
      <li> 1 : print initial and final defect and statistics </li>
      <li> 2 : print line for each iteration </li>
      </ul>
    */
    template <class L, class P>
    LoopSolver(L& op, P& prec,
               real_type reduction, int maxit, int verbose) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P have to have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L has to be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /**
       \brief Set up loop solver

       \param op The operator we solve.
       \param sp The scalar product to use, e. g. SeqScalarproduct.
       \param prec The preconditioner to apply in each iteration of the loop.
       Has to inherit from Preconditioner.
       \param reduction The relative defect reduction to achieve when applying
       the operator.
       \param maxit The maximum number of iteration steps allowed when applying
       the operator.
       \param verbose The verbosity level.

       Verbose levels are:
       <ul>
       <li> 0 : print nothing </li>
       <li> 1 : print initial and final defect and statistics </li>
       <li> 2 : print line for each iteration </li>
       </ul>
    */
    template <class L, class S, class P>
    LoopSolver(L& op, S& sp, P& prec,
               real_type reduction, int maxit, int verbose) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }


    //! \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        // clear solver statistics
        res.clear();

        // start a timer
        Dune::Timer watch;

        // prepare preconditioner
        _prec.pre(x, b);

        // overwrite b with defect
        _op.applyscaleadd(-1, x, b);

        this->convergenceCriterion().setInitial(x, b);
        if (_verbose > 0) {
            std::cout << "=== LoopSolver" << std::endl << std::flush;
            if (_verbose > 1) {
                this->convergenceCriterion().printInitial();
            }
        }
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        // allocate correction vector
        X v(x);

        // iteration loop
        int i = 1;
        for (; i <= _maxit; i++) {
            v = 0; // clear correction
            _prec.apply(v, b); // apply preconditioner
            x += v; // update solution
            _op.applyscaleadd(-1, v, b); // update defect

            this->convergenceCriterion().update(x, b);
            if (_verbose > 1) // print
                this->convergenceCriterion().print(i);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }
        }

        //correct i which is wrong if convergence was not achieved.
        i = std::min(_maxit, i);

        // print
        if (_verbose == 1)
            this->convergenceCriterion().print(i);

        // postprocess preconditioner
        _prec.post(x);

        // fill statistics
        res.iterations = i;
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();

        // final print
        if (_verbose > 0)
        {
            std::cout << "=== rate=" << res.conv_rate
                      << ", T=" << res.elapsed
                      << ", TIT=" << res.elapsed/i
                      << ", IT=" << i << std::endl;
        }
    }

private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X> &_op;
    Dune::Preconditioner<X, X> &_prec;
    Dune::ScalarProduct<X> &_sp;
    std::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
};

// all these solvers are taken from the SUMO library
//! gradient method
template <class X>
class GradientSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator that we do the inverse for.
    typedef X domain_type;
    //! \brief The range type of the operator  that we do the inverse for.
    typedef X range_type;
    //! \brief The field type of the operator  that we do the inverse for.
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif


    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
    */
    template <class L, class P>
    GradientSolver(L& op, P& prec,
                   real_type reduction, int maxit, int verbose) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P have to have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L has to be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }
    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
    */
    template <class L, class S, class P>
    GradientSolver(L& op, S& sp, P& prec,
                   real_type reduction, int maxit, int verbose) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P have to have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S have to have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        res.clear(); // clear solver statistics
        Dune::Timer watch; // start a timer
        _prec.pre(x, b); // prepare preconditioner
        _op.applyscaleadd(-1, x, b); // overwrite b with defect

        X p(x); // create local vectors
        X q(b);

        this->convergenceCriterion().setInitial(x, b);
        if (_verbose > 0) {
            std::cout << "=== GradientSolver" << std::endl << std::flush;
            if (_verbose > 1)
                this->convergenceCriterion().printInitial();
        }
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        int i = 1; // loop variables
        field_type lambda;
        for (; i <=_maxit; i++)
        {
            p = 0; // clear correction
            _prec.apply(p, b); // apply preconditioner
            _op.apply(p, q); // q=Ap
            lambda = _sp.dot(p, b)/_sp.dot(q, p); // minimization
            x.axpy(lambda, p); // update solution
            b.axpy(-lambda, q); // update defect

            this->convergenceCriterion().update(x, b);
            if (_verbose > 1) // print
                this->convergenceCriterion().print(i);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }
        }

        //correct i which is wrong if convergence was not achieved.
        i = std::min(_maxit, i);

        if (_verbose == 1) // printing for non verbose
            this->convergenceCriterion().print(i);

        _prec.post(x); // postprocess preconditioner
        res.iterations = i; // fill statistics
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();
    }

private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X> &_op;
    Dune::Preconditioner<X, X> &_prec;
    Dune::ScalarProduct<X> &_sp;
    std::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
};



//! \brief conjugate gradient method
template <class X>
class CGSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;

    /*!
      \brief Set up conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
    */
    template <class L, class P>
    CGSolver(L& op, P& prec, real_type reduction, int maxit, int verbose) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L must be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }
    /*!
      \brief Set up conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
    */
    template <class L, class S, class P>
    CGSolver(L& op, S& sp, P& prec, real_type reduction, int maxit, int verbose) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        res.clear(); // clear solver statistics
        Dune::Timer watch; // start a timer
        _prec.pre(x, b); // prepare preconditioner
        _op.applyscaleadd(-1, x, b); // overwrite b with defect

        X p(x); // the search direction
        X q(x); // a temporary vector

        this->convergenceCriterion().setInitial(x, b);
        if (_verbose > 0)
        {
            std::cout << "=== CGSolver" << std::endl << std::flush;
            if (_verbose > 1)
                this->convergenceCriterion().printInitial();
        }
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        // some local variables
        field_type rho, rholast, lambda, alpha, beta;

        // determine initial search direction
        p = 0; // clear correction
        _prec.apply(p, b); // apply preconditioner
        rholast = _sp.dot(p, b); // orthogonalization

        // the loop
        int i = 1;
        for (; i <= _maxit; i++) {
            // minimize in given search direction p
            _op.apply(p, q); // q=Ap
            alpha = _sp.dot(p, q); // scalar product
            lambda = rholast/alpha; // minimization
            x.axpy(lambda, p); // update solution
            b.axpy(-lambda, q); // update defect

            // convergence test
            this->convergenceCriterion().update(x, b);
            if (_verbose > 1) // print
                this->convergenceCriterion().print(i);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }

            // determine new search direction
            q = 0; // clear correction
            _prec.apply(q, b); // apply preconditioner
            rho = _sp.dot(q, b); // orthogonalization
            beta = rho/rholast; // scaling factor
            p *= beta; // scale old search direction
            p += q; // orthogonalization with correction
            rholast = rho; // remember rho for recurrence
        }

        //correct i which is wrong if convergence was not achieved.
        i = std::min(_maxit, i);

        if (_verbose == 1) // printing for non verbose
            this->convergenceCriterion().print(i);

        _prec.post(x); // postprocess preconditioner
        res.iterations = i; // fill statistics
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();
    }

private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X> &_op;
    Dune::Preconditioner<X, X> &_prec;
    Dune::ScalarProduct<X> &_sp;
    std::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
};


// Ronald Kriemanns BiCG-STAB implementation from Sumo
//! \brief Bi-conjugate Gradient Stabilized (BiCG-STAB)
template <class X>
class BiCGSTABSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
    */
    template <class L, class P>
    BiCGSTABSolver(L& op, P& prec,
                   real_type reduction, int maxit, int verbose) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must be of the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L must be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }
    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
    */
    template <class L, class S, class P>
    BiCGSTABSolver(L& op, S& sp, P& prec,
                   real_type reduction, int maxit, int verbose) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        const real_type EPSILON = 1e-80;
        double it;
        field_type rho, rho_new, alpha, beta, h, omega;

        //
        // get vectors and matrix
        //
        X& r=b;
        X p(x);
        X v(x);
        X t(x);
        X y(x);
        X rt(x);

        //
        // begin iteration
        //

        // r = r - Ax; rt = r
        res.clear(); // clear solver statistics
        Dune::Timer watch; // start a timer
        _prec.pre(x, r); // prepare preconditioner
        _op.applyscaleadd(-1, x, r); // overwrite b with defect

        rt = r;

        p = 0;
        v = 0;

        rho = 1;
        alpha = 1;
        omega = 1;

        this->convergenceCriterion().setInitial(x, r);
        if (_verbose > 0)
        {
            std::cout << "=== BiCGSTABSolver" << std::endl << std::flush;
            if (_verbose > 1)
                this->convergenceCriterion().printInitial();
        }
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        //
        // iteration
        //
        for (it = 0.5; it < _maxit; it += .5)
        {
            //
            // preprocess, set vecsizes etc.
            //

            // rho_new = < rt, r >
            rho_new = _sp.dot(rt, r);

            // look if breakdown occured
            if (std::abs(rho) <= EPSILON)
                DUNE_THROW(Dune::ISTLError, "breakdown in BiCGSTAB - rho "
                           << rho << " <= EPSILON " << EPSILON
                           << " after " << it << " iterations");
            if (std::abs(omega) <= EPSILON)
                DUNE_THROW(Dune::ISTLError, "breakdown in BiCGSTAB - omega "
                           << omega << " <= EPSILON " << EPSILON
                           << " after " << it << " iterations");

            if (it < 1)
                p = r;
            else {
                beta = (rho_new / rho) * (alpha / omega);
                p.axpy(-omega, v); // p = r + beta (p - omega*v)
                p *= beta;
                p += r;
            }

            // y = W^-1 * p
            y = 0;
            _prec.apply(y, p); // apply preconditioner

            // v = A * y
            _op.apply(y, v);

            // alpha = rho_new / < rt, v >
            h = _sp.dot(rt, v);

            if (std::abs(h) < EPSILON)
                DUNE_THROW(Dune::ISTLError, "h=0 in BiCGSTAB");

            alpha = rho_new / h;

            // apply first correction to x
            // x <- x + alpha y
            x.axpy(alpha, y);

            // r = r - alpha*v
            r.axpy(-alpha, v);

            //
            // test stop criteria
            //
            this->convergenceCriterion().update(x, r);
            if (_verbose > 1) // print
                this->convergenceCriterion().print(it);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }
            it += .5;

            // y = W^-1 * r
            y = 0;
            _prec.apply(y, r);

            // t = A * y
            _op.apply(y, t);

            // omega = < t, r > / < t, t >
            omega = _sp.dot(t, r) / _sp.dot(t, t);

            // apply second correction to x
            // x <- x + omega y
            x.axpy(omega, y);

            // r = s - omega*t (remember : r = s)
            r.axpy(-omega, t);

            rho = rho_new;

            //
            // test stop criteria
            //
            this->convergenceCriterion().update(x, r);
            if (_verbose > 1) // print
                this->convergenceCriterion().print(it);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }
        } // end for

        //correct i which is wrong if convergence was not achieved.
        it = std::min(static_cast<double>(_maxit), it);

        if (_verbose == 1) // printing for non verbose
            this->convergenceCriterion().print(it);

        _prec.post(x); // postprocess preconditioner
        res.iterations = static_cast<int>(std::ceil(it)); // fill statistics
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();
    }

private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X> &_op;
    Dune::Preconditioner<X, X> &_prec;
    Dune::ScalarProduct<X> &_sp;
    std::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
};

/*! \brief Minimal Residual Method (MINRES)

  Symmetrically Preconditioned MINRES as in A. Greenbaum, 'Iterative Methods for Solving Linear Systems', pp. 121
  Iterative solver for symmetric indefinite operators.
  Note that in order to ensure the (symmetrically) preconditioned system to remain symmetric, the preconditioner has to be spd.
*/
template <class X>
class MINRESSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif

    /*!
      \brief Set up MINRES solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
    */
    template <class L, class P>
    MINRESSolver(L& op, P& prec, real_type reduction, int maxit, int verbose) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L must be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }
    /*!
      \brief Set up MINRES solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
    */
    template <class L, class S, class P>
    MINRESSolver(L& op, S& sp, P& prec, real_type reduction, int maxit, int verbose) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        // clear solver statistics
        res.clear();
        // start a timer
        Dune::Timer watch;
        watch.reset();
        // prepare preconditioner
        _prec.pre(x, b);
        // overwrite rhs with defect
        _op.applyscaleadd(-1, x, b);

        this->convergenceCriterion().setInitial(x, b);
        if (_verbose > 0)
        {
            std::cout << "=== MINRESSolver" << std::endl << std::flush;
            if (_verbose > 1) {
                this->convergenceCriterion().printInitial();
            }
        }
        // check for convergence
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        // recurrence coefficients as computed in Lanczos algorithm
        field_type alpha, beta;
        // diagonal entries of givens rotation
        std::array<real_type, 2> c{{0.0, 0.0}};
        // off-diagonal entries of givens rotation
        std::array<field_type, 2> s{{0.0, 0.0}};

        // recurrence coefficients (column k of tridiag matrix T_k)
        std::array<field_type, 3> T{{0.0, 0.0, 0.0}};

        // the rhs vector of the min problem
        std::array<field_type, 2> xi{{1.0, 0.0}};

        // some temporary vectors
        X z(b), dummy(b);

        // initialize and clear correction
        z = 0.0;
        _prec.apply(z, b);

        // beta is real and positive in exact arithmetic
        // since it is the norm of the basis vectors (in unpreconditioned case)
        beta = std::sqrt(_sp.dot(z, b));
        field_type beta0 = beta;

        // the search directions
        std::array<X, 3> p{{b, b, b}};
        p[0] = 0.0;
        p[1] = 0.0;
        p[2] = 0.0;

        // orthonormal basis vectors (in unpreconditioned case)
        std::array<X, 3> q{{b, b, b}};
        q[0] = 0.0;
        q[1] *= 1.0/beta;
        q[2] = 0.0;

        z *= 1.0/beta;

        // the loop
        int i = 1;
        for (; i <= _maxit; i++) {

            dummy = z;
            int i1 = i%3,
                i0 = (i1+2)%3,
                i2 = (i1+1)%3;

            // symmetrically preconditioned Lanczos algorithm (see Greenbaum p.121)
            _op.apply(z, q[i2]); // q[i2] = Az
            q[i2].axpy(-beta, q[i0]);
            // alpha is real since it is the diagonal entry of the hermitian tridiagonal matrix
            // from the Lanczos Algorithm
            // so the order in the scalar product doesn't matter even for the complex case
            alpha = _sp.dot(z, q[i2]);
            q[i2].axpy(-alpha, q[i1]);

            z = 0.0;
            _prec.apply(z, q[i2]);

            // beta is real and positive in exact arithmetic
            // since it is the norm of the basis vectors (in unpreconditioned case)
            beta = std::sqrt(_sp.dot(q[i2], z));

            q[i2] *= 1.0/beta;
            z *= 1.0/beta;

            // QR Factorization of recurrence coefficient matrix
            // apply previous givens rotations to last column of T
            T[1] = T[2];
            if (i > 2) {
                T[0] = s[i%2]*T[1];
                T[1] = c[i%2]*T[1];
            }
            if (i > 1) {
                T[2] = c[(i+1)%2]*alpha - s[(i+1)%2]*T[1];
                T[1] = c[(i+1)%2]*T[1] + s[(i+1)%2]*alpha;
            }
            else
                T[2] = alpha;

            // update QR factorization
            generateGivensRotation(T[2], beta, c[i%2], s[i%2]);
            // to last column of T_k
            T[2] = c[i%2]*T[2] + s[i%2]*beta;
            // and to the rhs xi of the min problem
            xi[i%2] = -s[i%2]*xi[(i+1)%2];
            xi[(i+1)%2] *= c[i%2];

            // compute correction direction
            p[i2] = dummy;
            p[i2].axpy(-T[1], p[i1]);
            p[i2].axpy(-T[0], p[i0]);
            p[i2] *= 1.0/T[2];

            // apply correction/update solution
            x.axpy(beta0*xi[(i+1)%2], p[i2]);

            // remember beta_old
            T[2] = beta;

            // check for convergence
            this->convergenceCriterion().update(x, b);
            if (_verbose > 1)
                this->convergenceCriterion().print(i);
            if (this->convergenceCriterion().converged()) {
                res.converged = true;
                break;
            }
        } // end for

        if (_verbose == 1) // printing for non verbose
            this->convergenceCriterion().print(i);

        // postprocess preconditioner
        _prec.post(x);
        // fill statistics
        res.iterations = i;
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();
    }

private:

    void generateGivensRotation(field_type& dx, field_type& dy, real_type& cs, field_type& sn)
    {
        real_type norm_dx = std::abs(dx);
        real_type norm_dy = std::abs(dy);
        if (norm_dy < 1e-15) {
            cs = 1.0;
            sn = 0.0;
        } else if (norm_dx < 1e-15) {
            cs = 0.0;
            sn = 1.0;
        } else if (norm_dy > norm_dx) {
            real_type temp = norm_dx/norm_dy;
            cs = 1.0/std::sqrt(1.0 + temp*temp);
            sn = cs;
            cs *= temp;
            sn *= dx/norm_dx;
            // dy is real in exact arithmetic
            // so we don't need to conjugate here
            sn *= dy/norm_dy;
        } else {
            real_type temp = norm_dy/norm_dx;
            cs = 1.0/std::sqrt(1.0 + temp*temp);
            sn = cs;
            sn *= dy/dx;
            // dy and dx is real in exact arithmetic
            // so we don't have to conjugate both of them
        }
    }

    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X>& _op;
    Dune::Preconditioner<X, X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    std::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
};

/**
   \brief implements the Generalized Minimal Residual (GMRes) method

   GMRes solves the unsymmetric linear system Ax = b using the
   Generalized Minimal Residual method as described the SIAM Templates
   book (http://www.netlib.org/templates/templates.pdf).

   \tparam X trial vector, vector type of the solution
   \tparam Y test vector, vector type of the RHS
   \tparam F vector type for orthonormal basis of Krylov space

*/

template <class X, class Y=X, class F = Y>
class RestartedGMResSolver : public InverseOperator<X, Y>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef Y range_type;
    //! \brief The field type of the operator to be inverted
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif
    //! \brief The field type of the basis vectors
    typedef F basis_type;

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
      \param restart number of GMRes cycles before restart
    */
    template <class L, class P>
    RestartedGMResSolver(L& op, P& prec, real_type reduction, int restart, int maxit, int verbose) :
        _A(op), _W(prec),
        ssp(), _sp(ssp), _restart(restart),
        _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(P::category) == static_cast<int>(L::category),
                      "P and L must be the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
                      "L must be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
      \param restart number of GMRes cycles before restart
    */
    template <class L, class S, class P>
    RestartedGMResSolver(L& op, S& sp, P& prec, real_type reduction, int restart, int maxit, int verbose) :
        _A(op), _W(prec),
        _sp(sp), _restart(restart),
        _maxit(maxit), _verbose(verbose)
    {
        static_assert(static_cast<int>(P::category) == static_cast<int>(L::category),
                      "P and L must have the same category!");
        static_assert(static_cast<int>(P::category) == static_cast<int>(S::category),
                      "P and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, double, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        const real_type EPSILON = 1e-80;
        const int m = _restart;
        real_type norm, norm_0;
        int j = 1;
        std::vector<field_type> s(m+1), sn(m);
        std::vector<real_type> cs(m);
        // need copy of rhs if GMRes has to be restarted
        Y b2(b);
        // helper vector
        Y w(b);
        std::vector< std::vector<field_type> > H(m+1, s);
        std::vector<F> v(m+1, b);

        // start timer
        Dune::Timer watch;
        watch.reset();

        // clear solver statistics and set res.converged to false
        res.clear();
        _W.pre(x, b);

        // calculate defect and overwrite rhs with it
        _A.applyscaleadd(-1.0, x, b); // b -= Ax
        // calculate preconditioned defect
        v[0] = 0.0; _W.apply(v[0], b); // r = W^-1 b
        norm_0 = _sp.norm(v[0]);
        norm = norm_0;

        this->convergenceCriterion().setInitial(x, b, norm_0);
        if (_verbose > 0) {
            std::cout << "=== RestartedGMResSolver" << std::endl << std::flush;
            if (_verbose > 1)
                this->convergenceCriterion().printInitial();
        }
        if (this->convergenceCriterion().converged())
            res.converged = true;

        while (j <= _maxit && res.converged != true) {

            int i = 0;
            v[0] *= 1.0/norm;
            s[0] = norm;
            for (i=1; i < m+1; i++)
                s[i] = 0.0;

            for (i=0; i < m && j <= _maxit && res.converged != true; i++, j++) {
                w = 0.0;
                // use v[i+1] as temporary vector
                v[i+1] = 0.0;
                // do Arnoldi algorithm
                _A.apply(v[i], v[i+1]);
                _W.apply(w, v[i+1]);
                for (int k=0; k < i+1; k++) {
                    // notice that _sp.dot (v[k], w) = v[k]\adjoint w
                    // so one has to pay attention to the order
                    // the in scalar product for the complex case
                    // doing the modified Gram-Schmidt algorithm
                    H[k][i] = _sp.dot(v[k], w);
                    // w -= H[k][i] * v[k]
                    w.axpy(-H[k][i], v[k]);
                }
                H[i+1][i] = _sp.norm(w);
                if (std::abs(H[i+1][i]) < EPSILON)
                    DUNE_THROW(Dune::ISTLError,
                               "breakdown in GMRes - |w| == 0.0 after " << j << " iterations");

                // normalize new vector
                v[i+1] = w; v[i+1] *= 1.0/H[i+1][i];

                // update QR factorization
                for (int k=0; k < i; k++)
                    applyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);

                // compute new givens rotation
                generatePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
                // finish updating QR factorization
                applyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
                applyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

                // norm of the defect is the last component the vector s
                norm = std::abs(s[i+1]);

                this->convergenceCriterion().update(x, b);
                if (_verbose > 1)
                    this->convergenceCriterion().print(j);
                if (this->convergenceCriterion().converged()) {
                    res.converged = true;
                }
            } // end for

            // calculate update vector
            w = 0.0;
            update(w, i, H, s, v);
            // and current iterate
            x += w;

            // restart GMRes if convergence was not achieved,
            // i.e. linear defect has not reached desired reduction
            // and if j < _maxit
            if (res.converged != true && j <= _maxit ) {

                if (_verbose > 0)
                    std::cout << "=== GMRes::restart" << std::endl;
                // get saved rhs
                b = b2;
                // calculate new defect
                _A.applyscaleadd(-1.0, x, b); // b -= Ax;
                // calculate preconditioned defect
                v[0] = 0.0;
                _W.apply(v[0], b);
                norm = _sp.norm(v[0]);
            }

        } //end while

        // postprocess preconditioner
        _W.post(x);

        // save solver statistics
        res.iterations = j-1; // it has to be j-1!!!
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();

        if (_verbose > 0)
            this->convergenceCriterion().print(j);
    }

private :
    void update(X& w, int i,
                const std::vector<std::vector<field_type> >& H,
                const std::vector<field_type>& s,
                const std::vector<X>& v) {
        // solution vector of the upper triangular system
        std::vector<field_type> y(s);

        // backsolve
        for (int a=i-1; a >=0; a--) {
            field_type rhs(s[a]);
            for (int b=a+1; b < i; b++)
                rhs -= H[a][b]*y[b];
            y[a] = rhs/H[a][a];

            // compute update on the fly
            // w += y[a]*v[a]
            w.axpy(y[a], v[a]);
        }
    }

    template <typename T>
    typename std::enable_if<std::is_same<field_type, real_type>::value, T>::type conjugate(const T& t) {
        return t;
    }

    template <typename T>
    typename std::enable_if<!std::is_same<field_type, real_type>::value, T>::type conjugate(const T& t) {
        return conj(t);
    }

    void
    generatePlaneRotation(field_type& dx, field_type& dy, real_type& cs, field_type& sn)
    {
        real_type norm_dx = std::abs(dx);
        real_type norm_dy = std::abs(dy);
        if (norm_dy < 1e-15) {
            cs = 1.0;
            sn = 0.0;
        } else if (norm_dx < 1e-15) {
            cs = 0.0;
            sn = 1.0;
        } else if (norm_dy > norm_dx) {
            real_type temp = norm_dx/norm_dy;
            cs = 1.0/std::sqrt(1.0 + temp*temp);
            sn = cs;
            cs *= temp;
            sn *= dx/norm_dx;
            sn *= conjugate(dy)/norm_dy;
        } else {
            real_type temp = norm_dy/norm_dx;
            cs = 1.0/std::sqrt(1.0 + temp*temp);
            sn = cs;
            sn *= conjugate(dy/dx);
        }
    }


    void
    applyPlaneRotation(field_type& dx, field_type& dy, real_type& cs, field_type& sn)
    {
        field_type temp = cs * dx + sn * dy;
        dy = -conjugate(sn) * dx + cs * dy;
        dx = temp;
    }

    Dune::LinearOperator<X, X> &_A;
    Dune::Preconditioner<X, X> &_W;
    Dune::SeqScalarProduct<X> ssp;
    Dune::ScalarProduct<X> &_sp;
    int _restart;
    int _maxit;
    int _verbose;
};


/**
 * @brief Generalized preconditioned conjugate gradient solver.
 *
 * A preconditioned conjugate gradient that allows
 * the preconditioner to change between iterations.
 *
 * One example for such preconditioner is AMG when used without
 * a direct coarse solver. In this case the number of iterations
 * performed on the coarsest level might change between applications.
 *
 * In contrast to CGSolver the search directions are stored and
 * the orthogonalization is done explicitly.
 */
template <class X>
class GeneralizedPCGSolver : public InverseOperator<X, X>
{
    typedef Ewoms::ConvergenceCriterion<X> ConvergenceCriterion;

public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
    //! \brief The real type of the field type (is the same if using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
#else
    typedef field_type real_type;
#endif

    /*!
      \brief Set up nonlinear preconditioned conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&, P&, double, int, int)
      \param restart When to restart the construction of
      the Krylov search space.
    */
    template <class L, class P>
    GeneralizedPCGSolver(L& op, P& prec, real_type reduction, int maxit, int verbose,
                         int restart = 10) :
        ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit),
        _verbose(verbose), _restart(std::min(maxit, restart))
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P have to have the same category!");
        static_assert(static_cast<int>(L::category) ==
                      static_cast<int>(Dune::SolverCategory::sequential),
                      "L has to be sequential!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }

    /*!
      \brief Set up nonlinear preconditioned conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&, S&, P&, double, int, int)
      \param restart When to restart the construction of
      the Krylov search space.
    */
    template <class L, class P, class S>
    GeneralizedPCGSolver(L& op, S& sp, P& prec,
                         real_type reduction, int maxit, int verbose, int restart=10) :
        _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose),
        _restart(std::min(maxit, restart))
    {
        static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
                      "L and P must have the same category!");
        static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
                      "L and S must have the same category!");

        auto crit = std::make_shared<ResidReductionCriterion<X>>(_sp, reduction);
        this->setConvergenceCriterion(crit);
    }
    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&, Y&, InverseOperatorResult&)
    */
    virtual void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
        res.clear(); // clear solver statistics
        Dune::Timer watch; // start a timer
        _prec.pre(x, b); // prepare preconditioner
        _op.applyscaleadd(-1, x, b); // overwrite b with defect

        std::vector<std::shared_ptr<X> > p(_restart);
        std::vector<typename X::field_type> pp(_restart);
        X q(x); // a temporary vector
        X prec_res(x); // a temporary vector for preconditioner output

        p[0].reset(new X(x));

        this->convergenceCriterion().setInitial(x, b);
        if (_verbose > 0) {
            std::cout << "=== GeneralizedPCGSolver" << std::endl << std::flush;
            if (_verbose > 1)
                this->convergenceCriterion().printInitial();
        }
        if (this->convergenceCriterion().converged()) {
            res.converged = true;
            res.iterations = 0;
            res.reduction = this->convergenceCriterion().accuracy();
            res.conv_rate = 0;
            res.elapsed = 0;
            return;
        }

        // some local variables
        field_type rho, lambda;

        int i = 0;
        int ii = 0;
        // determine initial search direction
        *(p[0]) = 0; // clear correction
        _prec.apply(*(p[0]), b); // apply preconditioner
        rho = _sp.dot(*(p[0]), b); // orthogonalization
        _op.apply(*(p[0]), q); // q=Ap
        pp[0] = _sp.dot(*(p[0]), q); // scalar product
        lambda = rho/pp[0]; // minimization
        x.axpy(lambda, *(p[0])); // update solution
        b.axpy(-lambda, q); // update defect

        // convergence test
        this->convergenceCriterion().update(x, b);
        if (_verbose > 1) // print
            this->convergenceCriterion().print(i);
        if (this->convergenceCriterion().converged()) {
            // fill statistics
            res.converged = true;
        }

        while (i < _maxit) {
            // the loop
            int end = std::min(_restart, _maxit-i+1);
            for (ii=1; ii < end;++ii )
            {
                //std::cout<<" ii="<<ii<<" i="<<i<<std::endl;
                // compute next conjugate direction
                prec_res = 0; // clear correction
                _prec.apply(prec_res, b); // apply preconditioner

                p[ii].reset(new X(prec_res));
                _op.apply(prec_res, q);

                for (int j=0; j < ii;++j) {
                    rho = _sp.dot(q, *(p[j]))/pp[j];
                    p[ii]->axpy(-rho, *(p[j]));
                }

                // minimize in given search direction
                _op.apply(*(p[ii]), q); // q=Ap
                pp[ii] = _sp.dot(*(p[ii]), q); // scalar product
                rho = _sp.dot(*(p[ii]), b); // orthogonalization
                lambda = rho/pp[ii]; // minimization
                x.axpy(lambda, *(p[ii])); // update solution
                b.axpy(-lambda, q); // update defect

                // convergence test
                this->convergenceCriterion().update(x, b);
                if (_verbose > 1) // print
                    this->convergenceCriterion().print(i);
                if (this->convergenceCriterion().converged()) {
                    res.converged = true;
                    break;
                }
            }
            if (res.converged)
                break;
            if (end == _restart) {
                *(p[0])=*(p[_restart-1]);
                pp[0]=pp[_restart-1];
            }
        }

        // postprocess preconditioner
        _prec.post(x);

        // fill statistics
        res.iterations = i;
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate = std::pow(res.reduction, 1.0/res.iterations);
        res.elapsed = watch.elapsed();

        if (_verbose == 1) // printing in non-verbose mode
            this->convergenceCriterion().print(i);
    }

private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X, X> &_op;
    Dune::Preconditioner<X, X> &_prec;
    Dune::ScalarProduct<X> &_sp;
    int _maxit;
    int _verbose;
    int _restart;
};

/** \} end documentation */

} // end namespace

#endif
