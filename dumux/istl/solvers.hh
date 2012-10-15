// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
/*****************************************************************************
 *   Copyright (C) 2004-2012 by Christian Engwer                             *
 *   Copyright (C) 2005-2011 by Markus Blatt                                 *
 *   Copyright (C) 2007 by Robert Kloefkorn                                  *
 *   Copyright (C) 2004-2006 by Peter Bastian                                *
 *   Copyright (C) 2010 by Jorrit Fahlke                                     *
 *   Copyright (C) 2007-2012 by Oliver Sander                                *
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, version 2.                                *
 *                                                                           *
 *   As a special exception, you may use the DUNE library without            *
 *   restriction.  Specifically, if other files instantiate templates or     *
 *   use macros or inline functions from one or more of the DUNE source      *
 *   files, or you compile one or more of the DUNE source files and link     *
 *   them with other files to produce an executable, this does not by        *
 *   itself cause the resulting executable to be covered by the GNU          *
 *   General Public License.  This exception does not however invalidate     *
 *   any other reasons why the executable file might be covered by the       *
 *   GNU General Public License.                                             *
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
 * \brief Copy of dune-istl's linear solvers with support for
 *        pluggable convergence criteria added.
 *
 * For eWoms, pluggable convergence criteria for the linear solvers
 * are an important feature. Unfortunatly, the DUNE developers don't
 * seem to care, so this could not go directly into ISTL. For the sad
 * details, see
 *
 * http://www.dune-project.org/flyspray/index.php?do=details&task_id=1018
 *
 * Because we don't want to let ourselfs bugged down by dune (pardon the pun),
 * we just provide a patched copy of their linear solvers.
 */
#ifndef DUMUX_SOLVERS_HH
#define DUMUX_SOLVERS_HH

#include<cmath>
#include<complex>
#include<iostream>
#include<iomanip>
#include<string>

#include "residreductioncriterion.hh"
#include "weightedresidreductioncriterion.hh"
#include "fixpointcriterion.hh"

#include <dune/istl/istlexception.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>

#include <dune/common/timer.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/static_assert.hh>

#include <memory>

namespace Dumux {
  /** @defgroup ISTL_Solvers Iterative Linear Solvers
      @ingroup ISTL
  */
  /** @addtogroup ISTL_Solvers
      \{
  */

  /*!
      \brief Statistics about the application of an inverse operator

      The return value of an application of the inverse
      operator delivers some important information about
      the iteration.
  */
  struct InverseOperatorResult
  {
    /** \brief Default constructor */
    InverseOperatorResult ()
    {
      clear();
    }

    /** \brief Resets all data */
    void clear ()
    {
      iterations = 0;
      reduction = 0;
      converged = false;
      conv_rate = 1;
      elapsed = 0;
    }

    /** \brief Number of iterations */
    int iterations;

    /** \brief Reduction achieved: \f$ \|b-A(x^n)\|/\|b-A(x^0)\|\f$ */
    double reduction;

    /** \brief True if convergence criterion has been met */
    bool converged;

    /** \brief Convergence rate (average reduction per step) */
    double conv_rate;

    /** \brief Elapsed time in seconds */
    double elapsed;
  };


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
  template<class X, class Y>
  class InverseOperator {
  public:
    //! \brief Type of the domain of the operator to be inverted.
    typedef X domain_type;

    //! \brief Type of the range of the operator to be inverted.
    typedef Y range_type;

    /** \brief The field type of the operator. */
    typedef typename X::field_type field_type;

    /*!
     * \brief Return the criterion to be used to check for convergence of the linear solver.
     */
    virtual Dumux::ConvergenceCriterion<X> &convergenceCriterion()
    { return *convergenceCriterion_; }

    /*!
     * \copydoc convergenceCriterion()
     */
    virtual const Dumux::ConvergenceCriterion<X> &convergenceCriterion() const
    { return *convergenceCriterion_; }

    /*!
     * \brief Set the criterion to be used to check for convergence of the linear solver.
     */
    virtual void setConvergenceCriterion(Dune::shared_ptr<Dumux::ConvergenceCriterion<X> > convCrit)
    { convergenceCriterion_ = convCrit; }

    /*!
        \brief Apply inverse operator,

        \warning Note: right hand side b may be overwritten!

        \param x The left hand side to store the result in.
        \param b The right hand side
        \param res Object to store the statistics about applying the operator.
    */
    virtual void apply (X& x, Y& b, InverseOperatorResult& res) = 0;

    /*!
      \brief apply inverse operator, with given convergence criteria.

      \warning Right hand side b may be overwritten!

      \param x The left hand side to store the result in.
      \param b The right hand side
      \param reduction The minimum defect reduction to achieve.
      \param res Object to store the statistics about applying the operator.
    */
    virtual void apply (X& x, Y& b, double reduction, InverseOperatorResult& res) = 0;

    //! \brief Destructor
    virtual ~InverseOperator () {}

  private:
    Dune::shared_ptr<ConvergenceCriterion<X> > convergenceCriterion_;
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
  template<class X>
  class LoopSolver : public InverseOperator<X,X> {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
  public:
    //! \brief The domain type of the operator that we do the inverse for.
    typedef X domain_type;
    //! \brief The range type of the operator that we do the inverse for.
    typedef X range_type;
    //! \brief The field type of the operator that we do the inverse for.
    typedef typename X::field_type field_type;

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
    template<class L, class P>
    LoopSolver (L& op, P& prec,
                double reduction, int maxit, int verbose) :
      ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P have to have the same category!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
        "L has to be sequential!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }

    /*!
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
    template<class L, class S, class P>
    LoopSolver (L& op, S& sp, P& prec,
      double reduction, int maxit, int verbose) :
      _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(S::category),
        "L and S must have the same category!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }

    //! \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      // clear solver statistics
      res.clear();

      // start a timer
      Dune::Timer watch;

      // prepare preconditioner
      _prec.pre(x,b);

      // overwrite b with defect
      _op.applyscaleadd(-1,x,b);

      // compute norm, \todo parallelization
      this->convergenceCriterion().setInitial(x, b);

      // printing
      if (_verbose>0)
      {
        std::cout << "=== LoopSolver" << std::endl;
        if (_verbose>1)
        {
          this->convergenceCriterion().printInitial();
        }
      }

      // allocate correction vector
      X v(x);

      // iteration loop
      int i=1;
      for ( ; i<=_maxit; i++ )
      {
        v = 0;                      // clear correction
        _prec.apply(v,b);           // apply preconditioner
        x += v;                     // update solution
        _op.applyscaleadd(-1,v,b);  // update defect

        field_type lastAccuracy = this->convergenceCriterion().accuracy();
        this->convergenceCriterion().update(x, b);

        if (_verbose>1)             // print
          this->convergenceCriterion().print(i);

        if (this->convergenceCriterion().converged())
        {
          res.converged  = true;
          break;
        }
      }

      // print
      if (_verbose==1)
        this->convergenceCriterion().print(i);

      // postprocess preconditioner
      _prec.post(x);

      // fill statistics
      res.iterations = i;
      res.reduction = this->convergenceCriterion().accuracy();
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();

      // final print
      if (_verbose>0)
      {
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/i
                  << ", IT=" << i << std::endl;
      }
    }

    //! \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res)
    {
      double origTol = this->convergenceCriterion().tolerance();
      this->convergenceCriterion().setTolerance(reduction);
      (*this).apply(x,b,res);
      this->convergenceCriterion().setTolerance(origTol);
    }

  private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X,X>& _op;
    Dune::Preconditioner<X,X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    Dune::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
  };


  // all these solvers are taken from the SUMO library
  //! gradient method
  template<class X >
  class GradientSolver : public InverseOperator<X,X> {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
  public:
    //! \brief The domain type of the operator that we do the inverse for.
    typedef X domain_type;
    //! \brief The range type of the operator  that we do the inverse for.
    typedef X range_type;
    //! \brief The field type of the operator  that we do the inverse for.
    typedef typename X::field_type field_type;

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
    */
    template<class L, class P>
    GradientSolver (L& op, P& prec,
      double reduction, int maxit, int verbose) :
      ssp(), _op(op), _prec(prec), _sp(ssp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P have to have the same category!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
        "L has to be sequential!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }
    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
    */
    template<class L, class S, class P>
    GradientSolver (L& op, S& sp, P& prec,
      double reduction, int maxit, int verbose) :
      _op(op), _prec(prec), _sp(sp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P have to have the same category!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(S::category),
        "L and S have to have the same category!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      res.clear();                  // clear solver statistics
      Dune::Timer watch;                // start a timer
      _prec.pre(x,b);             // prepare preconditioner
      _op.applyscaleadd(-1,x,b);  // overwrite b with defect

      X p(x);                     // create local vectors
      X q(b);

      this->convergenceCriterion().setInitial(x, b);

      if (_verbose>0)             // printing
      {
        std::cout << "=== GradientSolver" << std::endl;
        if (_verbose>1)
        {
          this->convergenceCriterion().printInitial();
        }
      }

      int i=1;  // loop variables
      field_type lambda;
      for ( ; i<=_maxit; i++ )
      {
        p = 0;                      // clear correction
        _prec.apply(p,b);           // apply preconditioner
        _op.apply(p,q);             // q=Ap
        lambda = _sp.dot(p,b)/_sp.dot(q,p);// minimization
        x.axpy(lambda,p);           // update solution
        b.axpy(-lambda,q);          // update defect

        field_type lastAccuracy = this->convergenceCriterion().accuracy();
        this->convergenceCriterion().update(x, b);

        if (_verbose>1)             // print
          this->convergenceCriterion().print(i);

        if (this->convergenceCriterion().converged())
        {
          res.converged  = true;
          break;
        }
      }

      if (_verbose==1)                // printing for non verbose
        this->convergenceCriterion().print(i);

      _prec.post(x);                  // postprocess preconditioner
      res.iterations = i;               // fill statistics
      res.reduction = this->convergenceCriterion().accuracy();
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();
      if (_verbose>0)                 // final print
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/i
                  << ", IT=" << i << std::endl;
    }

    /*!
      \brief Apply inverse operator with given reduction factor.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res)
    {
      std::swap(_reduction,reduction);
      (*this).apply(x,b,res);
      std::swap(_reduction,reduction);
    }

  private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X,X>& _op;
    Dune::Preconditioner<X,X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    Dune::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    double _reduction;
    int _maxit;
    int _verbose;
  };



  //! \brief conjugate gradient method
  template<class X>
  class CGSolver : public InverseOperator<X,X> {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
  public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;

    /*!
      \brief Set up conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
    */
    template<class L, class P>
    CGSolver (L& op, P& prec, double reduction, int maxit, int verbose) :
      ssp(), _op(op), _prec(prec), _sp(ssp), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
        "L must be sequential!");
      
      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }
    /*!
      \brief Set up conjugate gradient solver.

      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
    */
    template<class L, class S, class P>
    CGSolver (L& op, S& sp, P& prec, double reduction, int maxit, int verbose) :
      _op(op), _prec(prec), _sp(sp), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
                          "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(S::category),
                          "L and S must have the same category!");
      
      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
     }
 
    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      res.clear();                  // clear solver statistics
      Dune::Timer watch;                // start a timer
      _prec.pre(x,b);             // prepare preconditioner
      _op.applyscaleadd(-1,x,b);  // overwrite b with defect

      X p(x);              // the search direction
      X q(x);              // a temporary vector

      this->convergenceCriterion().setInitial(x, b);

      if (_verbose>0)             // printing
      {
        std::cout << "=== CGSolver" << std::endl;
        if (_verbose>1) {
          this->convergenceCriterion().printInitial();
        }
      }

      // some local variables
      field_type rho,rholast,lambda,alpha,beta;

      // determine initial search direction
      p = 0;                          // clear correction
      _prec.apply(p,b);               // apply preconditioner
      rholast = _sp.dot(p,b);         // orthogonalization

      // the loop
      int i=1;
      for ( ; i<=_maxit; i++ )
      {
        // minimize in given search direction p
        _op.apply(p,q);             // q=Ap
        alpha = _sp.dot(p,q);       // scalar product
        lambda = rholast/alpha;     // minimization
        x.axpy(lambda,p);           // update solution
        b.axpy(-lambda,q);          // update defect

        // convergence test
        field_type lastAccuracy = this->convergenceCriterion().accuracy();

        this->convergenceCriterion().update(x, b);

        if (_verbose>1)             // print
          this->convergenceCriterion().print(i);

        if (this->convergenceCriterion().converged())
        {
          res.converged  = true;
          break;
        }

        // determine new search direction
        q = 0;                      // clear correction
        _prec.apply(q,b);           // apply preconditioner
        rho = _sp.dot(q,b);         // orthogonalization
        beta = rho/rholast;         // scaling factor
        p *= beta;                  // scale old search direction
        p += q;                     // orthogonalization with correction
        rholast = rho;              // remember rho for recurrence
      }

      if (_verbose==1)                // printing for non verbose
        this->convergenceCriterion().print(i);

      _prec.post(x);                  // postprocess preconditioner
      res.iterations = i;               // fill statistics
      res.reduction = this->convergenceCriterion().accuracy();
      res.conv_rate  = pow(res.reduction,1.0/i);
      res.elapsed = watch.elapsed();

      if (_verbose>0)                 // final print
      {
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/i
                  << ", IT=" << i << std::endl;
      }
    }

    /*!
      \brief Apply inverse operator with given reduction factor.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, double reduction,
      InverseOperatorResult& res)
    {
      double origTol = this->convergenceCriterion().tolerance();
      this->convergenceCriterion().setTolerance(reduction);
      (*this).apply(x,b,res);
      this->convergenceCriterion().setTolerance(origTol);
    }

  private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X,X>& _op;
    Dune::Preconditioner<X,X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    Dune::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
  };

// Ronald Kriemanns BiCG-STAB implementation from Sumo
//! \brief Bi-conjugate Gradient Stabilized (BiCG-STAB)
template<class X>
class BiCGSTABSolver : public InverseOperator<X,X> {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted
    typedef typename X::field_type field_type;

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
    */
    template<class L, class P>
    BiCGSTABSolver(L& op, 
                   P& prec,
                   double reduction,
                   int maxit, 
                   int verbose) 
      : ssp()
      , _op(op)
      , _prec(prec)
      , _sp(ssp)
      , _maxit(maxit)
      , _verbose(verbose)
    {
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(P::category), "L and P must be of the same category!");
      dune_static_assert(static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential), "L must be sequential!");
      
      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }
    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
    */
    template<class L, class S, class P>
    BiCGSTABSolver(L& op, 
                   S& sp,
                   P& prec,
                   double reduction,
                   int maxit,
                   int verbose) 
      : _op(op)
      , _prec(prec)
      , _sp(sp)
      , _maxit(maxit)
      , _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(S::category),
        "L and S must have the same category!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
     }
  
    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      const double EPSILON=1e-80;

      double               it;
      field_type           rho, rho_new, alpha, beta, h, omega;

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
      res.clear();                // clear solver statistics
      Dune::Timer watch;                // start a timer
      _prec.pre(x,r);             // prepare preconditioner
      _op.applyscaleadd(-1,x,r);  // overwrite b with defect

      rt=r;

      p=0;
      v=0;

      rho   = 1;
      alpha = 1;
      omega = 1;

      this->convergenceCriterion().setInitial(x, r);

      if (_verbose>0)             // printing
      {
        std::cout << "=== BiCGSTABSolver" << std::endl;
        if (_verbose>1)
          this->convergenceCriterion().printInitial();
      }

      //
      // iteration
      //

      for (it = 0.5; it < _maxit; it+=.5)
      {
        //
        // preprocess, set vecsizes etc.
        //

        // rho_new = < rt , r >
        rho_new = _sp.dot(rt,r);

        // look if breakdown occured
        if (std::abs(rho) <= EPSILON)
          DUNE_THROW(Dune::ISTLError,"breakdown in BiCGSTAB - rho "
            << rho << " <= EPSILON " << EPSILON
            << " after " << it << " iterations");
        if (std::abs(omega) <= EPSILON)
          DUNE_THROW(Dune::ISTLError,"breakdown in BiCGSTAB - omega "
            << omega << " <= EPSILON " << EPSILON
            << " after " << it << " iterations");


        if (it<1)
          p = r;
        else
        {
          beta = ( rho_new / rho ) * ( alpha / omega );
          p.axpy(-omega,v); // p = r + beta (p - omega*v)
          p *= beta;
          p += r;
        }

        // y = W^-1 * p
        y = 0;
        _prec.apply(y,p);           // apply preconditioner

        // v = A * y
        _op.apply(y,v);

        // alpha = rho_new / < rt, v >
        h = _sp.dot(rt,v);

        if ( std::abs(h) < EPSILON )
          DUNE_THROW(Dune::ISTLError,"h=0 in BiCGSTAB");

        alpha = rho_new / h;

        // apply first correction to x
        // x <- x + alpha y
        x.axpy(alpha,y);

        // r = r - alpha*v
        r.axpy(-alpha,v);

        //
        // test stop criteria
        //

        field_type lastAccuracy = this->convergenceCriterion().accuracy();
        this->convergenceCriterion().update(x, r);

        if (_verbose>1) // print
          this->convergenceCriterion().print(it);

        if (this->convergenceCriterion().converged())
        {
          res.converged = 1;
          break;

        }
        it+=.5;

        // y = W^-1 * r
        y = 0;
        _prec.apply(y,r);

        // t = A * y
        _op.apply(y,t);

        // omega = < t, r > / < t, t >
        omega = _sp.dot(t,r)/_sp.dot(t,t);

        // apply second correction to x
        // x <- x + omega y
        x.axpy(omega,y);

        // r = s - omega*t (remember : r = s)
        r.axpy(-omega,t);

        rho = rho_new;

        //
        // test stop criteria
        //

        lastAccuracy = this->convergenceCriterion().accuracy();
        this->convergenceCriterion().update(x, r);

        if (_verbose > 1)             // print
          this->convergenceCriterion().print(it);

        if (this->convergenceCriterion().converged())
        {
          res.converged = 1;
          break;
        }
      } // end for

      if (_verbose==1)                // printing for non verbose
        this->convergenceCriterion().print(it);

      _prec.post(x);                  // postprocess preconditioner
      res.iterations = static_cast<int>(std::ceil(it));              // fill statistics
      res.reduction = this->convergenceCriterion().accuracy();
      res.conv_rate  = pow(res.reduction,1.0/it);
      res.elapsed = watch.elapsed();
      if (_verbose>0)                 // final print
        std::cout << "=== rate=" << res.conv_rate
                  << ", T=" << res.elapsed
                  << ", TIT=" << res.elapsed/it
                  << ", IT=" << it << std::endl;
    }

    /*!
      \brief Apply inverse operator with given reduction factor.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res)
    {
      double origTol = this->convergenceCriterion().tolerance();
      this->convergenceCriterion().setTolerance(reduction);
      (*this).apply(x,b,res);
      this->convergenceCriterion().setTolerance(origTol);
    }

  private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X,X>& _op;
    Dune::Preconditioner<X,X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    Dune::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    int _maxit;
    int _verbose;
  };

  /*! \brief Minimal Residual Method (MINRES)

    Symmetrically Preconditioned MINRES as in A. Greenbaum, 'Iterative Methods for Solving Linear Systems', pp. 121
    Iterative solver for symmetric indefinite operators.
    Note that in order to ensure the (symmetrically) preconditioned system to remain symmetric, the preconditioner has to be spd.
  */
  template<class X>
  class MINRESSolver : public InverseOperator<X,X> {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
  public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef X range_type;
    //! \brief The field type of the operator to be inverted.
    typedef typename X::field_type field_type;
    //! \brief The real type of the field type (is the same of using real numbers, but differs for std::complex)
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;

    /*!
      \brief Set up MINRES solver.

      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
    */
    template<class L, class P>
    MINRESSolver (L& op, P& prec, double reduction, int maxit, int verbose) :
      ssp(), _op(op), _prec(prec), _sp(ssp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
        "L must be sequential!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }
    /*!
      \brief Set up MINRES solver.

      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
    */
    template<class L, class S, class P>
    MINRESSolver (L& op, S& sp, P& prec, double reduction, int maxit, int verbose) :
      _op(op), _prec(prec), _sp(sp), _reduction(reduction), _maxit(maxit), _verbose(verbose)
    {
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(P::category),
        "L and P must have the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(S::category),
        "L and S must have the same category!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }
   
    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      res.clear();                // clear solver statistics
      Dune::Timer watch;                // start a timer
      _prec.pre(x,b);             // prepare preconditioner
      _op.applyscaleadd(-1,x,b);  // overwrite b with defect/residual

      this->convergenceCriterion().setInitial(x, b);

      if (_verbose>0)             // printing
      {
        std::cout << "=== MINRESSolver" << std::endl;
        if (_verbose>1) {
          this->convergenceCriterion().printInitial();
        }
      }

      // some local variables
      field_type alpha,                   // recurrence coefficients as computed in the Lanczos alg making up the matrix T
        c[2]={0.0, 0.0},         // diagonal entry of Givens rotation
        s[2]={0.0, 0.0};         // off-diagonal entries of Givens rotation
        real_type beta;

        field_type T[3]={0.0, 0.0, 0.0};    // recurrence coefficients (column k of Matrix T)

        X   z(b.size()),    // some temporary vectors
          dummy(b.size());

        field_type xi[2]={1.0, 0.0};

        // initialize
        z = 0.0;                // clear correction

        _prec.apply(z,b);       // apply preconditioner z=M^-1*b

        beta = std::sqrt(std::abs(_sp.dot(z,b)));
        real_type beta0 = beta;

        X p[3];     // the search directions
        X q[3];     // Orthonormal basis vectors (in unpreconditioned case)

        q[0].resize(b.size());
        q[1].resize(b.size());
        q[2].resize(b.size());
        q[0] = 0.0;
        q[1] = b;
        q[1] /= beta;
        q[2] = 0.0;

        p[0].resize(b.size());
        p[1].resize(b.size());
        p[2].resize(b.size());
        p[0] = 0.0;
        p[1] = 0.0;
        p[2] = 0.0;


        z /= beta;      // this is w_current

        // the loop
        int i=1;
        for ( ; i<=_maxit; i++)
        {
          dummy = z; // remember z_old for the computation of the search direction p in the next iteration

          int i1 = i%3,
            i0 = (i1+2)%3,
            i2 = (i1+1)%3;

          // Symmetrically Preconditioned Lanczos (Greenbaum p.121)
          _op.apply(z,q[i2]);             // q[i2] = Az
          q[i2].axpy(-beta, q[i0]);
          alpha = _sp.dot(q[i2],z);
          q[i2].axpy(-alpha, q[i1]);

          z=0.0;
          _prec.apply(z,q[i2]);

          beta = std::sqrt(std::abs(_sp.dot(q[i2],z)));

          q[i2] /= beta;
          z /= beta;

          // QR Factorization of recurrence coefficient matrix
          // apply previous Givens rotations to last column of T
          T[1] = T[2];
          if (i>2)
          {
            T[0] = s[i%2]*T[1];
            T[1] = c[i%2]*T[1];
          }
          if (i>1)
          {
            T[2] = c[(i+1)%2]*alpha - s[(i+1)%2]*T[1];
            T[1] = c[(i+1)%2]*T[1] + s[(i+1)%2]*alpha;
          }
          else
            T[2] = alpha;

          // recompute c, s -> current Givens rotation \TODO use BLAS-routine drotg instead for greater robustness
//          cblas_drotg (a, b, c, s);
          c[i%2] = 1.0/std::sqrt(T[2]*T[2] + beta*beta);
          s[i%2] = beta*c[i%2];
          c[i%2] *= T[2];

          // apply current Givens rotation to T eliminating the last entry...
          T[2] = c[i%2]*T[2] + s[i%2]*beta;

          // ...and to xi, the right hand side of the least squares problem min_y||beta*xi-T*y||
          xi[i%2] = -s[i%2]*xi[(i+1)%2];
          xi[(i+1)%2] *= c[i%2];

          // compute correction direction
          p[i2] = dummy;
          p[i2].axpy(-T[1],p[i1]);
          p[i2].axpy(-T[0],p[i0]);
          p[i2] /= T[2];

          // apply correction/update solution
          x.axpy(beta0*xi[(i+1)%2], p[i2]);

          // remember beta_old
          T[2] = beta;

          // update residual - not necessary if in the preconditioned case we are content with the residual norm of the
          // preconditioned system as convergence test
//          _op.apply(p[i2],dummy);
//          b.axpy(-beta0*xi[(i+1)%2],dummy);

//          convergence test
          field_type lastAccuracy = this->convergenceCriterion().accuracy();
          this->convergenceCriterion().update(x, b);

          if (_verbose>1)             // print
            this->convergenceCriterion().print(i);

          if (this->convergenceCriterion().converged())
          {
            res.converged  = true;
            break;
          }
        }

        if (_verbose==1)                // printing for non verbose
          this->convergenceCriterion().print(i);

        _prec.post(x);                  // postprocess preconditioner
        res.iterations = i;               // fill statistics
        res.reduction = this->convergenceCriterion().accuracy();
        res.conv_rate  = pow(res.reduction,1.0/i);
        res.elapsed = watch.elapsed();

        if (_verbose>0)                 // final print
        {
          std::cout << "=== rate=" << res.conv_rate
                    << ", T=" << res.elapsed
                    << ", TIT=" << res.elapsed/i
                    << ", IT=" << i << std::endl;
        }

    }

    /*!
      \brief Apply inverse operator with given reduction factor.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res)
    {
      double origTol = this->convergenceCriterion().tolerance();
      this->convergenceCriterion().setTolerance(reduction);
      (*this).apply(x,b,res);
      this->convergenceCriterion().setTolerance(origTol);
    }

  private:
    Dune::SeqScalarProduct<X> ssp;
    Dune::LinearOperator<X,X>& _op;
    Dune::Preconditioner<X,X>& _prec;
    Dune::ScalarProduct<X>& _sp;
    Dune::shared_ptr<ConvergenceCriterion> convergenceCriterion_;
    double _reduction;
    int _maxit;
    int _verbose;
  };

  /*!
     \brief implements the Generalized Minimal Residual (GMRes) method

     GMRes solves the unsymmetric linear system Ax = b using the
     Generalized Minimal Residual method as described the SIAM Templates
     book (http://www.netlib.org/templates/templates.pdf).

     \todo construct F via rebind and an appropriate field_type

  */

  template<class X, class Y=X, class F = Y>
  class RestartedGMResSolver : public InverseOperator<X,Y>
  {
    typedef Dumux::ConvergenceCriterion<X> ConvergenceCriterion;
  public:
    //! \brief The domain type of the operator to be inverted.
    typedef X domain_type;
    //! \brief The range type of the operator to be inverted.
    typedef Y range_type;
    //! \brief The field type of the operator to be inverted
    typedef typename X::field_type field_type;
    //! \brief The field type of the basis vectors
    typedef F basis_type;

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,P&,double,int,int)
      \param restart number of GMRes cycles before restart
      \param recalc_defect recalculate the defect after everey restart or not [default=false]
    */
    template<class L, class P>
    RestartedGMResSolver (L& op, P& prec, double reduction, int restart, int maxit, int verbose, bool recalc_defect = false) :
      _A_(op), _M(prec),
      ssp(), _sp(ssp), _restart(restart),
      _reduction(reduction), _maxit(maxit), _verbose(verbose),
      _recalc_defect(recalc_defect)
    {
      dune_static_assert(static_cast<int>(P::category) == static_cast<int>(L::category),
        "P and L must be the same category!");
      dune_static_assert( static_cast<int>(L::category) == static_cast<int>(Dune::SolverCategory::sequential),
        "L must be sequential!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }

    /*!
      \brief Set up solver.

      \copydoc LoopSolver::LoopSolver(L&,S&,P&,double,int,int)
      \param restart number of GMRes cycles before restart
      \param recalc_defect recalculate the defect after everey restart or not [default=false]
    */
    template<class L, class S, class P>
    RestartedGMResSolver (L& op, S& sp, P& prec, double reduction, int restart, int maxit, int verbose, bool recalc_defect = false) :
      _A_(op), _M(prec),
      _sp(sp), _restart(restart),
      _reduction(reduction), _maxit(maxit), _verbose(verbose),
      _recalc_defect(recalc_defect)
    {
      dune_static_assert(static_cast<int>(P::category) == static_cast<int>(L::category),
        "P and L must have the same category!");
      dune_static_assert(static_cast<int>(P::category) == static_cast<int>(S::category),
        "P and S must have the same category!");

      this->setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(new ResidReductionCriterion<X>(_sp, reduction)));
    }  

    //! \copydoc InverseOperator::apply(X&,Y&,InverseOperatorResult&)
    virtual void apply (X& x, X& b, InverseOperatorResult& res)
    {
      apply(x,b,_reduction,res);
    }

    /*!
      \brief Apply inverse operator.

      \copydoc InverseOperator::apply(X&,Y&,double,InverseOperatorResult&)
    */
    virtual void apply (X& x, Y& b, double reduction, InverseOperatorResult& res)
    {
      int m = _restart;
      field_type norm;
      field_type norm_old = 0.0;
      field_type norm_0;
      field_type beta;
      int i, j = 1, k;
      std::vector<field_type> s(m+1), cs(m), sn(m);
      // helper vector
      X w(b);
      std::vector< std::vector<field_type> > H(m+1,s);
      std::vector<F> v(m+1,b);

      // start timer
      Dune::Timer watch;                // start a timer

      // clear solver statistics
      res.clear();
      _M.pre(x,b);
      if (_recalc_defect)
      {
        // norm_0 = norm(M^-1 b)
        w = 0.0; _M.apply(w,b); // w = M^-1 b
        norm_0 = _sp.norm(w);
        // r = _M.solve(b - A * x);
        w = b;
        _A_.applyscaleadd(-1,x, /* => */ w); // w = b - Ax;
        v[0] = 0.0; _M.apply(v[0],w); // r = M^-1 w
        beta = _sp.norm(v[0]);
      }
      else
      {
        // norm_0 = norm(M^-1 b)
        w = 0.0; _M.apply(w,b); // w = M^-1 b
        // r = _M.solve(b - A * x);
        _A_.applyscaleadd(-1,x, /* => */ b); // b = b - Ax;
        norm_0 = _sp.norm(b);
        v[0] = 0.0; _M.apply(v[0],b); // r = M^-1 b
        beta = _sp.norm(v[0]);
      }

      // avoid division by zero
      if (norm_0 == 0.0)
        norm_0 = 1.0;
      norm = norm_old = _sp.norm(v[0]);

      // check convergence
      this->convergenceCriterion().setInitial(x, b, norm_0);

      // print header
      if (_verbose > 0)
      {
        std::cout << "=== RestartedGMResSolver" << std::endl;
        if (_verbose > 1)
        {
          this->convergenceCriterion().printInitial();
        }
      }

      while (j <= _maxit && res.converged != true) {
        v[0] *= (1.0 / beta);
        for (i=1; i<=m; i++) s[i] = 0.0;
        s[0] = beta;

        for (i = 0; i < m && j <= _maxit && res.converged != true; i++, j++) {
          w = 0.0;
          v[i+1] = 0.0; // use v[i+1] as temporary vector
          _A_.apply(v[i], /* => */ v[i+1]);
          _M.apply(w, v[i+1]);
          for (k = 0; k <= i; k++) {
            H[k][i] = _sp.dot(w, v[k]);
            // w -= H[k][i] * v[k];
            w.axpy(-H[k][i], v[k]);
          }
          H[i+1][i] = _sp.norm(w);
          if (H[i+1][i] == 0.0)
            DUNE_THROW(Dune::ISTLError,"breakdown in GMRes - |w| "
              << w << " == 0.0 after " << j << " iterations");
          // v[i+1] = w * (1.0 / H[i+1][i]);
          v[i+1] = w; v[i+1] *= (1.0 / H[i+1][i]);

          for (k = 0; k < i; k++)
            applyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);

          generatePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
          applyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
          applyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

          norm = std::abs(s[i+1]);
          norm_old = norm;

          field_type lastAccuracy = this->convergenceCriterion().accuracy();
          this->convergenceCriterion().update(x, b, norm);
          if (_verbose > 1)             // print
          {
            this->convergenceCriterion().print(i);
          }

          if (this->convergenceCriterion().converged()) {
            res.converged = true;
          }
        }

        if (_recalc_defect)
        {
          // update x
          update(x, i - 1, H, s, v);

          // update residuum
          // r = M^-1 (b - A * x);
          w = b; _A_.applyscaleadd(-1,x, /* => */ w);
          _M.apply(v[0], w);
          beta = _sp.norm(v[0]);
          norm = beta;
        }
        else
        {
          // calc update vector
          w = 0;
          update(w, i - 1, H, s, v);

          // update x
          x += w;

          // r = M^-1 (b - A * x);
          // update defect
          _A_.applyscaleadd(-1,w, /* => */ b);
          // r = M^-1 (b - A * x);
          v[0] = 0.0; _M.apply(v[0],b); // r = M^-1 b
          beta = _sp.norm(v[0]);
          norm = beta;

          res.converged = false;
        }

        norm_old = norm;

        field_type lastAccuracy = this->convergenceCriterion().accuracy();
        this->convergenceCriterion().update(x, b, norm);

        if (_verbose > 1)             // print
        {
          this->convergenceCriterion().print(i);
        }

        if (this->convergenceCriterion().converged()) {
          // fill statistics
          res.converged = true;
        }

        if (!res.converged && _verbose > 0)
          std::cout << "=== GMRes::restart\n";
      }

      _M.post(x);                  // postprocess preconditioner

      res.iterations = j;
      res.reduction = this->convergenceCriterion().accuracy();
      res.conv_rate  = pow(res.reduction,1.0/j);
      res.elapsed = watch.elapsed();

      if (_verbose>0)
        print_result(res);
    }
  private:

    void
    print_result (const InverseOperatorResult & res) const
    {
      int j = res.iterations>0?res.iterations:1;
      std::cout << "=== rate=" << res.conv_rate
                << ", T=" << res.elapsed
                << ", TIT=" << res.elapsed/j
                << ", IT=" << res.iterations
                << std::endl;
    }

    static void
    update(X &x, int k,
      std::vector< std::vector<field_type> > & h,
      std::vector<field_type> & s, std::vector<F> v)
    {
      std::vector<field_type> y(s);

      // Backsolve:
      for (int i = k; i >= 0; i--) {
        y[i] /= h[i][i];
        for (int j = i - 1; j >= 0; j--)
          y[j] -= h[j][i] * y[i];
      }

      for (int j = 0; j <= k; j++)
        // x += v[j] * y[j];
        x.axpy(y[j],v[j]);
    }

    void
    generatePlaneRotation(field_type &dx, field_type &dy, field_type &cs, field_type &sn)
    {
      if (dy == 0.0) {
        cs = 1.0;
        sn = 0.0;
      } else if (std::abs(dy) > std::abs(dx)) {
        field_type temp = dx / dy;
        sn = 1.0 / std::sqrt( 1.0 + temp*temp );
        cs = temp * sn;
      } else {
        field_type temp = dy / dx;
        cs = 1.0 / std::sqrt( 1.0 + temp*temp );
        sn = temp * cs;
      }
    }


    void
    applyPlaneRotation(field_type &dx, field_type &dy, field_type &cs, field_type &sn)
    {
      field_type temp  =  cs * dx + sn * dy;
      dy = -sn * dx + cs * dy;
      dx = temp;
    }

    Dune::LinearOperator<X,X>& _A_;
    Dune::Preconditioner<X,X>& _M;
    Dune::SeqScalarProduct<X> ssp;
    Dune::ScalarProduct<X>& _sp;
    int _restart;
    double _reduction;
    int _maxit;
    int _verbose;
    bool _recalc_defect;
  };

  /** \} end documentation */

} // end namespace

#endif
