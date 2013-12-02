// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 *
 * \copydoc Ewoms::VcfvNewtonMethod
 */
#ifndef EWOMS_VCFV_NEWTON_METHOD_HH
#define EWOMS_VCFV_NEWTON_METHOD_HH

#include "vcfvnewtonconvergencewriter.hh"

#include <ewoms/nonlinear/newtonmethod.hh>
#include <opm/core/utility/PropertySystem.hpp>

namespace Ewoms {

template <class TypeTag>
class VcfvNewtonMethod;

template <class TypeTag>
class VcfvNewtonConvergenceWriter;
} // namespace Ewoms

namespace Opm {
namespace Properties {
//! create a type tag for the VCFV specific Newton method
NEW_TYPE_TAG(VcfvNewtonMethod, INHERITS_FROM(NewtonMethod));

//! The class dealing with the balance equations
NEW_PROP_TAG(Model);

//! The assembler for the Jacobian matrix
NEW_PROP_TAG(JacobianAssembler);

//! The class storing primary variables plus pseudo primary variables
NEW_PROP_TAG(PrimaryVariables);

//! The number of balance equations.
NEW_PROP_TAG(NumEq);

//! Specifies whether the Jacobian matrix should only be reassembled
//! if the current solution deviates too much from the evaluation point
NEW_PROP_TAG(EnablePartialReassemble);

/*!
 * \brief Specifies whether the update should be done using the line search
 *        method instead of the plain Newton method.
 *
 * Whether this property has any effect depends on whether the line
 * search method is implemented for the actual model's Newton
 * method's update_() method. By default line search is not used.
 */
NEW_PROP_TAG(NewtonEnableLineSearch);

//! Enable Jacobian recycling?
NEW_PROP_TAG(EnableJacobianRecycling);

//! Enable partial reassembly?
NEW_PROP_TAG(EnablePartialReassemble);

// set default values
SET_TYPE_PROP(VcfvNewtonMethod, NewtonMethod, Ewoms::VcfvNewtonMethod<TypeTag>);
SET_TYPE_PROP(VcfvNewtonMethod, NewtonConvergenceWriter,
              Ewoms::VcfvNewtonConvergenceWriter<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonEnableLineSearch, false);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvModel
 * \ingroup Newton
 *
 * \brief A Newton method for models using the VCVF discretization.
 *
 * This class is sufficient for most models which use the VCVF discretization.
 */
template <class TypeTag>
class VcfvNewtonMethod : public NewtonMethod<TypeTag>
{
    typedef Ewoms::NewtonMethod<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    VcfvNewtonMethod(Problem &problem) : ParentType(problem)
    {}

    /*!
     * \brief Register all run-time parameters of the Newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonEnableLineSearch,
                             "Use the line-search update method for the Newton "
                             "method (warning: slow!)");
    }

protected:
    friend class Ewoms::NewtonMethod<TypeTag>;

    /*!
     * \brief Update the relative error of the solution compared to
     *        the previous iteration.
     *
     * The relative error can be seen as a norm of the difference
     * between the current and the next iteration. For the VCFV Newton
     * method, this is the maxiumum of the difference weighted by the
     * primary variable weight.
     *
     * \param uCurrentIter The current iterative solution
     * \param uLastIter The solution of the last iteration
     * \param deltaU The difference between the current and the next solution
     */
    void updateRelError_(const SolutionVector &uCurrentIter,
                         const SolutionVector &uLastIter,
                         const GlobalEqVector &deltaU)
    {
        if (!this->enableRelativeCriterion_() && !enablePartialReassemble_())
            return;

        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        this->relError_ = 0;
        for (int i = 0; i < int(uLastIter.size()); ++i) {
            PrimaryVariables uNewI = uLastIter[i];
            uNewI -= deltaU[i];

            Scalar vertError
                = model_().relativeErrorVertex(i, uLastIter[i], uNewI);
            this->relError_ = std::max(this->relError_, vertError);
        }

        // take the other processes into account
        this->relError_ = this->comm_.max(this->relError_);

        Scalar maxError
            = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxRelativeError);
        if (this->relError_ > maxError)
            OPM_THROW(Opm::NumericalProblem,
                      "Newton: Relative error "
                      << this->relError_
                      << " is larger than maximum allowed error of "
                      << maxError);
    }

    /*!
     * \brief Update the absolute error of the solution compared to
     *        the previous iteration.
     *
     * \param uCurrentIter The current iterative solution
     * \param uLastIter The solution of the last iteration
     * \param deltaU The difference between the current and the next solution
     */
    void updateAbsError_(const SolutionVector &uCurrentIter,
                         const SolutionVector &uLastIter,
                         const GlobalEqVector &deltaU)
    {
        if (!this->enableAbsoluteCriterion_())
            return;
        if (enableLineSearch_())
            // the absolute error has already been calculated by
            // updateLineSearch()
            return;

        // we actually have to do the heavy lifting...
        GlobalEqVector tmp(uLastIter.size());
        this->absError_ = model_().globalResidual(tmp, uCurrentIter);
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void update_(SolutionVector &uCurrentIter, const SolutionVector &uLastIter,
                 const GlobalEqVector &deltaU)
    {
        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(deltaU.two_norm2()))
            OPM_THROW(Opm::NumericalProblem, "Non-finite update!");

        // compute the vertex and element colors for partial reassembly
        if (enablePartialReassemble_()) {
            Scalar minReasmTol = 10 * this->relTolerance_();
            Scalar maxReasmTol = 1e-4;

            // rationale: the newton method has quadratic convergene1
            Scalar reassembleTol = this->relError_ * this->relError_;
            reassembleTol
                = std::max(minReasmTol, std::min(maxReasmTol, reassembleTol));
            // Scalar reassembleTol = minReasmTol;

            model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (enableLineSearch_())
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        else {
            for (unsigned i = 0; i < uLastIter.size(); ++i) {
                uCurrentIter[i] = uLastIter[i];
                uCurrentIter[i] -= deltaU[i];
            }
        }
    }

    /*!
     * \brief Update using the line search algorithm.
     */
    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const GlobalEqVector &deltaU)
    {
        Scalar lambda = 1.0;
        GlobalEqVector tmp(uLastIter.size());

        while (true) {
            for (unsigned i = 0; i < uCurrentIter.size(); ++i) {
                for (int j = 0; j < numEq; ++j) {
                    uCurrentIter[i][j] = uLastIter[i][j] - lambda * deltaU[i][j];
                }
            }

            // calculate the residual of the current solution
            updateAbsError_(uCurrentIter, uLastIter, deltaU);
            if (this->absError_ < this->lastAbsError_ || lambda <= 1.0 / 8) {
                this->endIterMsg() << ", defect " << this->lastAbsError_ << "->"
                                   << this->absError_ << "@lambda=" << lambda;
                return;
            }

            // try with a smaller update
            lambda /= 2.0;
        }
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ end_()
     */
    void failed_()
    {
        ParentType::failed_();

        model_().jacobianAssembler().reassembleAll();
    }

    /*!
     * \brief Called when the Newton method was successful.
     *
     * This method is called _after_ end_()
     */
    void succeeded_()
    {
        ParentType::succeeded_();

        if (enableJacobianRecycling_())
            model_().jacobianAssembler().setMatrixReuseable(true);
        else
            model_().jacobianAssembler().reassembleAll();
    }

    /*!
     * \brief Returns a reference to the problem.
     */
    Model &model_()
    { return ParentType::model(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Model &model_() const
    { return ParentType::model(); }

    /*!
     * \brief Returns true iff the Jacobian assembler uses partial
     *        reassembly.
     */
    bool enablePartialReassemble_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, EnablePartialReassemble); }

    /*!
     * \brief Returns true iff the Jacobian assembler recycles the matrix
     *        possible.
     */
    bool enableJacobianRecycling_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableJacobianRecycling); }

    /*!
     * \brief Returns true iff line search update proceedure should be
     *        used instead of the normal one.
     */
    bool enableLineSearch_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, NewtonEnableLineSearch); }
};
} // namespace Ewoms

#endif
