/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2011-2012 by Bernd Flemisch
  Copyright (C) 2011 by Klaus Mosthaf

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
 * \copydoc Ewoms::FvBaseLocalJacobian
 */
#ifndef EWOMS_FV_BASE_LOCAL_JACOBIAN_HH
#define EWOMS_FV_BASE_LOCAL_JACOBIAN_HH

#include "fvbaseproperties.hh"

#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {

/*!
 * \ingroup Discretization
 *
 * \brief Calculates the Jacobian of the local residual for ECFV discretizations
 *
 * The default behavior is to use numeric differentiation, i.e.
 * forward or backward differences (2nd order), or central
 * differences (3rd order). The method used is determined by the
 * "NumericDifferenceMethod" property:
 *
 * - if the value of this property is smaller than 0, backward
 *   differences are used, i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - if the value of this property is 0, central
 *   differences are used, i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx
 *          \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward
 *   differences are used, i.e.:
 *   \f[
 *     \frac{\partial f(x)}{\partial x} \approx
 *          \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
 * is the value of a sub-control volume's primary variable at the
 * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
 */
template<class TypeTag>
class FvBaseLocalJacobian
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),

        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize)
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;
    typedef Dune::BlockVector<MatrixBlock> LocalStorageMatrix;

    // copying a local jacobian is not a good idea
    FvBaseLocalJacobian(const FvBaseLocalJacobian &);

public:
    FvBaseLocalJacobian()
    {
        internalElemContext_ = 0;
        newtonTolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonRelativeTolerance);
    }

    ~FvBaseLocalJacobian()
    {
        delete internalElemContext_;
    }

    /*!
     * \brief Register all run-time parameters for the local jacobian.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, NumericDifferenceMethod,
                             "The method used for numeric differentiation (-1: backward "
                             "differences, 0: central differences, 1: forward differences)");
    }

    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param prob The problem which we want to simulate.
     */
    void init(Problem &prob)
    {
        problemPtr_ = &prob;
        modelPtr_ = &prob.model();
        internalElemContext_ = new ElementContext(prob);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * This assembles the 'grad f(x^k)' and 'f(x^k)' part of the
     * newton update for a single element.
     *
     * \param element The grid element for which the local residual
     *                and its gradient should be calculated.
     */
    void assemble(const Element &element)
    {
        internalElemContext_->updateAll(element);

        assemble(*internalElemContext_);
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect, given all secondary variables for the element.
     *
     * After calling this method the ElementContext are in undefined
     * state, so do not use it anymore!
     *
     * \param elemCtx The element execution context for which the
     *                local residual and its gradient should be
     *                calculated.
     */
    void assemble(ElementContext &elemCtx)
    {
        // update the hints for the element's volume variables
        for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            int globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            for (int timeIdx = 0; timeIdx < historySize; ++timeIdx)
                model_().setHint(elemCtx.volVars(dofIdx, timeIdx),
                                 globalIdx,
                                 timeIdx);
        }

        // update the weights of the primary variables using the
        // current element variables
        model_().updatePVWeights(elemCtx);

        resize_(elemCtx);
        reset_(elemCtx);

        // calculate the local residual
        localResidual_.eval(residual_, residualStorage_, elemCtx);

        // save all flux variables calculated using the unmodified
        // primary variables. This automatically makes these flux
        // variables the evaluation point.
        elemCtx.saveFluxVars();

        // calculate the local jacobian matrix
        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < numDof; dofIdx++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().evalPartialDerivative_(elemCtx, dofIdx, pvIdx);

                // update the local stiffness matrix with the current
                // partial derivatives
                updateLocalJacobian_(elemCtx, dofIdx, pvIdx);
            }
        }

        // restore flux variables.
        //elemCtx.restoreScvfVars(); // not necessary
    }

    /*!
     * \brief Returns the unweighted epsilon value used to calculate
     *        the local derivatives
     */
    static Scalar baseEpsilon()
    { return GET_PROP_VALUE(TypeTag, BaseEpsilon); }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemCtx The element execution context for which the
     *                local residual and its gradient should be
     *                calculated.
     * \param dofIdx     The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon(const ElementContext &elemCtx,
                          int dofIdx,
                          int pvIdx) const
    {
        int globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        Scalar pvWeight = elemCtx.model().primaryVarWeight(globalIdx, pvIdx);
        Valgrind::CheckDefined(pvWeight);

        return baseEpsilon()/pvWeight;
    }

    /*!
     * \brief Return reference to the local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Return reference to the local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns the local Jacobian matrix of the residual of a sub-control volume.
     *
     * \param domainScvIdx The local index of the sub control volume
     *                     which contains the independents
     * \param rangeScvIdx The local index of the sub control volume
     *                    which contains the local residual
     */
    const MatrixBlock &jacobian(int domainScvIdx, int rangeScvIdx) const
    { return jacobian_[domainScvIdx][rangeScvIdx]; }

    /*!
     * \brief Returns the local Jacobian matrix the storage term of a sub-control volume.
     *
     * \param dofIdx The local index of sub control volume
     */
    const MatrixBlock &jacobianStorage(int dofIdx) const
    { return jacobianStorage_[dofIdx]; }

    /*!
     * \brief Returns the local residual of a sub-control volume.
     *
     * \param dofIdx The local index of the sub control volume
     */
    const VectorBlock &residual(int dofIdx) const
    { return residual_[dofIdx]; }

    /*!
     * \brief Returns the local storage term of a sub-control volume.
     *
     * \param dofIdx The local index of the sub control volume
     */
    const VectorBlock &residualStorage(int dofIdx) const
    { return residualStorage_[dofIdx]; }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return *modelPtr_; }

    /*!
     * \brief Returns the numeric difference method which is applied.
     */
    static int numericDifferenceMethod_()
    { return EWOMS_GET_PARAM(TypeTag, int, NumericDifferenceMethod); }

    /*!
     * \brief Resize all internal attributes to the size of the
     *        element.
     */
    void resize_(const ElementContext &elemCtx)
    {
        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);

        jacobian_.setSize(numPrimaryDof, numDof);
        jacobianStorage_.resize(numPrimaryDof);

        residual_.resize(numPrimaryDof);
        residualStorage_.resize(numPrimaryDof);

        derivResidual_.resize(numPrimaryDof);
        derivStorage_.resize(numPrimaryDof);
    }

    /*!
     * \brief Reset the all relevant internal attributes to 0
     */
    void reset_(const ElementContext &elemCtx)
    {
        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < numPrimaryDof; ++ dofIdx) {
            residual_[dofIdx] = 0.0;
            residualStorage_[dofIdx] = 0.0;

            jacobianStorage_[dofIdx] = 0.0;
            for (int dof2Idx = 0; dof2Idx < numDof; ++ dof2Idx) {
                jacobian_[dofIdx][dof2Idx] = 0.0;
            }
        }
    }

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than numerical differentiation is available.
     *
     * The default implementation of this method uses numeric
     * differentiation, i.e. forward or backward differences (2nd
     * order), or central differences (3rd order). The method used is
     * determined by the "NumericDifferenceMethod" property:
     *
     * - if the value of this property is smaller than 0, backward
     *   differences are used, i.e.:
     *   \f[
         \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
     *   \f]
     *
     * - if the value of this property is 0, central
     *   differences are used, i.e.:
     *   \f[
     *     \frac{\partial f(x)}{\partial x} \approx
     *          \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
     *   \f]
     *
     * - if the value of this property is larger than 0, forward
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
     *   \f]
     *
     * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
     * is the value of a sub-control volume's primary variable at the
     * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
     *
     * \param elemCtx The element context for which the local partial
     *                derivative ought to be calculated
     * \param dofIdx The sub-control volume index of the current
     *               finite element for which the partial derivative
     *               ought to be calculated
     * \param pvIdx The index of the primary variable at the dofIdx'
     *              sub-control volume of the current finite element
     *              for which the partial derivative ought to be
     *              calculated
     */
    void evalPartialDerivative_(ElementContext &elemCtx,
                                int dofIdx,
                                int pvIdx)
    {
        // save all quantities which depend on the specified primary
        // variable at the given sub control volume
        elemCtx.saveVolVars(dofIdx);

        PrimaryVariables priVars(elemCtx.primaryVars(dofIdx, /*timeIdx=*/0));
        Scalar eps = asImp_().numericEpsilon(elemCtx, dofIdx, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_() >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the residual
            elemCtx.updateVolVars(priVars, dofIdx, /*timeIdx=*/0);
            elemCtx.updateAllFluxVars();
            localResidual_.eval(derivResidual_, derivStorage_, elemCtx);
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ = residual_;
            derivStorage_ = residualStorage_;
        }

        if (numericDifferenceMethod_() <= 0) {
            // we are not using forward differences, i.e. we don't
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // calculate residual again, this time we use the local
            // residual's internal storage.
            elemCtx.updateVolVars(priVars, dofIdx, /*timeIdx=*/0);
            elemCtx.updateAllFluxVars();
            localResidual_.eval(elemCtx);

            derivResidual_ -= localResidual_.residual();
            derivStorage_ -= localResidual_.storageTerm();
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            derivResidual_ -= residual_;
            derivStorage_ -= residualStorage_;
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        derivResidual_ /= delta;
        derivStorage_ /= delta;

        // restore the original state of the element's volume
        // variables
        elemCtx.restoreVolVars(dofIdx);

#ifndef NDEBUG
        for (unsigned i = 0; i < derivResidual_.size(); ++i)
            Valgrind::CheckDefined(derivResidual_[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex 'dofIdx' .
     */
    void updateLocalJacobian_(const ElementContext &elemCtx,
                              int dofIdx,
                              int pvIdx)
    {
        // store the derivative of the storage term
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        if (dofIdx < numPrimaryDof) {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                jacobianStorage_[dofIdx][eqIdx][pvIdx] = derivStorage_[dofIdx][eqIdx];
            }
        }

        for (int primaryDofIdx = 0; primaryDofIdx < numPrimaryDof; primaryDofIdx++)
        {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[primaryDofIdx][dofIdx][eqIdx][pvIdx] is the rate of
                // change of the residual of equation 'eqIdx' at
                // vertex 'primaryDofIdx' depending on the primary variable
                // 'pvIdx' at vertex 'dofIdx'.
                jacobian_[primaryDofIdx][dofIdx][eqIdx][pvIdx] =
                    derivResidual_[primaryDofIdx][eqIdx];
                Valgrind::CheckDefined(jacobian_[primaryDofIdx][dofIdx][eqIdx][pvIdx]);
            }
        }
    }

    Problem *problemPtr_;
    Model *modelPtr_;

    static Scalar newtonTolerance_;
    ElementContext *internalElemContext_;

    LocalBlockMatrix jacobian_;
    LocalStorageMatrix jacobianStorage_;

    LocalBlockVector residual_;
    LocalBlockVector residualStorage_;

    LocalBlockVector derivResidual_;
    LocalBlockVector derivStorage_;

    LocalResidual localResidual_;
};

template <class TypeTag>
typename FvBaseLocalJacobian<TypeTag>::Scalar FvBaseLocalJacobian<TypeTag>::newtonTolerance_;
} // namespace Ewoms

#endif
