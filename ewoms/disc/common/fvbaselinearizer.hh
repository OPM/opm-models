// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2015 by Andreas Lauser
  Copyright (C) 2009-2011 by Bernd Flemisch

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
 * \copydoc Ewoms::FvBaseLinearizer
 */
#ifndef EWOMS_FV_BASE_LINEARIZER_HH
#define EWOMS_FV_BASE_LINEARIZER_HH

#include "fvbaseproperties.hh"

#include <ewoms/parallel/gridcommhandles.hh>
#include <ewoms/parallel/threadmanager.hh>
#include <ewoms/parallel/threadedentityiterator.hh>
#include <ewoms/aux/baseauxiliarymodule.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <type_traits>
#include <iostream>
#include <vector>
#include <set>

namespace Ewoms {
// forward declarations
template<class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief The common code for the linearizers of non-linear systems of equations
 *
 * This class assumes that these system of equations to be linearized are stemming from
 * models that use an finite volume scheme for spatial discretization and an Euler
 * scheme for time discretization.
 */
template<class TypeTag>
class FvBaseLinearizer
{
//! \cond SKIP_THIS
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridCommHandleFactory) GridCommHandleFactory;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef GlobalEqVector Vector;
    typedef JacobianMatrix Matrix;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize) };

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    static const bool linearizeNonLocalElements = GET_PROP_VALUE(TypeTag, LinearizeNonLocalElements);

    // copying the linearizer is not a good idea
    FvBaseLinearizer(const FvBaseLinearizer&);
//! \endcond

public:
    /*!
     * \brief The colors of elements and degrees of freedom required
     *        for partial relinearization.
     */
    enum EntityColor {
        /*!
         * Degree of freedom/element that needs to be relinearized
         * because some error is above the tolerance
         */
        Red = 2,

        /*!
         * Degree of freedom/element that needs to be relinearized because a
         * neighboring element/degree of freedom is red
         */
        Yellow = 1,

        /*!
         * Degree of freedom/element that does not need to be relinearized
         */
        Green = 0
    };

    FvBaseLinearizer()
    {
        simulatorPtr_ = 0;

        matrix_ = 0;

        // set relinearization accuracies to 0, so that if partial relinearization of the
        // system of equations is disabled, the relinearization accuracy is always
        // smaller than the current tolerance
        relinearizationAccuracy_ = 0.0;
    }

    ~FvBaseLinearizer()
    {
        delete matrix_;
        auto it = elementCtx_.begin();
        const auto &endIt = elementCtx_.end();
        for (; it != endIt; ++it)
            delete *it;
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian linearizer.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableLinearizationRecycling,
                             "Re-use of the linearized system of equations at the first "
                             "iteration of the next time step");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnablePartialRelinearization,
                             "relinearize only those degrees of freedom that have changed "
                             "'sufficiently' between two Newton iterations");
    }

    /*!
     * \brief Initialize the linearizer.
     *
     * At this point we can assume that all objects in the simulator
     * have been allocated. We cannot assume that they are fully
     * initialized, though.
     *
     * \copydetails Doxygen::simulatorParam
     */
    void init(Simulator& simulator)
    {
        simulatorPtr_ = &simulator;
        delete matrix_; // <- note that this even works for nullpointers!
        matrix_ = 0;
    }

    /*!
     * \brief Causes the Jacobian matrix to be recreated in the next iteration.
     */
    void recreateMatrix()
    {
        delete matrix_; // <- note that this even works for nullpointers!
        matrix_ = 0;
    }

    /*!
     * \brief Linearize the global non-linear system of equations
     *
     * That means that the global Jacobian of the residual is assembled and the residual
     * is evaluated for the current solution.
     *
     * The current state of affairs (esp. the previous and the current solutions) is
     * represented by the model object.
     */
    void linearize()
    {
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (!matrix_)
            initFirstIteration_();

        // we need to store whether the linearization was recycled
        // here because the linearize_ method modifies the
        // reuseLinearization_ attribute!
        bool linearizationReused = reuseLinearization_;

        // store the data required for the end-of-iteration message here because the
        // linearize_() method modifies them for the next iteration...
        int curNumRelin = numTotalElems_ - numGreenElems_;
        Scalar curRelAcc = relinearizationAccuracy_;

        int succeeded;
        try {
            linearize_();
            succeeded = 1;
            succeeded = gridView_().comm().min(succeeded);
        }
        catch (Opm::NumericalIssue &e)
        {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
            succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            OPM_THROW(Opm::NumericalIssue,
                       "A process did not succeed in linearizing the system");
        }

        if (!linearizationReused && enablePartialRelinearization_()) {
            model_().newtonMethod().endIterMsg()
                << ", relinearized " << curNumRelin << " of " << numTotalElems_
                << " elements (" << 100*Scalar(curNumRelin)/numTotalElems_ << "%)"
                << " and achieved a accuracy of " << curRelAcc;
        }
    }

    /*!
     * \brief If linearization recycling is enabled, this method
     *        specifies whether the next call to linearize() just
     *        rescales the storage term or does a full relinearization
     *
     * \param yesno If true, only rescale; else always do a full relinearization.
     */
    void setLinearizationReusable(bool yesno = true)
    {
        if (enableLinearizationRecycling_())
            reuseLinearization_ = yesno;
    }

    /*!
     * \brief If partial relinearization is enabled, this method causes all
     *        elements to be relinearized in the next linearize() call.
     */
    void relinearizeAll()
    {
        // do not reuse the current linearization
        reuseLinearization_ = false;

        // do not use partial relinearization for the next iteration
        relinearizationAccuracy_ = 0.0;
        relinearizationTolerance_ = 0.0;
        numGreenElems_ = 0;
        if (enablePartialRelinearization_()) {
            std::fill(dofError_.begin(), dofError_.end(), 1e100);
            std::fill(dofColor_.begin(), dofColor_.end(), Red);
            std::fill(elementColor_.begin(), elementColor_.end(), Red);
        }

        // set the intensive quantities cache
        auto &model = model_();
        int numGridDof = model.numGridDof();
        for (int timeIdx = 0; timeIdx < historySize; ++timeIdx) {
            for (int dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                model.setIntensiveQuantitiesCacheEntryValidity(dofIdx, timeIdx, false);
            }
        }
    }

    /*!
     * \brief Returns the largest error of a "green" degree of freedom for the
     *        most recent call of the linearize() method.
     *
     * This only has an effect if partial Jacobian relinearization is enabled. If it is
     * disabled, then this method always returns 0.
     *
     * This returns the _actual_ error computed as seen by computeColors(), not the
     * tolerance which it was given.
     */
    Scalar relinearizationAccuracy() const
    { return relinearizationAccuracy_; }

    /*!
     * \brief The maximum deflection seen for any DOF after an Newton update.
     */
    Scalar maxDofError() const
    { return maxDofError_; }

    /*!
     * \brief Update the distance where the non-linear system was
     *        originally insistently linearized and the point where it
     *        will be linerized the next time.
     *
     * This only has an effect if partial relinearize is enabled.
     *
     * \param uDelta The negative difference between the last and the next iterative solution.
     * \param resid The residual (right-hand side) for the current Newton iteration
     */
    void updateRelinearizationErrors(const GlobalEqVector &uDelta, const GlobalEqVector &resid)
    {
        if (!enablePartialRelinearization_())
            return;

        unsigned numGridDof = model_().numGridDof();

        // update the vector which stores the error for partial
        // relinearization for each degree of freedom
        maxDofError_ = 0;
        for (unsigned globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx) {
            if (!model_().isLocalDof(globalDofIdx)) {
                // ignore degrees of freedom of overlap and ghost degrees of freedom
                dofError_[globalDofIdx] = 0;
                continue;
            }

            // compute the error of the solution at a given DOF. we use the weighted
            // magnitude of the deviation between the solution for two consecutive
            // iterations.
            const auto &d = uDelta[globalDofIdx];
            Scalar distRel = 0;
            for (unsigned pvIdx = 0; pvIdx < d.size(); ++pvIdx) {
                Scalar tmp = std::abs(d[pvIdx]*model_().primaryVarWeight(globalDofIdx, pvIdx));
                distRel = std::max(distRel, tmp);
            }
            Valgrind::CheckDefined(distRel);
            dofError_[globalDofIdx] += distRel;
            maxDofError_ = std::max(maxDofError_, dofError_[globalDofIdx]);
        }
    }

    /*!
     * \brief Ensure that a given degree of freedom is relinarized in the next iteration.
     *
     * Calling this method usually means that the interpretation of the primary variables
     * for the DOF has changed.
     */
    void markDofRed(int dofIdx)
    {
        if (enablePartialRelinearization_())
            dofError_[dofIdx] = 1e100;
        this->model_().setIntensiveQuantitiesCacheEntryValidity(dofIdx, /*timeIdx=*/0, false);
    }

    /*!
     * \brief Returns the "relinearization color" of a degree of freedom
     */
    EntityColor dofColor(int dofIdx) const
    {
        if (!enablePartialRelinearization_())
            return Red;
        return dofColor_[dofIdx];
    }

    /*!
     * \brief Returns the "relinearization color" of an element
     */
    EntityColor elementColor(int elemIdx) const
    {
        if (!enablePartialRelinearization_())
            return Red;
        return elementColor_[elemIdx];
    }

    /*!
     * \brief Returns the maximum error for which a degree of freedom is not relinearized.
     */
    Scalar relinearizationTolerance() const
    { return relinearizationTolerance_; }

    /*!
     * \brief Sets the maximum error for which a degree of freedom is not relinearized.
     */
    void setRelinearizationTolerance(Scalar tolerance)
    { relinearizationTolerance_ = tolerance; }

    /*!
     * \brief Returns the error for a given degree of freedom after the last iteration.
     */
    Scalar dofError(int dofIdx) const
    { return dofError_[dofIdx]; }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const Matrix &matrix() const
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const GlobalEqVector &residual() const
    { return residual_; }

private:
    /*!
     * \brief Determine the colors of the degrees of freedom and of
     *        the elements for partial re-linarization for a given
     *        tolerance.
     */
    void computeColors_()
    {
        if (!enablePartialRelinearization_())
            return;

        // mark the red degrees of freedom and update the tolerance of
        // the linearization which actually will get achieved
        for (unsigned dofIdx = 0; dofIdx < dofColor_.size(); ++dofIdx) {
            if (dofError_[dofIdx] > relinearizationTolerance_) {
                // mark the degree of freedom 'red' if discrepancy is
                // larger than the given tolerance
                dofColor_[dofIdx] = Red;
                dofError_[dofIdx] = 0.0;
            }
            else if (!model_().isLocalDof(dofIdx)) {
                // always mark non-local degrees of freedom as red
                dofColor_[dofIdx] = Red;
                dofError_[dofIdx] = 0.0;
            }
            else
                dofColor_[dofIdx] = Green;
        }

        // mark all elements that connect to red ones as yellow
        int numGridDof = this->model_().numGridDof();
        for (int dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
            if (dofColor_[dofIdx] != Red)
                continue;

            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[dofIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[dofIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int dof2Idx = colIt.index();
                if (dof2Idx >= numGridDof)
                    break; // auxiliary equation

                if (dofColor_[dof2Idx] != Red)
                    dofColor_[dof2Idx] = Yellow;
            }
        }

        // calculate the element colors from the DOF colors. for the element centered
        // finite volume discretization, we can just copy the values. For other
        // discretizations, we need to loop over the grid and check if an element's
        // stencil contains a non-green DOFs.
        if (std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value)
            elementColor_ = dofColor_;
        else {
            const auto& gridView = this->model_().gridView();
            const auto& elemMapper = this->model_().elementMapper();
            Stencil stencil(gridView);
            ElementIterator elemIt = gridView.template begin</*codim=*/0>();
            const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;
                stencil.updateTopology(elem);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
                int elemIdx = elemMapper.index(elem);
#else
                int elemIdx = elemMapper.map(elem);
#endif
                int numElemDof = stencil.numDof();

                // determine the element color: it is red if any of the DOFs in its
                // stencil are non-green, else it is green.
                elementColor_[elemIdx] = Green;
                for (int elemDofIdx = 0; elemDofIdx < numElemDof; ++elemDofIdx) {
                    int dofIdx = stencil.globalSpaceIndex(elemDofIdx);
                    if (dofColor_[dofIdx] != Green) {
                        elementColor_[elemIdx] = Red;
                        break;
                    }
                }
            }
        }

        // gather statistics
        relinearizationAccuracy_ = 0;
        numGreenElems_ = 0;
        for (unsigned dofIdx = 0; dofIdx < dofColor_.size(); ++dofIdx) {
            if (dofColor_[dofIdx] != Red) {
                relinearizationAccuracy_ =
                    std::max(relinearizationAccuracy_, dofError_[dofIdx]);

                if (dofColor_[dofIdx] == Green)
                    ++numGreenElems_;
                continue;
            }
        }

        relinearizationAccuracy_ = gridView_().comm().max(relinearizationAccuracy_);
        numGreenElems_ = gridView_().comm().sum(numGreenElems_);
    }

    static bool enableLinearizationRecycling_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableLinearizationRecycling); }
    static bool enablePartialRelinearization_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnablePartialRelinearization); }

    Simulator &simulator_()
    { return *simulatorPtr_; }
    const Simulator &simulator_() const
    { return *simulatorPtr_; }

    Problem &problem_()
    { return simulator_().problem(); }
    const Problem &problem_() const
    { return simulator_().problem(); }

    Model &model_()
    { return simulator_().model(); }
    const Model &model_() const
    { return simulator_().model(); }

    const GridView &gridView_() const
    { return problem_().gridView(); }

    const ElementMapper &elementMapper_() const
    { return model_().elementMapper(); }

    const DofMapper &dofMapper_() const
    { return model_().dofMapper(); }

    void initFirstIteration_()
    {
        // calculate the number of DOFs and elements for the local process and for the
        // whole simulation
        int numGridDof = model_().numGridDof();
        int numAllDof =  model_().numTotalDof();
        int numElems = gridView_().size(/*codim=*/0);
        numTotalElems_ = gridView_().comm().sum(numElems);

        // initialize the BCRS matrix for the Jacobian
        createMatrix_();

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        residual_.resize(numAllDof);
        residual_ = 0;

        // create the per-thread context objects
        elementCtx_.resize(ThreadManager::maxThreads());
        for (int threadId = 0; threadId != ThreadManager::maxThreads(); ++ threadId)
            elementCtx_[threadId] = new ElementContext(simulator_());

        // initialize the storage part of the Jacobian matrix. Since we only need this if
        // linearization recycling is enabled, we do not waste space if it is disabled
        if (enableLinearizationRecycling_()) {
            storageJacobian_.resize(numGridDof);
            storageTerm_.resize(numGridDof);
       }

        // initialize data needed for partial relinearization
        reuseLinearization_ = false;
        if (enablePartialRelinearization_()) {
            dofColor_.resize(numGridDof);
            dofError_.resize(numGridDof);
            elementColor_.resize(numElems);
            relinearizeAll();
        }
    }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int numAllDof =  model_().numTotalDof();

        // allocate raw matrix
        matrix_ = new Matrix(numAllDof, numAllDof, Matrix::random);

        Stencil stencil(gridView_());

        // for the main model, find out the global indices of the neighboring degrees of
        // freedom of each primary degree of freedom
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(numAllDof);
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element &elem = *elemIt;
            stencil.update(elem);

            for (int primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                int myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (int dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    int neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    neighbors[myIdx].insert(neighborIdx);
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        const auto& model = model_();
        int numAuxMod = model.numAuxiliaryModules();
        for (int auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model.auxiliaryModule(auxModIdx)->addNeighbors(neighbors);

        // allocate space for the rows of the matrix
        for (int dofIdx = 0; dofIdx < numAllDof; ++ dofIdx)
            matrix_->setrowsize(dofIdx, neighbors[dofIdx].size());
        matrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (int dofIdx = 0; dofIdx < numAllDof; ++ dofIdx) {
            typename NeighborSet::iterator nIt = neighbors[dofIdx].begin();
            typename NeighborSet::iterator nEndIt = neighbors[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                matrix_->addindex(dofIdx, *nIt);
        }
        matrix_->endindices();
    }

    // reset the global linear system of equations. if partial
    // relinearization is enabled, this means that the Jacobian matrix
    // must only be erased partially!
    void resetSystem_()
    {
        size_t numGridDof = model_().numGridDof();
        size_t numTotalDof = model_().numTotalDof();

        if (!enablePartialRelinearization_()) {
            // if partial re-linearization is not enabled, we can just reset everything!
            residual_ = 0.0;
            (*matrix_) = 0;

            // reset the parts needed for linearization recycling
            if (enableLinearizationRecycling_()) {
                for (unsigned dofIdx=0; dofIdx < numGridDof; ++ dofIdx) {
                    storageTerm_[dofIdx] = 0.0;
                    storageJacobian_[dofIdx] = 0.0;
                }
            }

            return;
        }

        // always reset the right hand side completely
        residual_ = 0.0;
        if (enableLinearizationRecycling_())
            for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx)
                storageTerm_[dofIdx] = 0.0;

        // reset the rows in the Jacobian which correspond to DOFs of auxiliary equations
        for (unsigned dofIdx = numGridDof; dofIdx < numTotalDof; ++dofIdx) {
            // reset the row of the Jacobian matrix
            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[dofIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[dofIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                (*colIt) = 0.0;
            }
        }

        // partially reset the current linearization for rows corresponding to grid DOFs
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
            if (dofColor(dofIdx) == Green) {
                // for green DOFs we keep the left hand side of the current linearization
                // except for the entries of the Jacobian matrix which connect the green
                // DOF with auxiliary equations. the implementation of this this is
                // slightly hacky because it implicitly assumes that all auxiliary DOFs
                // are placed after the grid DOFs and that auxiliary DOFs don't change
                // any entries of grid DOFs. For auxiliary equations where this is not
                // the case, we would need to mark all grid DOFs which have a connection
                // to an auxiliary DOF as red...
                typedef typename JacobianMatrix::ColIterator ColIterator;
                ColIterator colIt = (*matrix_)[dofIdx].beforeEnd();
                const ColIterator colBeginIt = (*matrix_)[dofIdx].beforeBegin();
                for (; colIt != colBeginIt; -- colIt) {
                    if (colIt.index() < numGridDof)
                        // the column corresponds to a grid DOF. we have reset
                        // everything.
                        break;

                    // the column corresponds to an auxiliary DOF
                    (*colIt) = 0.0;
                }
            }
            else {
                // red or yellow DOF
                if (enableLinearizationRecycling_())
                    storageJacobian_[dofIdx] = 0.0;

                // reset all entries in the row of the Jacobian which connect two non-green
                // degrees of freedom
                typedef typename JacobianMatrix::ColIterator ColIterator;
                ColIterator colIt = (*matrix_)[dofIdx].begin();
                const ColIterator &colEndIt = (*matrix_)[dofIdx].end();
                for (; colIt != colEndIt; ++colIt) {
                    if (colIt.index() >= numGridDof || dofColor_[colIt.index()] != Green)
                        // the entry either corresponds to an auxiliary DOF or to a non-green
                        // grid DOF.
                        (*colIt) = 0.0;
                }
            }
        }
    }

    // linearize the whole system
    void linearize_()
    {
        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().simulator().timeStepSize();
        if (reuseLinearization_) {
            int numGridDof = model_().numGridDof();
            for (int dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                // use the flux term plus the source term as the residual: the numerator
                // in the d(storage)/dt term is 0 for the first iteration of a time step
                // because the initial guess for the next solution is the value of the
                // last time step.
                residual_[dofIdx] -= storageTerm_[dofIdx];

                // rescale the contributions of the storage term to the Jacobian matrix
                MatrixBlock &J_ii = (*matrix_)[dofIdx][dofIdx];

                J_ii -= storageJacobian_[dofIdx];
                storageJacobian_[dofIdx] *= oldDt_/curDt;
                J_ii += storageJacobian_[dofIdx];
            }

            reuseLinearization_ = false;
            oldDt_ = curDt;

            linearizeAuxiliaryEquations_();

            problem_().newtonMethod().endIterMsg()
                << ", linear system of equations rescaled using previous time step";
            return;
        }

        oldDt_ = curDt;
        computeColors_();
        resetSystem_();

        // relinearize the elements...
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(gridView_());
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementIterator elemIt = gridView_().template begin</*codim=*/0>();
            for (threadedElemIt.beginParallel(elemIt);
                 !threadedElemIt.isFinished(elemIt);
                 threadedElemIt.increment(elemIt))
            {
                const Element &elem = *elemIt;

                if (!linearizeNonLocalElements && elem.partitionType() != Dune::InteriorEntity)
                    continue;

                linearizeElement_(elem);
            }
        }

        linearizeAuxiliaryEquations_();
    }

    // linearize an element in the interior of the process' grid partition
    void linearizeElement_(const Element &elem)
    {
        if (enablePartialRelinearization_()) {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int globalElemIdx = model_().elementMapper().index(elem);
#else
            int globalElemIdx = model_().elementMapper().map(elem);
#endif
            if (elementColor(globalElemIdx) == Green) {
                linearizeGreenElement_(elem);
                return;
            }
        }

        int threadId = ThreadManager::threadId();

        ElementContext *elementCtx = elementCtx_[threadId];
        auto &localLinearizer = model_().localLinearizer(threadId);

        // the actual work of linearization is done by the local linearizer class
        elementCtx->updateAll(elem);
        localLinearizer.linearize(*elementCtx);

        // update the right hand side and the Jacobian matrix
        ScopedLock addLock(globalMatrixMutex_);
        int numPrimaryDof = elementCtx->numPrimaryDof(/*timeIdx=*/0);
        for (int primaryDofIdx = 0; primaryDofIdx < numPrimaryDof; ++ primaryDofIdx) {
            int globI = elementCtx->globalSpaceIndex(/*spaceIdx=*/primaryDofIdx, /*timeIdx=*/0);

            // we only need to update the Jacobian matrix for entries which connect two
            // non-green DOFs. if the row DOF corresponds to a green one, we can skip the
            // whole row...
            if (dofColor(globI) == Green)
                continue;

            // update the right hand side
            residual_[globI] += localLinearizer.residual(primaryDofIdx);

            if (enableLinearizationRecycling_()) {
                storageTerm_[globI] += localLinearizer.residualStorage(primaryDofIdx);
                storageJacobian_[globI] += localLinearizer.jacobianStorage(primaryDofIdx);
            }

            // update the global Jacobian matrix
            for (int dofIdx = 0; dofIdx < elementCtx->numDof(/*timeIdx=*/0); ++ dofIdx) {
                int globJ = elementCtx->globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);


                // only update the jacobian matrix for non-green degrees of freedom
                if (dofColor(globJ) == Green)
                    continue;

                (*matrix_)[globJ][globI] += localLinearizer.jacobian(dofIdx, primaryDofIdx);
            }
        }
        addLock.unlock();
    }

    // "linearize" a green element. green elements only get the
    // residual updated, but the Jacobian is left alone...
    void linearizeGreenElement_(const Element &elem)
    {
        int threadId = ThreadManager::threadId();
        ElementContext *elementCtx = elementCtx_[threadId];

        elementCtx->updateAll(elem);
        auto& localResidual = model_().localResidual(threadId);
        localResidual.eval(*elementCtx);

        ScopedLock addLock(globalMatrixMutex_);
        for (int dofIdx=0; dofIdx < elementCtx->numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
            int globI = elementCtx->globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            // update the right hand side
            for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx)
                residual_[globI][eqIdx] += Toolbox::value(localResidual.residual(dofIdx)[eqIdx]);
            if (enableLinearizationRecycling_())
                for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx)
                    storageTerm_[globI][eqIdx] += Toolbox::value(localResidual.storageTerm(dofIdx)[eqIdx]);
        }
        addLock.unlock();
    }

    void linearizeAuxiliaryEquations_()
    {
        auto& model = model_();
        for (unsigned auxModIdx = 0; auxModIdx < model.numAuxiliaryModules(); ++auxModIdx)
            model.auxiliaryModule(auxModIdx)->linearize(*matrix_, residual_);
    }

    Simulator *simulatorPtr_;
    std::vector<ElementContext*> elementCtx_;

    // the jacobian matrix
    Matrix *matrix_;
    // the right-hand side
    GlobalEqVector residual_;

    // attributes required for jacobian matrix recycling
    bool reuseLinearization_;
    // The storage part of the local Jacobian
    std::vector<MatrixBlock> storageJacobian_;
    std::vector<VectorBlock> storageTerm_;
    // time step size used for the last linearization
    Scalar oldDt_;

    // data required for partial relinearization
    std::vector<EntityColor> dofColor_;
    std::vector<Scalar> dofError_;
    std::vector<EntityColor> elementColor_;

    int numTotalElems_;
    int numGreenElems_;

    Scalar relinearizationTolerance_;
    Scalar relinearizationAccuracy_;
    Scalar maxDofError_;

    OmpMutex globalMatrixMutex_;
};

} // namespace Ewoms

#endif
