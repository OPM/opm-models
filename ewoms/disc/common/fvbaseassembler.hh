/*
  Copyright (C) 2010-2013 by Andreas Lauser
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
 * \copydoc Ewoms::FvBaseAssembler
 */
#ifndef EWOMS_FV_BASE_ASSEMBLER_HH
#define EWOMS_FV_BASE_ASSEMBLER_HH

#include "fvbaseproperties.hh"

#include <ewoms/parallel/gridcommhandles.hh>
#include <ewoms/parallel/threadmanager.hh>
#include <ewoms/parallel/threadedentityiterator.hh>
#include <ewoms/aux/baseauxiliarymodule.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>
#include <vector>
#include <set>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The common code for the assemblers of the global Jacobian matrix of models
 *        using an implicit finite volumes discretization scheme.
 */
template<class TypeTag>
class FvBaseAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridCommHandleFactory) GridCommHandleFactory;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef GlobalEqVector Vector;
    typedef JacobianMatrix Matrix;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    // copying the jacobian assembler is not a good idea
    FvBaseAssembler(const FvBaseAssembler &);

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
        Red = 0,

        /*!
         * Degree of freedom/element that needs to be relinearized because a
         * neighboring element/degree of freedom is red
         */
        Yellow = 1,

        /*!
         * Yellow degree of freedom has only non-green neighbor elements.
         *
         * This means that its error is below the tolerance, but its
         * defect can be linearized without any additional cost. This
         * is just an "internal" color which is not used ouside of the
         * jacobian assembler.
         */
        Orange = 2,

        /*!
         * Degree of freedom/element that does not need to be relinearized
         */
        Green = 3
    };

    FvBaseAssembler()
    {
        simulatorPtr_ = 0;

        matrix_ = 0;

        // set relinearization accuracy to 0, so that if partial
        // relinearization of the system of equations is disabled, the
        // relinearization accuracy is always smaller than the current
        // tolerance
        relinearizationAccuracy_ = 0.0;

    }

    ~FvBaseAssembler()
    {
        delete matrix_;
        auto it = elementCtx_.begin();
        const auto &endIt = elementCtx_.end();
        for (; it != endIt; ++it)
            delete *it;
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian assembler.
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
     * \brief Initialize the jacobian assembler.
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
     * \brief Assemble the global Jacobian of the residual and the
     *        residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        // we defer the initialization of the Jacobian matrix until here because the
        // auxiliary modules usually assume the problem, model and grid to be fully
        // initialized...
        if (!matrix_)
            initFirstIteration_();

        // we need to store whether the linearization was recycled
        // here because the assemble_ method modifies the
        // reuseLinearization_ attribute!
        bool linearizationReused = reuseLinearization_;

        int succeeded;
        try {
            assemble_();
            succeeded = 1;
            succeeded = gridView_().comm().min(succeeded);
        }
        catch (Opm::NumericalProblem &e)
        {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
            succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            OPM_THROW(Opm::NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }

        if (!linearizationReused && enablePartialRelinearization_()) {
            greenElems_ = gridView_().comm().sum(greenElems_);
            relinearizationAccuracy_ = gridView_().comm().max(nextRelinearizationAccuracy_);

            model_().newtonMethod().endIterMsg()
                << ", relinearized " << totalElems_ - greenElems_ << " of " << totalElems_
                << " elements (" << 100*Scalar(totalElems_ - greenElems_)/totalElems_ << "%)"
                << " and achieved an accuracy of " << relinearizationAccuracy_;
        }

        // reset all degree of freedom colors to green
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            dofColor_[i] = Green;
        }
    }

    /*!
     * \brief If linearization recycling is enabled, this method
     *        specifies whether the next call to assemble() just
     *        rescales the storage term or does a full relinearization
     *
     * \param yesno If true, only rescale; else do full Jacobian assembly.
     */
    void setLinearizationReusable(bool yesno = true)
    {
        if (enableLinearizationRecycling_())
            reuseLinearization_ = yesno;
    }

    /*!
     * \brief If partial relinearization is enabled, this method causes all
     *        elements to be relinearized in the next assemble() call.
     */
    void relinearizeAll()
    {
        // do not reuse the current linearization
        reuseLinearization_ = false;

        // do not use partial relinearization for the next iteration
        nextRelinearizationAccuracy_ = 0.0;
        if (enablePartialRelinearization_()) {
            std::fill(dofError_.begin(), dofError_.end(), 0.0);
            std::fill(dofColor_.begin(), dofColor_.end(), Red);
            std::fill(elementColor_.begin(), elementColor_.end(), Red);
        }
    }

    /*!
     * \brief Returns the largest error of a "green" degree of freedom
     *        for the most recent call of the assemble() method.
     *
     * This only has an effect if partial Jacobian relinearization is
     * enabled. If it is disabled, then this method always returns 0.
     *
     * This returns the _actual_ error computed as seen by
     * computeColors(), not the tolerance which it was given.
     */
    Scalar relinearizationAccuracy() const
    { return relinearizationAccuracy_; }

    /*!
     * \brief Update the distance where the non-linear system was
     *        originally insistently linearized and the point where it
     *        will be linerized the next time.
     *
     * This only has an effect if partial relinearize is enabled.
     *
     * \param u The iterative solution of the last iteration of the Newton method
     * \param uDelta The negative difference between the last and the next iterative solution.
     */
    void updateDiscrepancy(const GlobalEqVector &previousResid)
    {
        if (!enablePartialRelinearization_())
            return;

        unsigned numGridDof = model_().numGridDof();

        // update the vector which stores the error for partial
        // relinearization for each degree of freedom
        for (unsigned globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx) {
            if (model_().dofTotalVolume(globalDofIdx) <= 0) {
                // ignore degrees of freedom of overlap and ghost elements
                dofError_[globalDofIdx] = 0;
                continue;
            }

            // we need to add the distance the solution was moved for
            // this degree of freedom
            const auto &r = previousResid[globalDofIdx];
            Scalar dist = 0;
            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx)
                dist = std::max(dist, std::abs(r[eqIdx]*model_().eqWeight(globalDofIdx, eqIdx)));
            dofError_[globalDofIdx] = dist;
        }

    }

    /*!
     * \brief Force to relinearize a given degree of freedom the next
     *        time the assemble() method is called.
     *
     * \param globalDofIdx The global index of the degree of freedom
     *                     which ought to be red.
     */
    void markDofRed(int globalDofIdx)
    {
        if (!enablePartialRelinearization_())
            return;

        dofColor_[globalDofIdx] = Red;
    }

    /*!
     * \brief Determine the colors of the degrees of freedom and of
     *        the elements for partial re-linarization for a given
     *        tolerance.
     *
     * The following approach is used:
     *
     * - Set all degrees of freedom and elements to 'green'
     * - Mark all degrees of freedom as 'red' which exhibit an error
     *   above the tolerance
     * - Mark all elements which contain 'red' degrees of freedom as 'red'
     * - Mark all degrees of freedom which are not 'red' and are part of a
     *   'red' element as 'yellow'
     * - Mark all elements which are not 'red' and contain a
     *   'yellow' degree of freedom as 'yellow'
     *
     * \param tolerance The error below which a degree of freedom
     *                  won't be relinearized.
     */
    void computeColors(Scalar tolerance)
    {
        if (!enablePartialRelinearization_())
            return;

        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();

        // mark the red degrees of freedom and update the tolerance of
        // the linearization which actually will get achieved
        nextRelinearizationAccuracy_ = 0;
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            if (dofError_[i] > tolerance)
                // mark the degree of freedom 'red' if discrepancy is
                // larger than the given tolerance
                dofColor_[i] = Red;
            else
                nextRelinearizationAccuracy_ =
                    std::max(nextRelinearizationAccuracy_, dofError_[i]);
        }

        Stencil stencil(gridView_());
        // Mark all red elements
        for (; elemIt != elemEndIt; ++elemIt) {
            stencil.update(*elemIt);

            // find out whether the current element contains a red
            // degree of freedom
            bool isRed = false;
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                if (dofColor_[globalIdx] == Red) {
                    isRed = true;
                    break;
                }
            }

            // if yes, the element color is also red, else it is not
            // red, i.e. green for the mean time
            int globalElemIdx = elementMapper_().map(*elemIt);
            if (isRed)
                elementColor_[globalElemIdx] = Red;
            else
                elementColor_[globalElemIdx] = Green;
        }

        // Mark yellow degrees of freedom (as orange for the mean time)
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Red)
                // non-red elements do not tint degrees of freedom yellow!
                continue;

            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                // if a degree of freedom is already red, don't
                // recolor it to yellow!
                if (dofColor_[globalIdx] != Red) {
                    dofColor_[globalIdx] = Orange;
                }
            }
        }

        // at this point, we communicate the yellow degrees of freedom
        // to the neighboring processes because a neigbor process may
        // not see the red degree of freedom for yellow border degrees
        // of freedom
        auto minHandle =
            GridCommHandleFactory::template minHandle<EntityColor>(dofColor_, dofMapper_());
        gridView_().communicate(*minHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // Mark yellow elements
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Red)
                // element is already red
                continue;

            // check whether the element features a yellow
            // (resp. orange at this point) degree of freedom
            bool isYellow = false;
            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                if (dofColor_[globalIdx] == Orange) {
                    isYellow = true;
                    break;
                }
            }

            if (isYellow)
                elementColor_[elemIdx] = Yellow;
        }

        // Demote orange degrees of freedom to yellow ones if it has
        // at least one green element as a neighbor.
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Green)
                // yellow and red elements do not make orange degrees
                // of freedom yellow!
                continue;

            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                // if a degree of freedom is orange, recolor it to yellow!
                if (dofColor_[globalIdx] == Orange)
                    dofColor_[globalIdx] = Yellow;
            }
        }

        // demote the border orange degrees of freedom
        const auto maxHandle =
            GridCommHandleFactory::template maxHandle<EntityColor>(dofColor_, dofMapper_());
        gridView_().communicate(*maxHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // promote the remaining orange degrees of freedom to red
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            // if a degree of freedom is green or yellow don't do anything!
            if (dofColor_[i] == Green || dofColor_[i] == Yellow)
                continue;

            // make sure the degree of freedom is red (this is a no-op
            // degrees of freedom which are already red!)
            dofColor_[i] = Red;

            // set the error of this degree of freedom to 0 because
            // the system will be relinearized at this dof
            dofError_[i] = 0.0;
        }
    }

    /*!
     * \brief Returns the relinearization color of a degree of freedom
     *
     * \copydetails Doxygen::elementParam
     * \copydetails Doxygen::dofIdxParam
     */
    int dofColor(const ElementContext &elemCtx, int dofIdx) const
    {
        if (!enablePartialRelinearization_())
            return Red;

        int globalIdx = elemCtx.globalSpaceIdx(dofIdx, /*timeIdx=*/0);
        return dofColor_[globalIdx];
    }

    /*!
     * \brief Returns the relinearization color of a degree of freedom
     *
     * \param globalDofIdx The global index of the degree of freedom.
     */
    int dofColor(int globalDofIdx) const
    {
        if (!enablePartialRelinearization_())
            return Red;

        return dofColor_[globalDofIdx];
    }

    /*!
     * \brief Returns the relinearization color of an element
     *
     * \copydetails Doxygen::elementParam
     */
    int elementColor(const Element &element) const
    {
        if (!enablePartialRelinearization_())
            return Red;

        return elementColor_[elementMapper_().map(element)];
    }

    /*!
     * \brief Returns the relinearization color of an element
     *
     * \param globalElementIdx The global index of the element.
     */
    int elementColor(int globalElementIdx) const
    {
        if (!enablePartialRelinearization_())
            return Red;

        return elementColor_[globalElementIdx];
    }

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
        // initialize the BCRS matrix
        createMatrix_();

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseLinearization_ = false;

        int numGridDof = model_().numGridDof();
        int numAllDof =  model_().numTotalDof();
        int numElems = gridView_().size(/*codim=*/0);

        // create the per-thread objects
        elementCtx_.resize(ThreadManager::maxThreads());
        for (int threadId = 0; threadId != ThreadManager::maxThreads(); ++ threadId)
            elementCtx_[threadId] = new ElementContext(simulator_());

        residual_.resize(numAllDof);

        // initialize the storage part of the Jacobian matrix. Since
        // we only need this if Jacobian matrix recycling is enabled,
        // we do not waste space if it is disabled
        if (enableLinearizationRecycling_()) {
            storageJacobian_.resize(numGridDof);
            storageTerm_.resize(numGridDof);
        }

        totalElems_ = gridView_().comm().sum(numElems);

        // initialize data needed for partial relinearization
        if (enablePartialRelinearization_()) {
            dofColor_.resize(numGridDof);
            dofError_.resize(numGridDof);
            elementColor_.resize(numElems);
        }
        relinearizeAll();
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
        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            stencil.update(*eIt);

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
        // do not do anything if we can re-use the current linearization
        if (reuseLinearization_)
            return;

        size_t numGridDof = model_().numGridDof();

        // reset the right hand side.
        residual_ = 0.0;

        if (!enablePartialRelinearization_()) {
            // If partial relinearization of the jacobian is not enabled,
            // we can just reset everything!
            (*matrix_) = 0;

            // reset the parts needed for Jacobian recycling
            if (enableLinearizationRecycling_()) {
                for (unsigned i=0; i < numGridDof; ++ i) {
                    storageJacobian_[i] = 0;
                    storageTerm_[i] = 0;
                }
            }

            return;
        }

        // reset all entries corrosponding to a red or yellow degree
        // of freedom
        for (unsigned rowIdx = 0; rowIdx < numGridDof; ++rowIdx) {
            if (dofColor_[rowIdx] == Green)
                // the equations for this control volume are already below the treshold
                continue;

            // here we have yellow or red degrees of freedom...

            // reset the parts needed for Jacobian recycling
            if (enableLinearizationRecycling_()) {
                storageJacobian_[rowIdx] = 0;
                storageTerm_[rowIdx] = 0;
            }

            // set all matrix entries in the row to 0
            typedef typename JacobianMatrix::ColIterator ColIterator;
            ColIterator colIt = (*matrix_)[rowIdx].begin();
            const ColIterator &colEndIt = (*matrix_)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                (*colIt) = 0.0;
            }
        }
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().simulator().timeStepSize();
        if (reuseLinearization_) {
            int numGridDof = model_().numGridDof();
            for (int dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                // rescale the mass term of the jacobian matrix
                MatrixBlock &J_ii = (*matrix_)[dofIdx][dofIdx];

                J_ii -= storageJacobian_[dofIdx];
                storageJacobian_[dofIdx] *= oldDt_/curDt;
                J_ii += storageJacobian_[dofIdx];

                // use the flux term plus the source term as the new residual (since the
                // delta in the d(storage)/dt is 0 for the first iteration and the
                // residual is approximately 0 in the last iteration, the flux term plus
                // the source term must be equal to the negative change of the storage
                // term of the last iteration of the last time step...)
                residual_[dofIdx] = storageTerm_[dofIdx];
                residual_[dofIdx] *= -1;
            }

            reuseLinearization_ = false;
            oldDt_ = curDt;

            assembleAuxiliaryEquations_();

            problem_().newtonMethod().endIterMsg()
                << ", linear system of equations reused from previous time step";
            return;
        }

        oldDt_ = curDt;
        greenElems_ = 0;

        // relinearize the elements...
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(gridView_());
#if HAVE_OPENMP
#pragma omp parallel
#endif
        {
            ElementIterator elemIt = gridView_().template begin</*codim=*/0>();
            for (threadedElemIt.beginParallel(elemIt);
                 !threadedElemIt.isFinished(elemIt);
                 threadedElemIt.increment(elemIt))
            {
                const Element &elem = *elemIt;

                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

                assembleElement_(elem);
            }
        }

        assembleAuxiliaryEquations_();
    }

    // assemble an element in the interior of the process' grid
    // partition
    void assembleElement_(const Element &elem)
    {
        if (enablePartialRelinearization_()) {
            int globalElemIdx = model_().elementMapper().map(elem);
            if (elementColor_[globalElemIdx] == Green) {
                ++greenElems_;

                assembleGreenElement_(elem);
                return;
            }
        }

        int threadId = ThreadManager::threadId();

        ElementContext *elementCtx = elementCtx_[threadId];
        auto &localJacobian = model_().localJacobian(threadId);

        elementCtx->updateAll(elem);
        localJacobian.assemble(*elementCtx);

        ScopedLock addLock(globalMatrixMutex_);
        for (int primaryDofIdx = 0;
             primaryDofIdx < elementCtx->numPrimaryDof(/*timeIdx=*/0);
             ++ primaryDofIdx)
        {
            int globI = elementCtx->globalSpaceIndex(/*spaceIdx=*/primaryDofIdx,
                                                     /*timeIdx=*/0);

            // update the right hand side
            residual_[globI] += localJacobian.residual(primaryDofIdx);

            if (enableLinearizationRecycling_()) {
                storageTerm_[globI] += localJacobian.residualStorage(primaryDofIdx);
            }

            // only update the jacobian matrix for non-green degrees of freedom
            if (dofColor(globI) != Green) {
                if (enableLinearizationRecycling_())
                    storageJacobian_[globI] += localJacobian.jacobianStorage(primaryDofIdx);

                // update the jacobian matrix
                for (int dofIdx = 0; dofIdx < elementCtx->numDof(/*timeIdx=*/0); ++ dofIdx) {
                    int globJ = elementCtx->globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);
                    (*matrix_)[globI][globJ] += localJacobian.jacobian(primaryDofIdx, dofIdx);
                }
            }
        }
        addLock.unlock();
    }

    // "assemble" a green element. green elements only get the
    // residual updated, but the jacobian is left alone...
    void assembleGreenElement_(const Element &elem)
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
            residual_[globI] += localResidual.residual(dofIdx);
            if (enableLinearizationRecycling_())
                storageTerm_[globI] += localResidual.storageTerm(dofIdx);
        }
        addLock.unlock();
    }

    void assembleAuxiliaryEquations_()
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
    // time step size of last assembly
    Scalar oldDt_;

    // data required for partial relinearization
    std::vector<EntityColor> dofColor_;
    std::vector<Scalar> dofError_;
    std::vector<EntityColor> elementColor_;

    int totalElems_;
    int greenElems_;

    Scalar nextRelinearizationAccuracy_;
    Scalar relinearizationAccuracy_;

    OmpMutex globalMatrixMutex_;
};

} // namespace Ewoms

#endif
