// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <iostream>
#include <vector>
#include <set>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The common code for the assemblers of the global Jacobian
 *        matrix of models using an implicit finite volumes
 *        discretization scheme.
 */
template<class TypeTag>
class FvBaseAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
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
     * \brief The colors of elements and DOFs required for partial
     *        Jacobian reassembly.
     */
    enum EntityColor {
        /*!
         * DOF/element that needs to be reassembled because some
         * relative error is above the tolerance
         */
        Red = 0,

        /*!
         * DOF/element that needs to be reassembled because a
         * neighboring element/DOF is red
         */
        Yellow = 1,

        /*!
         * Yellow DOF has only non-green neighbor elements.
         *
         * This means that its relative error is below the tolerance,
         * but its defect can be linearized without any additional
         * cost. This is just an "internal" color which is not used
         * ouside of the jacobian assembler.
         */
        Orange = 2,

        /*!
         * DOF/element that does not need to be reassembled
         */
        Green = 3
    };

    FvBaseAssembler()
    {
        problemPtr_ = 0;
        elementCtx_ = 0;

        matrix_ = 0;

        // set reassemble accuracy to 0, so that if partial reassembly
        // of the jacobian matrix is disabled, the reassemble accuracy
        // is always smaller than the current relative tolerance
        reassembleAccuracy_ = 0.0;
    }

    ~FvBaseAssembler()
    {
        delete matrix_;
        delete elementCtx_;
    }

    /*!
     * \brief Register all run-time parameters for the Jacobian assembler.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableJacobianRecycling, "Re-use of the jacobian matrix at the first iteration of the next time step");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnablePartialReassemble, "Re-assemble only those degrees of freedom that have changed 'sufficiently' be changed between two Newton iterations");
    }

    /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated. We can not assume that they are
     * fully initialized, though.
     *
     * \copydetails Doxygen::problemParam
     */
    void init(Problem& problem)
    {
        problemPtr_ = &problem;
        delete elementCtx_;
        elementCtx_ = new ElementContext(problem);

        // initialize the BCRS matrix
        createMatrix_();

        // initialize the jacobian matrix and the right hand side
        // vector
        *matrix_ = 0;
        reuseMatrix_ = false;

        int numDof = problem.model().numDof();
        int numElems = gridView_().size(/*codim=*/0);

        residual_.resize(numDof);

        // initialize the storage part of the Jacobian matrix. Since
        // we only need this if Jacobian matrix recycling is enabled,
        // we do not waste space if it is disabled
        if (enableJacobianRecycling_()) {
            storageJacobian_.resize(numDof);
            storageTerm_.resize(numDof);
        }

        totalElems_ = gridView_().comm().sum(numElems);

        // initialize data needed for partial reassembly
        if (enablePartialReassemble_()) {
            dofColor_.resize(numDof);
            dofDelta_.resize(numDof);
            elementColor_.resize(numElems);
        }
        reassembleAll();
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool printReassembleStatistics = enablePartialReassemble_() && !reuseMatrix_;
        int succeeded;
        try {
            assemble_();
            succeeded = 1;
            succeeded = gridView_().comm().min(succeeded);
        }
        catch (Opm::NumericalProblem &e)
        {
            std::cout << "rank " << problem_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = 0;
            succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            OPM_THROW(Opm::NumericalProblem,
                       "A process did not succeed in linearizing the system");
        };

        if (printReassembleStatistics)
        {
            greenElems_ = gridView_().comm().sum(greenElems_);
            reassembleAccuracy_ = gridView_().comm().max(nextReassembleAccuracy_);

            problem_().newtonMethod().endIterMsg()
                << ", reassembled "
                << totalElems_ - greenElems_ << "/" << totalElems_
                << " (" << 100*Scalar(totalElems_ - greenElems_)/totalElems_ << "%) elems @accuracy="
                << reassembleAccuracy_;
        }

        // reset all DOF colors to green
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            dofColor_[i] = Green;
        }
    }

    /*!
     * \brief If Jacobian matrix recycling is enabled, this method
     *        specifies whether the next call to assemble() just
     *        rescales the storage term or does a full reassembly
     *
     * \param yesno If true, only rescale; else do full Jacobian assembly.
     */
    void setMatrixReuseable(bool yesno = true)
    {
        if (enableJacobianRecycling_())
            reuseMatrix_ = yesno;
    }

    /*!
     * \brief If partial Jacobian matrix reassembly is enabled, this
     *        method causes all elements to be reassembled in the next
     *        assemble() call.
     */
    void reassembleAll()
    {
        // do not reuse the current linearization
        reuseMatrix_ = false;

        // do not use partial reassembly for the next iteration
        nextReassembleAccuracy_ = 0.0;
        if (enablePartialReassemble_()) {
            std::fill(dofDelta_.begin(),
                      dofDelta_.end(),
                      0.0);
            std::fill(dofColor_.begin(),
                      dofColor_.end(),
                      Red);
            std::fill(elementColor_.begin(),
                      elementColor_.end(),
                      Red);
        }
    }

    /*!
     * \brief Returns the largest relative error of a "green" DOF
     *        for the most recent call of the assemble() method.
     *
     * This only has an effect if partial Jacobian reassembly is
     * enabled. If it is disabled, then this method always returns 0.
     *
     * This returns the _actual_ relative computed seen by
     * computeColors(), not the tolerance which it was given.
     */
    Scalar reassembleAccuracy() const
    { return reassembleAccuracy_; }

    /*!
     * \brief Update the distance where the non-linear system was
     *        originally insistently linearized and the point where it
     *        will be linerized the next time.
     *
     * This only has an effect if partial reassemble is enabled.
     *
     * \param u The iterative solution of the last iteration of the Newton method
     * \param uDelta The negative difference between the last and the next iterative solution.
     */
    void updateDiscrepancy(const SolutionVector &u,
                           const GlobalEqVector &uDelta)
    {
        if (!enablePartialReassemble_())
            return;

        // update the vector with the distances of the current
        // evaluation point used for linearization from the original
        // evaluation point
        for (unsigned i = 0; i < dofDelta_.size(); ++i) {
            PrimaryVariables uCurrent(u[i]);
            PrimaryVariables uNext(uCurrent);
            uNext -= uDelta[i];

            // we need to add the distance the solution was moved for
            // this DOF
            Scalar dist = model_().relativeDofError(i,
                                                    uCurrent,
                                                    uNext);
            dofDelta_[i] += std::abs(dist);
        }

    }

    /*!
     * \brief Force to reassemble a given degree of freedom the next
     *        time the assemble() method is called.
     *
     * \param globalDofIdx The global index of the DOF which ought
     *                     to be red.
     */
    void markDofRed(int globalDofIdx)
    {
        if (!enablePartialReassemble_())
            return;

        dofColor_[globalDofIdx] = Red;
    }

    /*!
     * \brief Determine the colors of dofs and elements for partial
     *        reassembly given a relative tolerance.
     *
     * The following approach is used:
     *
     * - Set all DOFs and elements to 'green'
     * - Mark all DOFs as 'red' which exhibit an relative error above
     *   the tolerance
     * - Mark all elements which feature a 'red' DOF as 'red'
     * - Mark all DOFs which are not 'red' and are part of a
     *   'red' element as 'yellow'
     * - Mark all elements which are not 'red' and contain a
     *   'yellow' DOF as 'yellow'
     *
     * \param relTol The relative error below which a DOF won't be
     *               reassembled. Note that this specifies the
     *               worst-case relative error between the last
     *               linearization point and the current solution and
     *               _not_ the delta vector of the Newton iteration!
     */
    void computeColors(Scalar relTol)
    {
        if (!enablePartialReassemble_())
            return;

        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();

        // mark the red DOFs and update the tolerance of the
        // linearization which actually will get achieved
        nextReassembleAccuracy_ = 0;
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            if (dofDelta_[i] > relTol)
                // mark dof as red if discrepancy is larger than
                // the relative tolerance
                dofColor_[i] = Red;
            else
                nextReassembleAccuracy_ =
                    std::max(nextReassembleAccuracy_, dofDelta_[i]);
        };

        Stencil stencil(gridView_());
        // Mark all red elements
        for (; elemIt != elemEndIt; ++elemIt) {
            stencil.update(*elemIt);

            // find out whether the current element features a red
            // dof
            bool isRed = false;
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                if (dofColor_[globalIdx] == Red) {
                    isRed = true;
                    break;
                }
            };

            // if yes, the element color is also red, else it is not
            // red, i.e. green for the mean time
            int globalElemIdx = elementMapper_().map(*elemIt);
            if (isRed)
                elementColor_[globalElemIdx] = Red;
            else
                elementColor_[globalElemIdx] = Green;
        }

        // Mark yellow DOFs (as orange for the mean time)
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Red)
                continue; // non-red elements do not tint DOFs
                          // yellow!

            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                // if a DOF is already red, don't recolor it to
                // yellow!
                if (dofColor_[globalIdx] != Red) {
                    dofColor_[globalIdx] = Orange;
                }
            };
        }

        // at this point we communicate the yellow DOFs to the
        // neighboring processes because a neigbor process may not see
        // the red DOF for yellow border DOFs
        auto minHandle = GridCommHandleFactory::template minHandle<EntityColor>(dofColor_, dofMapper_());
        gridView_().communicate(*minHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // Mark yellow elements
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] == Red) {
                continue; // element is already red!
            }

            // check whether the element features a yellow
            // (resp. orange at this point) DOF
            bool isYellow = false;
            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                if (dofColor_[globalIdx] == Orange) {
                    isYellow = true;
                    break;
                }
            };

            if (isYellow)
                elementColor_[elemIdx] = Yellow;
        }

        // Demote orange DOFs to yellow ones if it has at least
        // one green element as a neighbor.
        elemIt = gridView_().template begin<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = this->elementMapper_().map(*elemIt);
            if (elementColor_[elemIdx] != Green)
                continue; // yellow and red elements do not make
                          // orange DOFs yellow!

            stencil.update(*elemIt);
            int numDof = stencil.numDof();
            for (int dofIdx=0; dofIdx < numDof; ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                // if a DOF is orange, recolor it to yellow!
                if (dofColor_[globalIdx] == Orange)
                    dofColor_[globalIdx] = Yellow;
            };
        }

        // demote the border orange DOFs
        const auto maxHandle = GridCommHandleFactory::template maxHandle<EntityColor>(dofColor_, dofMapper_());
        gridView_().communicate(*maxHandle,
                                Dune::InteriorBorder_InteriorBorder_Interface,
                                Dune::ForwardCommunication);

        // promote the remaining orange DOFs to red
        for (unsigned i = 0; i < dofColor_.size(); ++i) {
            // if a DOF is green or yellow don't do anything!
            if (dofColor_[i] == Green || dofColor_[i] == Yellow)
                continue;

            // make sure the DOF is red (this is a no-op DOFs
            // which are already red!)
            dofColor_[i] = Red;

            // set the error of this DOF to 0 because the system
            // will be consistently linearized at this dof
            dofDelta_[i] = 0.0;
        };
    }

    /*!
     * \brief Returns the reassemble color of a DOF
     *
     * \copydetails Doxygen::elementParam
     * \copydetails Doxygen::dofIdxParam
     */
    int dofColor(const ElementContext &elemCtx, int dofIdx) const
    {
        if (!enablePartialReassemble_())
            return Red;

        int globalIdx = elemCtx.globalSpaceIdx(dofIdx, /*timeIdx=*/0);
        return dofColor_[globalIdx];
    }

    /*!
     * \brief Returns the reassemble color of a DOF
     *
     * \param globalDofIdx The global index of the DOF.
     */
    int dofColor(int globalDofIdx) const
    {
        if (!enablePartialReassemble_())
            return Red;

        return dofColor_[globalDofIdx];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \copydetails Doxygen::elementParam
     */
    int elementColor(const Element &element) const
    {
        if (!enablePartialReassemble_())
            return Red;

        return elementColor_[elementMapper_().map(element)];
    }

    /*!
     * \brief Returns the Jacobian reassemble color of an element
     *
     * \param globalElementIdx The global index of the element.
     */
    int elementColor(int globalElementIdx) const
    {
        if (!enablePartialReassemble_())
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
    static bool enableJacobianRecycling_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableJacobianRecycling); }
    static bool enablePartialReassemble_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnablePartialReassemble); }

    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return problem_().model(); }
    Model &model_()
    { return problem_().model(); }
    const GridView &gridView_() const
    { return model_().gridView(); }
    const ElementMapper &elementMapper_() const
    { return model_().elementMapper(); }
    const DofMapper &dofMapper_() const
    { return model_().dofMapper(); }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        int numDof = model_().numDof();

        // allocate raw matrix
        matrix_ = new Matrix(numDof, numDof, Matrix::random);

        Stencil stencil(gridView_());

        // find out the global indices of the neighboring DOFs of
        // each primary DOF
        typedef std::set<int> NeighborSet;
        std::vector<NeighborSet> neighbors(numDof);
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
        };

        // allocate space for the rows of the matrix
        for (int i = 0; i < numDof; ++i) {
            matrix_->setrowsize(i, neighbors[i].size());
        }
        matrix_->endrowsizes();

        // fill the rows with indices. each DOF talks to all of its
        // neighbors. (it also talks to itself since DOFs are
        // sometimes quite egocentric.)
        for (int i = 0; i < numDof; ++i) {
            typename NeighborSet::iterator nIt = neighbors[i].begin();
            typename NeighborSet::iterator nEndIt = neighbors[i].end();
            for (; nIt != nEndIt; ++nIt) {
                matrix_->addindex(i, *nIt);
            }
        }
        matrix_->endindices();
    }

    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // do not do anything if we can re-use the current linearization
        if (reuseMatrix_)
            return;

        // reset the right hand side.
        residual_ = 0.0;

        if (!enablePartialReassemble_()) {
            // If partial reassembly of the jacobian is not enabled,
            // we can just reset everything!
            (*matrix_) = 0;

            // reset the parts needed for Jacobian recycling
            if (enableJacobianRecycling_()) {
                int numDof = matrix_->N();
                for (int i=0; i < numDof; ++ i) {
                    storageJacobian_[i] = 0;
                    storageTerm_[i] = 0;
                };
            }

            return;
        }

        // reset all entries corrosponding to a red or yellow DOF
        for (unsigned rowIdx = 0; rowIdx < matrix_->N(); ++rowIdx) {
            if (dofColor_[rowIdx] == Green)
                continue; // the equations for this control volume are
                          // already below the treshold

            // here we have yellow or red DOFs...

            // reset the parts needed for Jacobian recycling
            if (enableJacobianRecycling_()) {
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
        };
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // if we can "recycle" the current linearization, we do it
        // here and be done with it...
        Scalar curDt = problem_().timeManager().timeStepSize();
        if (reuseMatrix_) {
            int numDof = storageJacobian_.size();
            for (int i = 0; i < numDof; ++i) {
                // rescale the mass term of the jacobian matrix
                MatrixBlock &J_i_i = (*matrix_)[i][i];

                J_i_i -= storageJacobian_[i];
                storageJacobian_[i] *= oldDt_/curDt;
                J_i_i += storageJacobian_[i];

                // use the flux term plus the source term as the new
                // residual (since the delta in the d(storage)/dt is 0
                // for the first iteration and the residual is
                // approximately 0 in the last iteration, the flux
                // term plus the source term must be equal to the
                // negative change of the storage term of the last
                // iteration of the last time step...)
                residual_[i] = storageTerm_[i];
                residual_[i] *= -1;
            };

            reuseMatrix_ = false;
            oldDt_ = curDt;
            return;
        }

        oldDt_ = curDt;
        greenElems_ = 0;

        // reassemble the elements...
        ElementIterator elemIt = gridView_().template begin<0>();
        ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element &elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity  &&
                elem.partitionType() != Dune::BorderEntity) {
                assembleGhostElement_(elem);
            }
            else
                assembleElement_(elem);
        };
    }


    // assemble an element in the interior of the process' grid
    // partition
    void assembleElement_(const Element &elem)
    {
        if (enablePartialReassemble_()) {
            int globalElemIdx = model_().elementMapper().map(elem);
            if (elementColor_[globalElemIdx] == Green) {
                ++greenElems_;

                assembleGreenElement_(elem);
                return;
            }
        }

        elementCtx_->updateAll(elem);
        model_().localJacobian().assemble(*elementCtx_);

        for (int primaryDofIdx = 0; primaryDofIdx < elementCtx_->numPrimaryDof(/*timeIdx=*/0); ++ primaryDofIdx) {
            int globI = elementCtx_->globalSpaceIndex(/*spaceIdx=*/primaryDofIdx, /*timeIdx=*/0);

            // update the right hand side
            residual_[globI] += model_().localJacobian().residual(primaryDofIdx);

            if (enableJacobianRecycling_()) {
                storageTerm_[globI] +=
                    model_().localJacobian().residualStorage(primaryDofIdx);
            }

            // only update the jacobian matrix for non-green degrees of freedom
            if (dofColor(globI) != Green) {
                if (enableJacobianRecycling_())
                    storageJacobian_[globI] +=
                        model_().localJacobian().jacobianStorage(primaryDofIdx);

                // update the jacobian matrix
                for (int dofIdx = 0; dofIdx < elementCtx_->numDof(/*timeIdx=*/0); ++ dofIdx) {
                    int globJ = elementCtx_->globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);
                    (*matrix_)[globI][globJ] +=
                        model_().localJacobian().jacobian(primaryDofIdx, dofIdx);
                }
            }
        }
    }

    // "assemble" a green element. green elements only get the
    // residual updated, but the jacobian is left alone...
    void assembleGreenElement_(const Element &elem)
    {
        elementCtx_->updateAll(elem);
        model_().localResidual().eval(*elementCtx_);

        for (int dofIdx=0; dofIdx < elementCtx_->numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
            int globI = elementCtx_->globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            // update the right hand side
            residual_[globI] += model_().localResidual().residual(dofIdx);
            if (enableJacobianRecycling_())
                storageTerm_[globI] += model_().localResidual().storageTerm(dofIdx);
        }
    }

    // "assemble" a ghost element
    void assembleGhostElement_(const Element &elem)
    {
        const auto &stencil = elementCtx_->stencil(/*timeIdx=*/0);
        elementCtx_->updateStencilTopology(elem);

        for (int dofIdx = 0; dofIdx < stencil.numPrimaryDof(); ++dofIdx) {
            int partitionType = stencil.partitionType(dofIdx);

            if (partitionType == Dune::InteriorEntity ||
                partitionType == Dune::BorderEntity)
            {
                // do not change the non-ghost dofs
                continue;
            }

            // set main diagonal entries for the vertex
            int globalIdx = stencil.globalSpaceIndex(dofIdx);
            typedef typename Matrix::block_type BlockType;
            BlockType &J = (*matrix_)[globalIdx][globalIdx];
            for (int j = 0; j < BlockType::rows; ++j)
                J[j][j] = 1.0;

            // set residual for the dof
            residual_[globalIdx] = 0;
        }
     }

    Problem *problemPtr_;
    ElementContext *elementCtx_;

    // the jacobian matrix
    Matrix *matrix_;
    // the right-hand side
    GlobalEqVector residual_;

    // attributes required for jacobian matrix recycling
    bool reuseMatrix_;
    // The storage part of the local Jacobian
    std::vector<MatrixBlock> storageJacobian_;
    std::vector<VectorBlock> storageTerm_;
    // time step size of last assembly
    Scalar oldDt_;

    // attributes required for partial jacobian reassembly
    std::vector<EntityColor> dofColor_;
    std::vector<Scalar> dofDelta_;
    std::vector<EntityColor> elementColor_;

    int totalElems_;
    int greenElems_;

    Scalar nextReassembleAccuracy_;
    Scalar reassembleAccuracy_;
};

} // namespace Ewoms

#endif
