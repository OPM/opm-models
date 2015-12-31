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
    FvBaseLinearizer()
    {
        simulatorPtr_ = 0;

        matrix_ = 0;
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
    { }

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

        int succeeded;
        try {
            linearize_();
            succeeded = 1;
            succeeded = gridView_().comm().min(succeeded);
        }
        catch (Opm::NumericalProblem &e)
        {
            std::cout << "rank " << simulator_().gridView().comm().rank()
                      << " caught an exception while linearizing:" << e.what()
                      << "\n"  << std::flush;
            succeeded = 0;
            succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            OPM_THROW(Opm::NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }
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
        // initialize the BCRS matrix for the Jacobian
        createMatrix_();

        // initialize the jacobian matrix and the right hand side vector
        *matrix_ = 0;
        residual_.resize(model_().numTotalDof());
        residual_ = 0;

        // create the per-thread context objects
        elementCtx_.resize(ThreadManager::maxThreads());
        for (int threadId = 0; threadId != ThreadManager::maxThreads(); ++ threadId)
            elementCtx_[threadId] = new ElementContext(simulator_());
    }

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        unsigned numAllDof =  model_().numTotalDof();

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

            for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
                unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

                for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                    unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                    neighbors[myIdx].insert(neighborIdx);
                }
            }
        }

        // add the additional neighbors and degrees of freedom caused by the auxiliary
        // equations
        const auto& model = model_();
        unsigned numAuxMod = model.numAuxiliaryModules();
        for (unsigned auxModIdx = 0; auxModIdx < numAuxMod; ++auxModIdx)
            model.auxiliaryModule(auxModIdx)->addNeighbors(neighbors);

        // allocate space for the rows of the matrix
        for (unsigned dofIdx = 0; dofIdx < numAllDof; ++ dofIdx)
            matrix_->setrowsize(dofIdx, neighbors[dofIdx].size());
        matrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (unsigned dofIdx = 0; dofIdx < numAllDof; ++ dofIdx) {
            typename NeighborSet::iterator nIt = neighbors[dofIdx].begin();
            typename NeighborSet::iterator nEndIt = neighbors[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
                matrix_->addindex(dofIdx, *nIt);
        }
        matrix_->endindices();
    }

    // reset the global linear system of equations.
    void resetSystem_()
    {
        residual_ = 0.0;
        (*matrix_) = 0;
    }

    // linearize the whole system
    void linearize_()
    {
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
        int threadId = ThreadManager::threadId();

        ElementContext *elementCtx = elementCtx_[threadId];
        auto &localLinearizer = model_().localLinearizer(threadId);

        // the actual work of linearization is done by the local linearizer class
        elementCtx->updateAll(elem);
        localLinearizer.linearize(*elementCtx);

        // update the right hand side and the Jacobian matrix
        ScopedLock addLock(globalMatrixMutex_);
        unsigned numPrimaryDof = elementCtx->numPrimaryDof(/*timeIdx=*/0);
        for (unsigned primaryDofIdx = 0; primaryDofIdx < numPrimaryDof; ++ primaryDofIdx) {
            int globI = elementCtx->globalSpaceIndex(/*spaceIdx=*/primaryDofIdx, /*timeIdx=*/0);

            // update the right hand side
            residual_[globI] += localLinearizer.residual(primaryDofIdx);

            // update the global Jacobian matrix
            for (unsigned dofIdx = 0; dofIdx < elementCtx->numDof(/*timeIdx=*/0); ++ dofIdx) {
                int globJ = elementCtx->globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);

                (*matrix_)[globJ][globI] += localLinearizer.jacobian(dofIdx, primaryDofIdx);
            }
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

    OmpMutex globalMatrixMutex_;
};

} // namespace Ewoms

#endif
