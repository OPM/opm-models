/*
  Copyright (C) 2008-2013 by Andreas Lauser
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
 * \copydoc Ewoms::FvBaseDiscretization
 */
#ifndef EWOMS_FV_BASE_DISCRETIZATION_HH
#define EWOMS_FV_BASE_DISCRETIZATION_HH

#include "fvbaseproperties.hh"
#include "fvbaseassembler.hh"
#include "fvbaselocaljacobian.hh"
#include "fvbaselocalresidual.hh"
#include "fvbaseelementcontext.hh"
#include "fvbaseboundarycontext.hh"
#include "fvbaseconstraintscontext.hh"
#include "fvbaseconstraints.hh"
#include "fvbasediscretization.hh"
#include "fvbasegradientcalculator.hh"
#include "fvbasenewtonmethod.hh"
#include "fvbaseprimaryvariables.hh"
#include "fvbaseintensivequantities.hh"
#include "fvbaseextensivequantities.hh"

#include <ewoms/parallel/gridcommhandles.hh>
#include <ewoms/linear/nullborderlistmanager.hh>
#include <ewoms/common/simulator.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>

#include <limits>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace Ewoms {
template<class TypeTag>
class FvBaseDiscretization;
}

namespace Opm {
namespace Properties {
//! Set the default type for the time manager
SET_TYPE_PROP(FvBaseDiscretization, Simulator, Ewoms::Simulator<TypeTag>);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(FvBaseDiscretization, VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(FvBaseDiscretization, ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! marks the border indices (required for the algebraic overlap stuff)
SET_PROP(FvBaseDiscretization, BorderListCreator)
{
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Ewoms::Linear::NullBorderListCreator<GridView, DofMapper> type;
};

SET_TYPE_PROP(FvBaseDiscretization, DiscLocalResidual, Ewoms::FvBaseLocalResidual<TypeTag>);

SET_TYPE_PROP(FvBaseDiscretization, DiscIntensiveQuantities, Ewoms::FvBaseIntensiveQuantities<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, DiscExtensiveQuantities, Ewoms::FvBaseExtensiveQuantities<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(FvBaseDiscretization, GradientCalculator, Ewoms::FvBaseGradientCalculator<TypeTag>);

SET_TYPE_PROP(FvBaseDiscretization, DiscLocalJacobian, Ewoms::FvBaseLocalJacobian<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, LocalJacobian, typename GET_PROP_TYPE(TypeTag, DiscLocalJacobian));


//! Set the type of a global jacobian matrix from the solution types
SET_PROP(FvBaseDiscretization, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

//! The maximum allowed number of timestep divisions for the
//! Newton solver
SET_INT_PROP(FvBaseDiscretization, MaxTimeStepDivisions, 10);

/*!
 * \brief A vector of quanties, each for one equation.
 */
SET_TYPE_PROP(FvBaseDiscretization, EqVector,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief A vector for mass/energy rates.
 *
 * E.g. Neumann fluxes or source terms
 */
SET_TYPE_PROP(FvBaseDiscretization, RateVector,
              typename GET_PROP_TYPE(TypeTag, EqVector));

/*!
 * \brief Type of object for specifying boundary conditions.
 */
SET_TYPE_PROP(FvBaseDiscretization, BoundaryRateVector,
              typename GET_PROP_TYPE(TypeTag, RateVector));

/*!
 * \brief The class which represents constraints.
 */
SET_TYPE_PROP(FvBaseDiscretization, Constraints, Ewoms::FvBaseConstraints<TypeTag>);

/*!
 * \brief The type for storing a residual for an element.
 */
SET_TYPE_PROP(FvBaseDiscretization, ElementEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief The type for storing a residual for the whole grid.
 */
SET_TYPE_PROP(FvBaseDiscretization, GlobalEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief An object representing a local set of primary variables.
 */
SET_TYPE_PROP(FvBaseDiscretization, PrimaryVariables, Ewoms::FvBasePrimaryVariables<TypeTag>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(FvBaseDiscretization, SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The class representing intensive quantities.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(FvBaseDiscretization, IntensiveQuantities, Ewoms::FvBaseIntensiveQuantities<TypeTag>);

/*!
 * \brief The element context
 */
SET_TYPE_PROP(FvBaseDiscretization, ElementContext, Ewoms::FvBaseElementContext<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, BoundaryContext, Ewoms::FvBaseBoundaryContext<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, ConstraintsContext, Ewoms::FvBaseConstraintsContext<TypeTag>);

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(FvBaseDiscretization, JacobianAssembler, Ewoms::FvBaseAssembler<TypeTag>);

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(FvBaseDiscretization, MaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(FvBaseDiscretization, MaxTimeStepSize, 1e100);
#endif
//! By default, accept any time step larger than zero
SET_SCALAR_PROP(FvBaseDiscretization, MinTimeStepSize, 0.0);

//! The base epsilon value for finite difference calculations
SET_SCALAR_PROP(FvBaseDiscretization, BaseEpsilon, 0.9123e-10);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(FvBaseDiscretization, NumericDifferenceMethod, +1);

//! Enable the VTK output by default
SET_BOOL_PROP(FvBaseDiscretization, EnableVtkOutput, true);

// disable linearization recycling by default
SET_BOOL_PROP(FvBaseDiscretization, EnableLinearizationRecycling, false);

// disable partial relinearization by default
SET_BOOL_PROP(FvBaseDiscretization, EnablePartialRelinearization, false);

// disable constraints by default
SET_BOOL_PROP(FvBaseDiscretization, EnableConstraints, false);

// by default, disable the intensive quantity cache. If the intensive quantities are
// relatively cheap to calculate, the cache basically does not yield any performance
// impact because of the intensive quantity cache will cause additional pressure on the
// CPU caches...
SET_BOOL_PROP(FvBaseDiscretization, EnableIntensiveQuantityCache, false);

// do not use thermodynamic hints by default. If you enable this, make sure to also
// enable the intensive quantity cache above to avoid getting an exception...
SET_BOOL_PROP(FvBaseDiscretization, EnableThermodynamicHints, false);

// if the deflection of the newton method is large, we do not need to solve the linear
// approximation accurately. Assuming that the value for the current solution is quite
// close to the final value, a reduction of 3 orders of magnitude in the defect should be
// sufficient...
SET_SCALAR_PROP(FvBaseDiscretization, LinearSolverTolerance, 1e-3);

//! Set the history size of the time discretization to 2 (for implicit euler)
SET_INT_PROP(FvBaseDiscretization, TimeDiscHistorySize, 2);

//! Most models don't need the gradients at the center of the SCVs, so
//! we disable them by default.
SET_BOOL_PROP(FvBaseDiscretization, RequireScvCenterGradients, false);
}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The base class for the finite volume discretization schemes.
 */
template<class TypeTag>
class FvBaseDiscretization
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryContext) BoundaryContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, GradientCalculator) GradientCalculator;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, DiscBaseOutputModule) DiscBaseOutputModule;
    typedef typename GET_PROP_TYPE(TypeTag, GridCommHandleFactory) GridCommHandleFactory;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),
        dim = GridView::dimension
    };

    typedef std::vector<IntensiveQuantities> IntensiveQuantitiesVector;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying a discretization object is not a good idea
    FvBaseDiscretization(const FvBaseDiscretization &);

public:
    // this constructor required to be explicitly specified because
    // we've defined a constructor above which deletes all implicitly
    // generated constructors in C++.
    FvBaseDiscretization(Simulator &simulator)
        : simulator_(simulator)
        , gridView_(simulator.gridView())
        , newtonMethod_(simulator)
        , jacAsm_(new JacobianAssembler())
    {
        asImp_().updateBoundary_();

        int nDofs = asImp_().numDof();
        for (int timeIdx = 0; timeIdx < historySize; ++timeIdx) {
            solution_[timeIdx].resize(nDofs);

            if (storeIntensiveQuantities_()) {
                intensiveQuantityCache_[timeIdx].resize(nDofs);
                intensiveQuantityCacheUpToDate_[timeIdx].resize(nDofs, /*value=*/false);
            }
        }

        asImp_().registerOutputModules_();
    }

    ~FvBaseDiscretization()
    {
        // delete all output modules
        auto modIt = outputModules_.begin();
        const auto &modEndIt = outputModules_.end();
        for (; modIt != modEndIt; ++modIt)
            delete *modIt;

        delete jacAsm_;
    }

    /*!
     * \brief Register all run-time parameters for the model.
     */
    static void registerParameters()
    {
        JacobianAssembler::registerParameters();
        LocalJacobian::registerParameters();
        LocalResidual::registerParameters();
        GradientCalculator::registerParameters();
        IntensiveQuantities::registerParameters();
        ExtensiveQuantities::registerParameters();
        NewtonMethod::registerParameters();

        // register runtime parameters of the output modules
        Ewoms::VtkPrimaryVarsModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableVtkOutput, "Global switch for turing on writing VTK files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableThermodynamicHints, "Enable thermodynamic hints");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableIntensiveQuantityCache, "Turn on caching of intensive quantities");
    }

    /*!
     * \brief Apply the initial conditions to the model.
     */
    void finishInit()
    {
        // initialize the volume of the finite volumes to zero
        int nDofs = asImp_().numDof();
        dofTotalVolume_.resize(nDofs);
        std::fill(dofTotalVolume_.begin(), dofTotalVolume_.end(), 0.0);

        ElementContext elemCtx(simulator_);

        // iterate through the grid and evaluate the initial condition
        ElementIterator it = gridView_.template begin</*codim=*/0>();
        const ElementIterator &eendit = gridView_.template end</*codim=*/0>();
        for (; it != eendit; ++it) {
            // ignore everything which is not in the interior if the
            // current process' piece of the grid
            if (it->partitionType() != Dune::InteriorEntity)
                continue;

            // deal with the current element
            elemCtx.updateStencil(*it);

            // loop over all element vertices, i.e. sub control volumes
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); dofIdx++)
            {
                // map the local degree of freedom index to the global one
                unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                dofTotalVolume_[globalIdx] +=
                    elemCtx.stencil(/*timeIdx=*/0).subControlVolume(dofIdx).volume();
            }
        }

        const auto sumHandle =
            GridCommHandleFactory::template sumHandle<double>(dofTotalVolume_,
                                                              asImp_().dofMapper());
        gridView_.communicate(*sumHandle,
                              Dune::InteriorBorder_InteriorBorder_Interface,
                              Dune::ForwardCommunication);

        localJacobian_.init(simulator_);
        jacAsm_->init(simulator_);

        if (storeIntensiveQuantities_()) {
            // invalidate all cached intensive quantities
            for (int timeIdx = 0; timeIdx < historySize; ++ timeIdx) {
                std::fill(intensiveQuantityCacheUpToDate_[timeIdx].begin(),
                          intensiveQuantityCacheUpToDate_[timeIdx].end(),
                          false);
            }
        }
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom to which the model
     *        applies.
     */
    void applyInitialSolution()
    {
        // first set the whole domain to zero
        SolutionVector &uCur = asImp_().solution(/*timeIdx=*/0);
        uCur = Scalar(0.0);

        ElementContext elemCtx(simulator_);

        // iterate through the grid and evaluate the initial condition
        ElementIterator it = gridView_.template begin</*codim=*/0>();
        const ElementIterator &eendit = gridView_.template end</*codim=*/0>();
        for (; it != eendit; ++it) {
            // ignore everything which is not in the interior if the
            // current process' piece of the grid
            if (it->partitionType() != Dune::InteriorEntity)
                continue;

            // deal with the current element
            elemCtx.updateStencil(*it);

            // loop over all element vertices, i.e. sub control volumes
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); dofIdx++)
            {
                // map the local degree of freedom index to the global one
                unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                // let the problem do the dirty work of nailing down
                // the initial solution.
                simulator_.problem().initial(uCur[globalIdx], elemCtx, dofIdx, /*timeIdx=*/0);
                uCur[globalIdx].checkDefined();
            }
        }

        // synchronize the ghost DOFs (if necessary)
        asImp_().syncOverlap();

        // also set the solution of the "previous" time steps to the
        // initial solution.
        for (int timeIdx = 1; timeIdx < historySize; ++timeIdx)
            solution_[timeIdx] = solution_[/*timeIdx=*/0];

    }

    /*!
     * \brief Returns the newton method object
     */
    NewtonMethod &newtonMethod()
    { return newtonMethod_; }

    /*!
     * \copydoc newtonMethod()
     */
    const NewtonMethod &newtonMethod() const
    { return newtonMethod_; }

    /*!
     * \brief Return the thermodynamic hint for a entity on the grid at given time.
     *
     * The hint is defined as a IntensiveQuantities object which is supposed to be
     * "close" to the IntensiveQuantities of the current solution. It can be used as a
     * good starting point for non-linear solvers when having to solve non-linear
     * relations while updating the intensive quantities. (This may yield a major
     * performance boost depending on how the physical models require.)
     *
     * \attention If no up-to date intensive quantities are available, or if hints have been
     *            disabled, this method will return 0.
     *
     * \param globalIdx The global space index for the entity where a hint is requested.
     * \param timeIdx The index used by the time discretization.
     */
    const IntensiveQuantities *thermodynamicHint(int globalIdx, int timeIdx) const
    {
        if (!enableThermodynamicHints_())
            return 0;

        if (intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])
            return &intensiveQuantityCache_[timeIdx][globalIdx];

        // use the intensive quantities for the first up-to-date time index as hint
        for (int timeIdx2 = 0; timeIdx2 < historySize; ++timeIdx2)
            if (intensiveQuantityCacheUpToDate_[timeIdx2][globalIdx])
                return &intensiveQuantityCache_[timeIdx2][globalIdx];

        // no suitable up-to-date intensive quantities...
        return 0;
    }

    /*!
     * \brief Return the cached intensive quantities for a entity on the
     *        grid at given time.
     *
     * \attention If no up-to date intensive quantities are available,
     *            this method will return 0.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is requested.
     * \param timeIdx The index used by the time discretization.
     */
    const IntensiveQuantities *cachedIntensiveQuantities(int globalIdx, int timeIdx) const
    {
        if (!enableIntensiveQuantitiesCache_() ||
            !intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])
            return 0;

        return &intensiveQuantityCache_[timeIdx][globalIdx];
    }

    /*!
     * \brief Update the intensive quantity cache for a entity on the grid at given time.
     *
     * \param intQuants The IntensiveQuantities object hint for a given degree of freedom.
     * \param globalIdx The global space index for the entity where a
     *                  hint is to be set.
     * \param timeIdx The index used by the time discretization.
     */
    void updateCachedIntensiveQuantities(const IntensiveQuantities &intQuants,
                                         int globalIdx,
                                         int timeIdx) const
    {
        if (!storeIntensiveQuantities_())
            return;

        intensiveQuantityCache_[timeIdx][globalIdx] = intQuants;
        intensiveQuantityCacheUpToDate_[timeIdx][globalIdx] = true;
    }

    /*!
     * \brief Invalidate the cache for a given intensive quantities object.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is to be set.
     * \param timeIdx The index used by the time discretization.
     */
    void invalidateIntensiveQuantitiesCacheEntry(int globalIdx,
                                                 int timeIdx) const
    {
        if (!storeIntensiveQuantities_())
            return;

        intensiveQuantityCacheUpToDate_[timeIdx][globalIdx] = false;
    }

    /*!
     * \brief Move the intensive quantities for a given time index to the back.
     *
     * This method should only be called by the time discretization.
     *
     * \param numSlots The number of time step slots for which the
     *                 hints should be shifted.
     */
    void shiftIntensiveQuantityCache(int numSlots = 1)
    {
        if (!storeIntensiveQuantities_())
            return;

        for (int timeIdx = 0; timeIdx < historySize - numSlots; ++ timeIdx) {
            intensiveQuantityCache_[timeIdx + numSlots] = intensiveQuantityCache_[timeIdx];
            intensiveQuantityCacheUpToDate_[timeIdx + numSlots] = intensiveQuantityCacheUpToDate_[timeIdx];
        }

        // invalidate the cache for the most recent time index
        std::fill(intensiveQuantityCacheUpToDate_[/*timeIdx=*/0].begin(),
                  intensiveQuantityCacheUpToDate_[/*timeIdx=*/0].end(),
                  false);
    }

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param dest Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(GlobalEqVector &dest,
                          const SolutionVector &u) const
    {
        SolutionVector tmp(asImp_().solution(/*timeIdx=*/0));
        solution_[/*timeIdx=*/0] = u;
        Scalar res = asImp_().globalResidual(dest);
        solution_[/*timeIdx=*/0] = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param dest Stores the result
     */
    Scalar globalResidual(GlobalEqVector &dest) const
    {
        dest = 0;

        LocalBlockVector residual, storageTerm;

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateAll(*elemIt);
            residual.resize(elemCtx.numDof(/*timeIdx=*/0));
            storageTerm.resize(elemCtx.numDof(/*timeIdx=*/0));
            asImp_().localResidual().eval(residual, storageTerm, elemCtx);

            int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
            for (int dofIdx = 0; dofIdx < numPrimaryDof; ++dofIdx) {
                int globalI = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                dest[globalI] += residual[dofIdx];
            }
        };

        // add up the residuals on the process borders
        const auto sumHandle =
            GridCommHandleFactory::template sumHandle<EqVector>(dest, asImp_().dofMapper());
        gridView_.communicate(*sumHandle,
                              Dune::InteriorBorder_InteriorBorder_Interface,
                              Dune::ForwardCommunication);

        // calculate the square norm of the residual. this is not
        // entirely correct, since the residual for the finite volumes
        // which are on the boundary are counted once for every
        // process. As often in life: shit happens (, we don't care)...
        Scalar result2 = dest.two_norm2();
        result2 = asImp_().gridView().comm().sum(result2);

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \copydetails Doxygen::storageParam
     */
    void globalStorage(EqVector &storage, int timeIdx = 0) const
    {
        storage = 0;

        LocalBlockVector elemStorage;

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateIntensiveQuantities(timeIdx);

            int numPrimaryDof = elemCtx.numPrimaryDof(timeIdx);
            elemStorage.resize(numPrimaryDof);

            localResidual().evalStorage(elemStorage, elemCtx, timeIdx);

            for (int dofIdx = 0; dofIdx < numPrimaryDof; ++dofIdx)
                storage += elemStorage[dofIdx];
        };

        storage = gridView_.comm().sum(storage);
    }

    /*!
     * \brief Ensure that the difference between the storage terms of the last and of the
     *        current time step is consistent with the source and boundary terms.
     *
     * This method is purely intented for debugging purposes. If the program is compiled
     * with optimizations enabled, it becomes a no-op.
     */
    void checkConservativeness(Scalar tolerance = 1e-6) const
    {
#ifndef NDEBUG
        EqVector storageBeginTimeStep;
        EqVector storageEndTimeStep;

        Scalar totalBoundaryArea(0.0);
        Scalar totalVolume(0.0);
        EqVector totalRate(0.0);

        // we assume the implicit Euler time discretization for now...
        assert(historySize == 2);

        globalStorage(storageBeginTimeStep, /*timeIdx=*/1);
        globalStorage(storageEndTimeStep, /*timeIdx=*/0);

        // calculate the rate at the boundary and the source rate
        ElementContext elemCtx(simulator_);
        auto eIt = simulator_.gridView().template begin</*codim=*/0>();
        const auto &eEndIt = simulator_.gridView().template end</*codim=*/0>();
        for (; eIt != eEndIt; ++eIt) {
            if (eIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateAll(*eIt);

            // handle the boundary terms
            if (elemCtx.onBoundary()) {
                BoundaryContext boundaryCtx(elemCtx);

                for (int faceIdx = 0; faceIdx < boundaryCtx.numBoundaryFaces(/*timeIdx=*/0); ++faceIdx) {
                    BoundaryRateVector values;
                    simulator_.problem().boundary(values,
                                                         boundaryCtx,
                                                         faceIdx,
                                                         /*timeIdx=*/0);
                    Valgrind::CheckDefined(values);

                    int dofIdx = boundaryCtx.interiorScvIndex(faceIdx, /*timeIdx=*/0);
                    const auto &insideIntQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);

                    Scalar bfArea =
                        boundaryCtx.boundarySegmentArea(faceIdx, /*timeIdx=*/0)
                        * insideIntQuants.extrusionFactor();

                    values *= bfArea;

                    totalBoundaryArea += bfArea;
                    totalRate += values;
                }
            }

            // deal with the source terms
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
                BoundaryRateVector values;
                simulator_.problem().source(values,
                                            elemCtx,
                                            dofIdx,
                                            /*timeIdx=*/0);
                Valgrind::CheckDefined(values);

                const auto &intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                Scalar dofVolume =
                    elemCtx.dofVolume(dofIdx, /*timeIdx=*/0)
                    * intQuants.extrusionFactor();
                values *= dofVolume;
                totalVolume += dofVolume;
                totalRate += values;
            }
        }

        // summarize everything over all processes
        const auto &comm = simulator_.gridView().comm();
        totalRate = comm.sum(totalRate);
        totalBoundaryArea = comm.sum(totalBoundaryArea);
        totalVolume = comm.sum(totalVolume);

        if (comm.rank() == 0) {
            EqVector storageRate = storageBeginTimeStep;
            storageRate -= storageEndTimeStep;
            storageRate /= simulator_.timeStepSize();
            for (int eqIdx = 0; eqIdx < EqVector::dimension; ++eqIdx)
                assert(std::abs(storageRate[eqIdx] - totalRate[eqIdx]) <= tolerance);
        }
#endif // NDEBUG
    }

    /*!
     * \brief Returns the volume \f$\mathrm{[m^3]}\f$ of a given control volume.
     *
     * \param globalIdx The global index of the control volume's
     *                  associated vertex
     */
    Scalar dofTotalVolume(int globalIdx) const
    { return dofTotalVolume_[globalIdx]; }

    /*!
     * \brief Reference to the solution at a given history index as a block vector.
     *
     * \param timeIdx The index of the solution used by the time discretization.
     */
    const SolutionVector &solution(int timeIdx) const
    { return solution_[timeIdx]; }

    /*!
     * \copydoc solution(int) const
     */
    SolutionVector &solution(int timeIdx)
    { return solution_[timeIdx]; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianBaseAssembler() const
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }
    /*!
     * \copydoc localJacobian() const
     */
    LocalJacobian &localJacobian()
    { return localJacobian_; }

    /*!
     * \brief Returns the object to calculate the local residual function.
     */
    const LocalResidual &localResidual() const
    { return asImp_().localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual() const
     */
    LocalResidual &localResidual()
    { return asImp_().localJacobian().localResidual(); }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalDofIdx The global index of the degree of freedom
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        Scalar absPv = std::abs(asImp_().solution(/*timeIdx=*/1)[globalDofIdx][pvIdx]);
        return 1.0/std::max(absPv, 1.0);
    }

    /*!
     * \brief Returns the relative weight of an equation
     *
     * \param globalVertexIdx The global index of the vertex
     * \param eqIdx The index of the equation
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    { return 1.0; }

    /*!
     * \brief Returns the relative error between two vectors of
     *        primary variables.
     *
     * \param vertexIdx The global index of the control volume's
     *                  associated vertex
     * \param pv1 The first vector of primary variables
     * \param pv2 The second vector of primary variables
     */
    Scalar relativeDofError(int vertexIdx,
                            const PrimaryVariables &pv1,
                            const PrimaryVariables &pv2) const
    {
        Scalar result = 0.0;
        for (int j = 0; j < numEq; ++j) {
            Scalar weight = asImp_().primaryVarWeight(vertexIdx, j);
            Scalar eqErr = std::abs((pv1[j] - pv2[j])*weight);
            //Scalar eqErr = std::abs(pv1[j] - pv2[j]);
            //eqErr *= std::max<Scalar>(1.0, std::abs(pv1[j] + pv2[j])/2);

            result = std::max(result, eqErr);
        }
        return result;
    }

    /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver The non-linear solver
     */
    bool update(NewtonMethod &solver)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < asImp_().solution(/*timeIdx=*/0).size(); ++i)
            asImp_().solution(/*timeIdx=*/0)[i].checkDefined();
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.apply();
        if (converged) {
            asImp_().updateSuccessful();
        }
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < asImp_().solution(/*timeIdx=*/0).size(); ++i) {
            Valgrind::CheckDefined(asImp_().solution(/*timeIdx=*/0)[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }

    /*!
     * \brief Syncronize the values of the primary variables on the
     *        degrees of freedom that overlap with the neighboring
     *        processes.
     *
     * By default, this method does nothing...
     */
    void syncOverlap()
    { }

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primary a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    { updateBoundary_(); }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    { }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        intensiveQuantityCache_[/*timeIdx=*/0] = intensiveQuantityCache_[/*timeIdx=*/1];
        intensiveQuantityCacheUpToDate_[/*timeIdx=*/0] = intensiveQuantityCacheUpToDate_[/*timeIdx=*/1];

        solution_[/*timeIdx=*/0] = solution_[/*timeIdx=*/1];
        jacAsm_->relinearizeAll();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        // make the current solution the previous one.
        solution_[/*timeIdx=*/1] = solution_[/*timeIdx=*/0];

        // shift the intensive quantities cache by one position in the
        // history
        asImp_().shiftIntensiveQuantityCache(/*numSlots=*/1);
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        OPM_THROW(std::runtime_error,
                  "Not implemented: The discretization chosen for this problem does not support"
                  " restart files. (serialize() method unimplemented)");
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        OPM_THROW(std::runtime_error,
                  "Not implemented: The discretization chosen for this problem does not support"
                  " restart files. (deserialize() method unimplemented)");
    }

    /*!
     * \brief Write the current solution for a degree of freedom to a
     *        restart file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param dof The Dune entity which's data should be serialized
     */
    template <class DofEntity>
    void serializeEntity(std::ostream &outstream,
                         const DofEntity &dof)
    {
        int dofIdx = asImp_().dofMapper().map(dof);

        // write phase state
        if (!outstream.good()) {
            OPM_THROW(std::runtime_error, "Could not serialize degree of freedom " << dofIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << solution_[/*timeIdx=*/0][dofIdx][eqIdx] << " ";
        }
    }

    /*!
     * \brief Reads the current solution variables for a degree of
     *        freedom from a restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param dof The Dune entity which's data should be deserialized
     */
    template <class DofEntity>
    void deserializeEntity(std::istream &instream,
                           const DofEntity &dof)
    {
        int dofIdx = asImp_().dofMapper().map(dof);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                OPM_THROW(std::runtime_error,
                          "Could not deserialize degree of freedom " << dofIdx);
            instream >> solution_[/*timeIdx=*/0][dofIdx][eqIdx];
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDof() const
    { OPM_THROW(std::logic_error,
                "The discretization class must implement the numDof() method!"); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper &dofMapper() const
    { OPM_THROW(std::logic_error,
                "The discretization class must implement the dofMapper() method!"); }

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return simulator_.problem().vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return simulator_.problem().elementMapper(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler;
        jacAsm_->init(simulator_);
    }

    /*!
     * \brief Return whether a degree of freedom is located on the
     *        domain boundary.
     *
     * \param globalIdx The global space index of the degree of freedom of interest.
     */
    bool onBoundary(int globalIdx) const
    { return onBoundary_[globalIdx]; }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return ""; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     *
     * \param pvIdx The index of the primary variable of interest.
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;
        oss << "primary variable_" << pvIdx;
        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     *
     * \param eqIdx The index of the conservation equation of interest.
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        oss << "equation_" << eqIdx;
        return oss.str();
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of intensive quantities
     *
     * \copydetails Doxygen::ecfvElemCtxParam
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    { }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK writer.
     *
     * \param writer The writer object to which the fields should be added.
     * \param u The solution function
     * \param deltaU The delta of the solution function before and after the Newton update
     */
    void addConvergenceVtkFields(Ewoms::VtkMultiWriter<GridView> &writer,
                                 const SolutionVector &u,
                                 const GlobalEqVector &deltaU) const
    {
        typedef std::vector<double> ScalarBuffer;

        GlobalEqVector globalResid(u.size());
        asImp_().globalResidual(globalResid, u);

        // create the required scalar fields
        unsigned numDof = asImp_().numDof();

        // global defect of the two auxiliary equations
        ScalarBuffer* def[numEq];
        ScalarBuffer* delta[numEq];
        ScalarBuffer* priVars[numEq];
        ScalarBuffer* priVarWeight[numEq];
        ScalarBuffer* relError = writer.allocateManagedScalarBuffer(numDof);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            priVars[pvIdx] = writer.allocateManagedScalarBuffer(numDof);
            priVarWeight[pvIdx] = writer.allocateManagedScalarBuffer(numDof);
            delta[pvIdx] = writer.allocateManagedScalarBuffer(numDof);
            def[pvIdx] = writer.allocateManagedScalarBuffer(numDof);
        }

        for (unsigned globalIdx = 0; globalIdx < numDof; ++ globalIdx)
        {
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                (*priVars[pvIdx])[globalIdx] = u[globalIdx][pvIdx];
                (*priVarWeight[pvIdx])[globalIdx] = asImp_().primaryVarWeight(globalIdx, pvIdx);
                (*delta[pvIdx])[globalIdx] = - deltaU[globalIdx][pvIdx];
                (*def[pvIdx])[globalIdx] = globalResid[globalIdx][pvIdx];
            }

            PrimaryVariables uOld(u[globalIdx]);
            PrimaryVariables uNew(uOld);
            uNew -= deltaU[globalIdx];
            (*relError)[globalIdx] = asImp_().relativeDofError(globalIdx, uOld, uNew);
        }

        DiscBaseOutputModule::attachScalarDofData_(writer, *relError, "relErr");

        for (int i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss.str(""); oss << "priVar_" << asImp_().primaryVarName(i);
            DiscBaseOutputModule::attachScalarDofData_(writer,
                                                       *priVars[i],
                                                       oss.str());

            oss.str(""); oss << "delta_" << asImp_().primaryVarName(i);
            DiscBaseOutputModule::attachScalarDofData_(writer,
                                                       *delta[i],
                                                       oss.str());

            oss.str(""); oss << "weight_" << asImp_().primaryVarName(i);
            DiscBaseOutputModule::attachScalarDofData_(writer,
                                                       *priVarWeight[i],
                                                       oss.str());

            oss.str(""); oss << "defect_" << asImp_().eqName(i);
            DiscBaseOutputModule::attachScalarDofData_(writer,
                                                       *def[i],
                                                       oss.str());
        }

        asImp_().prepareOutputFields();
        asImp_().appendOutputFields(writer);
    }

    /*!
     * \brief Prepare the quantities relevant for the current solution
     *        to be appended to the output writers.
     */
    void prepareOutputFields() const
    {
        auto modIt = outputModules_.begin();
        const auto &modEndIt = outputModules_.end();
        for (; modIt != modEndIt; ++modIt) {
            (*modIt)->allocBuffers();
        }

        // iterate over grid
        ElementContext elemCtx(simulator_);

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateIntensiveQuantities(/*timeIdx=*/0);
            elemCtx.updateExtensiveQuantities(/*timeIdx=*/0);

            modIt = outputModules_.begin();
            for (; modIt != modEndIt; ++modIt)
                (*modIt)->processElement(elemCtx);
        }
    }

    /*!
     * \brief Append the quantities relevant for the current solution
     *        to an output writer.
     */
    void appendOutputFields(BaseOutputWriter &writer) const
    {
        auto modIt = outputModules_.begin();
        const auto &modEndIt = outputModules_.end();
        for (; modIt != modEndIt; ++modIt)
            (*modIt)->commitBuffers(writer);
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return gridView_; }

protected:
    /*!
     * \brief Finalize the initialization of the discretization
     *
     * This method requires the model to be constructed.
     */
    void finishInit_()
    {
        // initialize the volume of the finite volumes to zero
        int nDofs = asImp_().numDof();
        dofTotalVolume_.resize(nDofs);
        std::fill(dofTotalVolume_.begin(), dofTotalVolume_.end(), 0.0);

        ElementContext elemCtx(simulator_);

        // iterate through the grid and evaluate the initial condition
        ElementIterator it = gridView_.template begin</*codim=*/0>();
        const ElementIterator &eendit = gridView_.template end</*codim=*/0>();
        for (; it != eendit; ++it) {
            // ignore everything which is not in the interior if the
            // current process' piece of the grid
            if (it->partitionType() != Dune::InteriorEntity)
                continue;

            // deal with the current element
            elemCtx.updateStencil(*it);

            // loop over all element vertices, i.e. sub control volumes
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); dofIdx++)
            {
                // map the local degree of freedom index to the global one
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                dofTotalVolume_[globalDofIdx] += elemCtx.dofVolume(dofIdx, /*timeIdx=*/0);
            }
        }

        const auto sumHandle =
            GridCommHandleFactory::template sumHandle<double>(dofTotalVolume_,
                                                              asImp_().dofMapper());
        gridView_.communicate(*sumHandle,
                              Dune::InteriorBorder_InteriorBorder_Interface,
                              Dune::ForwardCommunication);

        localJacobian_.init(simulator_);
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(simulator_);
    }

    static bool storeIntensiveQuantities_()
    { return enableIntensiveQuantitiesCache_() || enableThermodynamicHints_(); }

    static bool enableIntensiveQuantitiesCache_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableIntensiveQuantityCache); }

    static bool enableThermodynamicHints_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableThermodynamicHints); }

    /*!
     * \brief Register all output modules which make sense for the model.
     *
     * This method is supposed to be overloaded by the actual models,
     * or else only the primary variables can be written to the result
     * files.
     */
    void registerOutputModules_()
    {
        // add the output modules available on all model
        auto *mod = new Ewoms::VtkPrimaryVarsModule<TypeTag>(simulator_);
        this->outputModules_.push_back(mod);
    }

    /*!
     * \brief Reference to the local residal object
     */
    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    /*!
     * \brief Find the degrees of freedoms adjacent to the grid boundary.
     */
    void updateBoundary_()
    {
        // resize the vectors and set everything to not being on the boundary
        onBoundary_.resize(asImp_().numDof());
        std::fill(onBoundary_.begin(), onBoundary_.end(), /*value=*/false);

        // loop over all elements of the grid
        Stencil stencil(gridView_);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            stencil.update(*elemIt);

            // do nothing if the element does not have boundary intersections
            if (stencil.numBoundaryFaces() == 0)
                continue;

            for (int dofIdx = 0; dofIdx < stencil.numPrimaryDof(); ++dofIdx) {
                int globalIdx = stencil.globalSpaceIndex(dofIdx);
                onBoundary_[globalIdx] = true;
            }
        }
    }

    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_.comm().rank() == 0; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Simulator &simulator_;

    // the representation of the spatial domain of the problem
    GridView gridView_;

    NewtonMethod newtonMethod_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    mutable SolutionVector solution_[historySize];
    mutable IntensiveQuantitiesVector intensiveQuantityCache_[historySize];
    mutable std::vector<bool> intensiveQuantityCacheUpToDate_[historySize];

    // all the index of the BoundaryTypes object for a vertex
    std::vector<bool> onBoundary_;

    std::list<BaseOutputModule<TypeTag>*> outputModules_;

    std::vector<double> dofTotalVolume_;
};
} // namespace Ewoms

#endif

