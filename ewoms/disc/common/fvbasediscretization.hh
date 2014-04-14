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
#include "fvbasevolumevariables.hh"
#include "fvbasefluxvariables.hh"

#include <ewoms/parallel/gridcommhandles.hh>
#include <ewoms/linear/nullborderlistcreator.hh>
#include <ewoms/common/timemanager.hh>

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
SET_TYPE_PROP(FvBaseDiscretization, TimeManager, Ewoms::TimeManager<TypeTag>);

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(FvBaseDiscretization, GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

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

SET_TYPE_PROP(FvBaseDiscretization, DiscVolumeVariables, Ewoms::FvBaseVolumeVariables<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, DiscFluxVariables, Ewoms::FvBaseFluxVariables<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(FvBaseDiscretization, GradientCalculator, Ewoms::FvBaseGradientCalculator<TypeTag>);

SET_TYPE_PROP(FvBaseDiscretization, DiscLocalJacobian, Ewoms::FvBaseLocalJacobian<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, LocalJacobian,
              typename GET_PROP_TYPE(TypeTag, DiscLocalJacobian));


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
SET_TYPE_PROP(FvBaseDiscretization, PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(FvBaseDiscretization, SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(FvBaseDiscretization, VolumeVariables, Ewoms::FvBaseVolumeVariables<TypeTag>);

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

// disable linearization recycling by default
SET_BOOL_PROP(FvBaseDiscretization, EnableLinearizationRecycling, false);

// disable partial relinearization by default
SET_BOOL_PROP(FvBaseDiscretization, EnablePartialRelinearization, false);

// disable constraints by default
SET_BOOL_PROP(FvBaseDiscretization, EnableConstraints, false);

// by default, disable the volume variable cache. If the volume
// variables are relatively cheap to calculate, the cache basically
// does not yield any performance impact because of the volume
// variable cache will cause additional pressure on the CPU caches...
SET_BOOL_PROP(FvBaseDiscretization, EnableVolumeVariablesCache, false);

// do not use thermodynamic hints by default. If you enable this, make
// sure to also enable the volume variable cache above to avoid
// getting an exception...
SET_BOOL_PROP(FvBaseDiscretization, EnableThermodynamicHints, false);

// if the deflection of the newton method is large, we do not need to
// solve the linear approximation accurately. Assuming that the value
// for the current solution is quite close to the final value, a
// reduction of 3 orders of magnitude in the defect should be
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
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GradientCalculator) GradientCalculator;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, DiscVtkBaseOutputModule) DiscVtkBaseOutputModule;
    typedef typename GET_PROP_TYPE(TypeTag, GridCommHandleFactory) GridCommHandleFactory;

    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),
        dim = GridView::dimension
    };

    typedef std::vector<VolumeVariables> VolumeVariablesVector;

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
    FvBaseDiscretization(Problem &problem)
        : problem_(problem)
        , gridView_(problem.gridView())
    {}

    ~FvBaseDiscretization()
    {
        // delete all VTK output modules
        auto modIt = vtkOutputModules_.begin();
        const auto &modEndIt = vtkOutputModules_.end();
        for (; modIt != modEndIt; ++modIt)
            delete *modIt;

        delete jacAsm_;
    }

    /*!
     * \brief Returns true iff a fluid phase is used by the model.
     *
     * \param phaseIdx The index of the fluid phase in question
     */
    bool phaseIsConsidered(int phaseIdx) const
    { return true; }

    /*!
     * \brief Register all run-time parameters for the model.
     */
    static void registerParameters()
    {
        JacobianAssembler::registerParameters();
        LocalJacobian::registerParameters();
        LocalResidual::registerParameters();
        GradientCalculator::registerParameters();
        VolumeVariables::registerParameters();
        FluxVariables::registerParameters();
        NewtonMethod::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkPrimaryVarsModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableThermodynamicHints, "Enable thermodynamic hints");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableVolumeVariablesCache, "Turn on caching of local quantities");
    }

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \copydetails Doxygen::problemParam
     */
    void init()
    {
        asImp_().updateBoundary_();

        int nDofs = asImp_().numDof();
        for (int timeIdx = 0; timeIdx < historySize; ++timeIdx) {
            solution_[timeIdx].resize(nDofs);

            if (storeVolumeVariables_()) {
                volVarsCache_[timeIdx].resize(nDofs);
                volVarsCacheUpToDate_[timeIdx].resize(nDofs, /*value=*/false);
            }
        }

        dofTotalVolume_.resize(nDofs);

        localJacobian_.init(problem_);
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem_);

        asImp_().applyInitialSolution_();
        asImp_().syncOverlap();

        // also set the solution of the "previous" time steps to the
        // initial solution.
        for (int timeIdx = 1; timeIdx < historySize; ++timeIdx)
            solution_[timeIdx] = solution_[/*timeIdx=*/0];

        asImp_().registerVtkModules_();
    }

    /*!
     * \brief Return the thermodynamic hint for a entity on the grid
     *        at given time.
     *
     * The hint is defined as a VolumeVariables object which is
     * supposed to be "close" to the VolumeVariables of the current
     * solution. It can be used as a good starting point for
     * non-linear solvers when having to solve non-linear relations
     * while updating the VolumeVariable. (This may yield a major
     * performance boost depending on how the physical models
     * require.)
     *
     * \attention If no up-to date volume variables are available, or
     *            if hints have been disabled, this method will return
     *            0.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is requested.
     * \param timeIdx The index used by the time discretization.
     */
    const VolumeVariables *thermodynamicHint(int globalIdx, int timeIdx) const
    {
        if (!enableThermodynamicHints_())
            return 0;

        if (volVarsCacheUpToDate_[timeIdx][globalIdx])
            return &volVarsCache_[timeIdx][globalIdx];

        // use the volume variables for the first up-to-date time index as hint
        for (int timeIdx2 = 0; timeIdx2 < historySize; ++timeIdx2)
            if (volVarsCacheUpToDate_[timeIdx2][globalIdx])
                return &volVarsCache_[timeIdx2][globalIdx];

        // no suitable up-to-date volume variables...
        return 0;
    }

    /*!
     * \brief Return the cached volume variables for a entity on the
     *        grid at given time.
     *
     * \attention If no up-to date volume variables are available,
     *            this method will return 0.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is requested.
     * \param timeIdx The index used by the time discretization.
     */
    const VolumeVariables *cachedVolumeVariables(int globalIdx, int timeIdx) const
    {
        if (!enableVolumeVariablesCache_() ||
            !volVarsCacheUpToDate_[timeIdx][globalIdx])
            return 0;

        return &volVarsCache_[timeIdx][globalIdx];
    }

    /*!
     * \brief Update the volume variable cache for a entity on the grid at given time.
     *
     * \param volVars The VolumeVariables object hint for a given degree of freedom.
     * \param globalIdx The global space index for the entity where a
     *                  hint is to be set.
     * \param timeIdx The index used by the time discretization.
     */
    void updateCachedVolumeVariables(const VolumeVariables &volVars,
                                     int globalIdx,
                                     int timeIdx) const
    {
        if (!storeVolumeVariables_())
            return;

        volVarsCache_[timeIdx][globalIdx] = volVars;
        volVarsCacheUpToDate_[timeIdx][globalIdx] = true;
    }

    /*!
     * \brief Invalidate the cache for a given volume variables object.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is to be set.
     * \param timeIdx The index used by the time discretization.
     */
    void invalidateVolumeVariablesCacheEntry(int globalIdx,
                                             int timeIdx) const
    {
        if (!storeVolumeVariables_())
            return;

        volVarsCacheUpToDate_[timeIdx][globalIdx] = false;
    }

    /*!
     * \brief Move the volume variables for a given time index to the back.
     *
     * This method should only be called by the time discretization.
     *
     * \param numSlots The number of time step slots for which the
     *                 hints should be shifted.
     */
    void shiftVolumeVariablesCache(int numSlots = 1)
    {
        if (!storeVolumeVariables_())
            return;

        for (int timeIdx = 0; timeIdx < historySize - numSlots; ++ timeIdx) {
            volVarsCache_[timeIdx + numSlots] = volVarsCache_[timeIdx];
            volVarsCacheUpToDate_[timeIdx + numSlots] = volVarsCacheUpToDate_[timeIdx];
        }

        // invalidate the cache for the most recent time index
        std::fill(volVarsCacheUpToDate_[/*timeIdx=*/0].begin(),
                  volVarsCacheUpToDate_[/*timeIdx=*/0].end(),
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

        ElementContext elemCtx(problem_);
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
    void globalStorage(EqVector &storage)
    {
        storage = 0;

        ElementContext elemCtx(problem_);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateVolVars(/*timeIdx=*/0);
            localResidual().evalStorage(elemCtx, /*timeIdx=*/0);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                storage += localResidual().storageTerm()[i];
        };

        storage = gridView_.comm().sum(storage);
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
     * \param globalVertexIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar absPv = std::abs(asImp_().solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx]);
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
            Valgrind::CheckDefined(asImp_().solution(/*timeIdx=*/0)[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.apply();
        if (converged) {
            asImp_().syncOverlap();
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
        volVarsCache_[/*timeIdx=*/0] = volVarsCache_[/*timeIdx=*/1];
        volVarsCacheUpToDate_[/*timeIdx=*/0] = volVarsCacheUpToDate_[/*timeIdx=*/1];

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

        // shift the volume variables cache by one position in the
        // history
        asImp_().shiftVolumeVariablesCache(/*numSlots=*/1);
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
    { return problem_.vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_.elementMapper(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler;
        jacAsm_->init(problem_);
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
     *        element given the complete set of volume variables
     *
     * \copydetails Doxygen::ecfvElemCtxParam
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    { }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param writer  The VTK multi writer object on which the fields should be added.
     * \param u       The solution function
     * \param deltaU  The delta of the solution function before and after the Newton update
     */
    template <class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
                                 const SolutionVector &u,
                                 const GlobalEqVector &deltaU) const
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        GlobalEqVector globalResid(u.size());
        asImp_().globalResidual(globalResid, u);

        // create the required scalar fields
        unsigned numDof = asImp_().numDof();

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* priVars[numEq];
        ScalarField* priVarWeight[numEq];
        ScalarField* relError = writer.allocateManagedBuffer(numDof);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            priVars[pvIdx] = writer.allocateManagedBuffer(numDof);
            priVarWeight[pvIdx] = writer.allocateManagedBuffer(numDof);
            delta[pvIdx] = writer.allocateManagedBuffer(numDof);
            def[pvIdx] = writer.allocateManagedBuffer(numDof);
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

        DiscVtkBaseOutputModule::attachDofData_(writer, *relError, "relErr", /*numComponents=*/1);

        for (int i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss.str(""); oss << "priVar_" << asImp_().primaryVarName(i);
            DiscVtkBaseOutputModule::attachDofData_(writer,
                                                    *priVars[i],
                                                    oss.str(),
                                                    /*numComponents=*/1);

            oss.str(""); oss << "delta_" << asImp_().primaryVarName(i);
            DiscVtkBaseOutputModule::attachDofData_(writer,
                                                    *delta[i],
                                                    oss.str(),
                                                    /*numComponents=*/1);

            oss.str(""); oss << "weight_" << asImp_().primaryVarName(i);
            DiscVtkBaseOutputModule::attachDofData_(writer,
                                                    *priVarWeight[i],
                                                    oss.str(),
                                                    /*numComponents=*/1);

            oss.str(""); oss << "defect_" << asImp_().eqName(i);
            DiscVtkBaseOutputModule::attachDofData_(writer,
                                                    *def[i],
                                                    oss.str(),
                                                    /*numComponents=*/1);
        }

        asImp_().addOutputVtkFields(writer);
    }

    /*!
     * \brief Append the quantities relevant for the current solution to a VTK multi writer.
     *
     * This should be overwritten by the actual model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param writer The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer) const
    {
        auto modIt = vtkOutputModules_.begin();
        const auto &modEndIt = vtkOutputModules_.end();
        for (; modIt != modEndIt; ++modIt)
            (*modIt)->allocBuffers(writer);

        // iterate over grid
        ElementContext elemCtx(problem_);

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateVolVars(/*timeIdx=*/0);
            elemCtx.updateFluxVars(/*timeIdx=*/0);

            modIt = vtkOutputModules_.begin();
            for (; modIt != modEndIt; ++modIt)
                (*modIt)->processElement(elemCtx);
        }

        modIt = vtkOutputModules_.begin();
        for (; modIt != modEndIt; ++modIt)
            (*modIt)->commitBuffers(writer);
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_.gridView(); }

protected:
    static bool storeVolumeVariables_()
    { return enableVolumeVariablesCache_() || enableThermodynamicHints_(); }

    static bool enableVolumeVariablesCache_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableVolumeVariablesCache); }

    static bool enableThermodynamicHints_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableThermodynamicHints); }

    /*!
     * \brief Register all VTK output modules which make sense for the model.
     *
     * This method is supposed to be overloaded by the actual models,
     * or else only the primary variables can be written to the result
     * files.
     */
    void registerVtkModules_()
    {
        // add the VTK output modules available on all model
        auto *vtkMod = new Ewoms::VtkPrimaryVarsModule<TypeTag>(problem_);
        this->vtkOutputModules_.push_back(vtkMod);
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
     * \brief Applies the initial solution for all vertices of the grid.
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        SolutionVector &uCur = asImp_().solution(/*timeIdx=*/0);
        uCur = Scalar(0.0);

        // initialize the volume of the FV boxes to zero
        std::fill(dofTotalVolume_.begin(), dofTotalVolume_.end(), 0.0);

        ElementContext elemCtx(problem_);

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
                problem_.initial(uCur[globalIdx], elemCtx, dofIdx, /*timeIdx=*/0);
                Valgrind::CheckDefined(uCur[globalIdx]);

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

        if (storeVolumeVariables_()) {
            // invalidate all cached volume variables
            for (int timeIdx = 0; timeIdx < historySize; ++ timeIdx) {
                std::fill(volVarsCacheUpToDate_[timeIdx].begin(),
                          volVarsCacheUpToDate_[timeIdx].end(),
                          false);
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
    Problem &problem_;

    // the representation of the spatial domain of the problem
    GridView gridView_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    mutable SolutionVector solution_[historySize];
    mutable VolumeVariablesVector volVarsCache_[historySize];
    mutable std::vector<bool> volVarsCacheUpToDate_[historySize];

    // all the index of the BoundaryTypes object for a vertex
    std::vector<bool> onBoundary_;

    std::list<VtkBaseOutputModule<TypeTag>*> vtkOutputModules_;

    std::vector<double> dofTotalVolume_;
};
} // namespace Ewoms

#endif

