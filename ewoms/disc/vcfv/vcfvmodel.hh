// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
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
 * \copydoc Ewoms::VcfvModel
 */
#ifndef EWOMS_VCFV_MODEL_HH
#define EWOMS_VCFV_MODEL_HH

#include "vcfvelementcontext.hh"
#include "vcfvlocaljacobian.hh"
#include "vcfvlocalresidual.hh"
#include "vcfvproperties.hh"
#include "vcfvpropertydefaults.hh"

#include <ewoms/vtk/vcfvvtkoutputmodule.hh>
#include <ewoms/vtk/vcfvvtkprimaryvarsmodule.hh>
#include <ewoms/parallel/gridcommhandles.hh>

#include <dune/common/fvector.hh>

#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace Ewoms {

/*!
 * \ingroup VcfvModel
 *
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class VcfvModel
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

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::ctype CoordScalar;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;
#else
    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;
#endif

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying a model is not a good idea
    VcfvModel(const VcfvModel &);

public:
    // this constructor required to be explicitly specified because
    // we've defined a constructor above which deletes all implicitly
    // generated constructors in C++.
    VcfvModel()
    {}

    ~VcfvModel()
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
        VolumeVariables::registerParameters();
        FluxVariables::registerParameters();
        NewtonMethod::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VcfvVtkPrimaryVarsModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableHints, "Enable thermodynamic hints");
    }

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \copydetails Doxygen::problemParam
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        updateBoundary_();

        int nDofs = asImp_().numDofs();
        for (int timeIdx = 0; timeIdx < historySize; ++timeIdx)
            solution_[timeIdx].resize(nDofs);
        boxVolume_.resize(nDofs);

        localJacobian_.init(problem_());
        jacAsm_ = new JacobianAssembler();
        jacAsm_->init(problem_());

        asImp_().applyInitialSolution_();

        // resize the hint vectors
        if (enableHints_()) {
            int nVerts = gridView_().size(dim);
            for (int timeIdx = 0; timeIdx < historySize; ++timeIdx) {
                hints_[timeIdx].resize(nVerts);
                hintsUsable_[timeIdx].resize(nVerts);
                std::fill(hintsUsable_[timeIdx].begin(),
                          hintsUsable_[timeIdx].end(),
                          false);
            }
        }

        // also set the solution of the "previous" time step to the
        // initial solution.
        solution_[/*timeIdx=*/1] = solution_[/*timeIdx=*/0];

        asImp_().registerVtkModules_();
    }

    /*!
     * \brief Return the hint for a entity on the grid at given time.
     *
     * The hint is defined as a VolumeVariables object which is
     * supposed to be "close" to the VolumeVariables of the current
     * solution. It can be used as a good starting point for
     * non-linear solvers when having to solve non-linear relations
     * while updating the VolumeVariable. (This may yield a major
     * performance boost depending on how the physical models
     * require.)
     *
     * \attention If no hint is available, this method will return 0.
     *
     * \param globalIdx The global space index for the entity where a
     *                  hint is requested.
     * \param timeIdx The index used by the time discretization.
     */
    const VolumeVariables *hint(int globalIdx, int timeIdx) const
    {
        if (!enableHints_() ||
            !hintsUsable_[timeIdx][globalIdx])
        {
            return 0;
        }

        return &hints_[timeIdx][globalIdx];
    }

    /*!
     * \brief Update the hint for a entity on the grid at given time.
     *
     * \param hint The VolumeVariables object which ought to serve as
     *             the hint for a given entity.
     * \param globalIdx The global space index for the entity where a
     *                  hint is to be set.
     * \param timeIdx The index used by the time discretization.
     */
    void setHint(const VolumeVariables &hint,
                 int globalIdx,
                 int timeIdx) const
    {
        if (!enableHints_())
            return;

        hints_[timeIdx][globalIdx] = hint;
        hintsUsable_[timeIdx][globalIdx] = true;
    }

    /*!
     * \brief Move the hints for a given time index to the back.
     *
     * This method should only be called by the time discretization.
     *
     * \param numSlots The number of time step slots for which the
     *                 hints should be shifted.
     */
    void shiftHints(int numSlots = 1)
    {
        for (int timeIdx = 0; timeIdx < historySize - numSlots; ++ timeIdx)
            hints_[timeIdx + numSlots] = hints_[timeIdx];
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
        SolutionVector tmp(solution(/*timeIdx=*/0));
        solution_[/*timeIdx=*/0] = u;
        Scalar res = globalResidual(dest);
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

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateAll(*elemIt);
            residual.resize(elemCtx.numScv());
            storageTerm.resize(elemCtx.numScv());
            localResidual().eval(residual, storageTerm, elemCtx);

            int numScv = elemCtx.numScv();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx) {
                int globalI = vertexMapper().map(*elemIt, scvIdx, dim);
                dest[globalI] += residual[scvIdx];
            }
        };

        // add up the residuals on the process borders
        GridCommHandleSum<EqVector, GlobalEqVector, VertexMapper, /*commCodim=*/dim>
            sumHandle(dest, vertexMapper());
        gridView().communicate(sumHandle,
                               Dune::InteriorBorder_InteriorBorder_Interface,
                               Dune::ForwardCommunication);

        // calculate the square norm of the residual. this is not
        // entirely correct, since the residual for the finite volumes
        // which are on the boundary are counted once for every
        // process. As often in life: shit happens (, we don't care)...
        Scalar result2 = dest.two_norm2();
        result2 = gridView().comm().sum(result2);

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

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);
            localResidual().evalStorage(elemCtx, /*timeIdx=*/0);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                storage += localResidual().storageTerm()[i];
        };

        storage = gridView_().comm().sum(storage);
    }

    /*!
     * \brief Returns the volume \f$\mathrm{[m^3]}\f$ of a given control volume.
     *
     * \param globalIdx The global index of the control volume's
     *                  associated vertex
     */
    Scalar boxVolume(int globalIdx) const
    { return boxVolume_[globalIdx][0]; }

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
     * \copydoc jacobianAssembler() const
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
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual() const
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar absPv = std::abs(this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx]);
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
    Scalar relativeErrorVertex(int vertexIdx,
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
        for (size_t i = 0; i < solution(/*timeIdx=*/0).size(); ++i)
            Valgrind::CheckDefined(solution(/*timeIdx=*/0)[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.apply();
        if (converged)
            asImp_().updateSuccessful();
        else
            asImp_().updateFailed();


#if HAVE_VALGRIND
        for (size_t i = 0; i < solution(/*timeIdx=*/0).size(); ++i) {
            Valgrind::CheckDefined(solution(/*timeIdx=*/0)[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }

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
        hints_[/*timeIdx=*/0] = hints_[/*timeIdx=*/1];

        solution_[/*timeIdx=*/0] = solution_[/*timeIdx=*/1];
        jacAsm_->reassembleAll();
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

        // shift the hints by one position in the history
        shiftHints();
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
    { res.template serializeEntities<dim>(asImp_(), this->gridView_()); }

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
        res.template deserializeEntities<dim>(asImp_(), this->gridView_());
        solution_[/*timeIdx=*/1] = solution_[/*timeIdx=*/0];
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param vert The DUNE Codim<dim> entity which's data should be
     *             serialized
     */
    void serializeEntity(std::ostream &outstream,
                         const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);

        // write phase state
        if (!outstream.good()) {
            OPM_THROW(std::runtime_error,
                       "Could not serialize vertex "
                       << vertIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << solution_[/*timeIdx=*/0][vertIdx][eqIdx] << " ";
        }
    }

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param vertex The DUNE Codim<dim> entity which's data should be
     *               deserialized
     */
    void deserializeEntity(std::istream &instream,
                           const Vertex &vertex)
    {
        int vertIdx = dofMapper().map(vertex);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                OPM_THROW(std::runtime_error,
                           "Could not deserialize vertex "
                           << vertIdx);
            instream >> solution_[/*timeIdx=*/0][vertIdx][eqIdx];
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    { return gridView_().size(dim); }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * This usually means a mapper for vertices.
     */
    const DofMapper &dofMapper() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        delete jacAsm_;
        jacAsm_ = new JacobianAssembler;
        jacAsm_->init(problem_());
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
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "vcfv"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     *
     * \param pvIdx The index of the primary variable of interest.
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;
        oss << pvIdx;
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
        oss << eqIdx;
        return oss.str();
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \copydetails Doxygen::vcfvElemCtxParam
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
        unsigned numVertices = this->gridView_().size(dim);
        //unsigned numElements = this->gridView_().size(0);

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* priVars[numEq];
        ScalarField* priVarWeight[numEq];
        ScalarField* relError = writer.allocateManagedBuffer(numVertices);
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            priVars[pvIdx] = writer.allocateManagedBuffer(numVertices);
            priVarWeight[pvIdx] = writer.allocateManagedBuffer(numVertices);
            delta[pvIdx] = writer.allocateManagedBuffer(numVertices);
            def[pvIdx] = writer.allocateManagedBuffer(numVertices);
        }

        VertexIterator vIt = this->gridView_().template begin<dim>();
        VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt != vEndIt; ++ vIt)
        {
            int globalIdx = vertexMapper().map(*vIt);
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                (*priVars[pvIdx])[globalIdx] = u[globalIdx][pvIdx];
                (*priVarWeight[pvIdx])[globalIdx] = asImp_().primaryVarWeight(globalIdx, pvIdx);
                (*delta[pvIdx])[globalIdx] =
                    - deltaU[globalIdx][pvIdx];
                (*def[pvIdx])[globalIdx] = globalResid[globalIdx][pvIdx];
            }

            PrimaryVariables uOld(u[globalIdx]);
            PrimaryVariables uNew(uOld);
            uNew -= deltaU[globalIdx];
            (*relError)[globalIdx] = asImp_().relativeErrorVertex(globalIdx,
                                                                  uOld,
                                                                  uNew);
        }

        writer.attachVertexData(*relError, "relative Error");

        for (int i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss.str(""); oss << "priVar_" << asImp_().primaryVarName(i);
            writer.attachVertexData(*priVars[i], oss.str());

            oss.str(""); oss << "delta_" << asImp_().primaryVarName(i);
            writer.attachVertexData(*delta[i], oss.str());

            oss.str(""); oss << "weight_" << asImp_().primaryVarName(i);
            writer.attachVertexData(*priVarWeight[i], oss.str());

            oss.str(""); oss << "defect_" << asImp_().eqName(i);
            writer.attachVertexData(*def[i], oss.str());
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
        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);
            elemCtx.updateScvfVars(/*timeIdx=*/0);

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
    { return problem_().gridView(); }

protected:
    static bool enableHints_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableHints); }

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
        auto *vtkMod = new Ewoms::VcfvVtkPrimaryVarsModule<TypeTag>(this->problem_());
        this->vtkOutputModules_.push_back(vtkMod);
    }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Reference to the local residal object
     */
    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    /*!
     * \brief Updates the stuff which determines a vertex' or
     *        element's boundary type
     */
    void updateBoundary_()
    {
        // resize the vectors
        onBoundary_.resize(numDofs());

        // loop over all elements of the grid
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            // do nothing if the element does not have boundary intersections
            if (!elemIt->hasBoundaryIntersections())
                continue;

            // retrieve the reference element for the current element
            const Element &elem = *elemIt;
            Dune::GeometryType geoType = elem.geometry().type();
            const ReferenceElement &refElement = ReferenceElements::general(geoType);

            // loop over all intersections of the element
            IntersectionIterator isIt = gridView_().ibegin(elem);
            const IntersectionIterator &endIt = gridView_().iend(elem);
            for (; isIt != endIt; ++isIt)
            {
                // do nothing if the face is _not_ on the boundary
                if (!isIt->boundary())
                    continue;

                // loop over all vertices of the intersection
                int faceIdx = isIt->indexInInside();
                int numFaceVerts = refElement.size(faceIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    // find the local element index of the face's
                    // vertex
                    int scvIdx = refElement.subEntity(/*entityIdx=*/faceIdx,
                                                   /*entityCodim=*/1,
                                                   /*subEntityIdx=*/faceVertIdx,
                                                   /*subEntityCodim=*/dim);
                    int globalIdx = vertexMapper().map(*elemIt, scvIdx, /*codim=*/dim);

                    onBoundary_[globalIdx] = true;
                } // loop over intersection's vertices
            } // loop over intersections
        } // loop over elements
    }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        SolutionVector &uCur = solution(/*timeIdx=*/0);
        uCur = Scalar(0.0);
        boxVolume_ = Scalar(0.0);

        ElementContext elemCtx(this->problem_());

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator it = gridView_().template begin<0>();
        const ElementIterator &eendit = gridView_().template end<0>();
        for (; it != eendit; ++it) {
            // deal with the current element
            elemCtx.updateFVElemGeom(*it);

            // loop over all element vertices, i.e. sub control volumes
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
            {
                // map the local vertex index to the global one
                int globalIdx = vertexMapper().map(*it,
                                                   scvIdx,
                                                   dim);

                // let the problem do the dirty work of nailing down
                // the initial solution.
                problem_().initial(uCur[globalIdx],
                                   elemCtx,
                                   scvIdx,
                                   /*timeIdx=*/0);
                Valgrind::CheckDefined(uCur[globalIdx]);

                boxVolume_[globalIdx] +=
                    elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[scvIdx].volume;
            }
        }

        // add the finite volumes which cross process borders
        GridCommHandleSum<Dune::FieldVector<Scalar, 1>,
                          Dune::BlockVector<Dune::FieldVector<Scalar, 1> >,
                          VertexMapper,
                          /*commCodim=*/dim> sumVolumeHandle(boxVolume_, vertexMapper());
        gridView().communicate(sumVolumeHandle,
                               Dune::InteriorBorder_InteriorBorder_Interface,
                               Dune::ForwardCommunication);
    }

    // the hint cache for the previous and the current volume
    // variables
    mutable std::vector<bool> hintsUsable_[historySize];
    mutable std::vector<VolumeVariables> hints_[historySize];

    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    JacobianAssembler *jacAsm_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    mutable SolutionVector solution_[historySize];

    // all the index of the BoundaryTypes object for a vertex
    std::vector<bool> onBoundary_;

    std::list<VcfvVtkOutputModule<TypeTag>*> vtkOutputModules_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;
};
}

#endif
