// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Base class for models using box discretization
 */
#ifndef DUMUX_BOX_MODEL_HH
#define DUMUX_BOX_MODEL_HH

#include "boxproperties.hh"
#include "boxpropertydefaults.hh"

#include "boxelementcontext.hh"
#include "boxlocaljacobian.hh"
#include "boxlocalresidual.hh"
#include <dumux/boxmodels/vtk/boxvtkoutputmodule.hh>
#include <dumux/boxmodels/vtk/boxvtkprimaryvarsmodule.hh>

#include <dumux/parallel/vertexhandles.hh>

#include <sstream>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class BoxModel
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

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::ctype CoordScalar;
    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    // copying a model is not a good idea
    BoxModel(const BoxModel &);

public:
    /*!
     * \brief The constructor.
     */
    BoxModel()
    {
    }

    ~BoxModel()
    {
        // delete all VTK output modules
        auto modIt = vtkOutputModules_.begin();
        const auto &modEndIt = vtkOutputModules_.end();
        for (; modIt != modEndIt; ++modIt)
            delete *modIt;

        delete jacAsm_;
    }

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param prob The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &prob)
    {
        problemPtr_ = &prob;

        updateBoundaryTypes_();

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

    const VolumeVariables *hint(int globalIdx, int timeIdx) const
    {
        if (!enableHints_() ||
            !hintsUsable_[timeIdx][globalIdx])
        {
            return 0;
        }

        return &hints_[timeIdx][globalIdx];
    }

    void setHint(const VolumeVariables &hint,
                 int globalIdx,
                 int timeIdx) const
    {
        if (!enableHints_())
            return;

        hints_[timeIdx][globalIdx] = hint;
        hintsUsable_[timeIdx][globalIdx] = true;
    };

    void shiftHints(int numSlots = 1)
    {
        for (int timeIdx = 0; timeIdx < historySize - numSlots; ++ timeIdx)
            hints_[timeIdx + numSlots] = hints_[timeIdx];
    };

    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param dest Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(GlobalEqVector &dest,
                          const SolutionVector &u)
    {
        SolutionVector tmp(solution(/*timeIdx=*/0));
        solution(/*timeIdx=*/0) = u;
        Scalar res = globalResidual(dest);
        solution(/*timeIdx=*/0) = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param dest Stores the result
     */
    Scalar globalResidual(GlobalEqVector &dest)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateAll(*elemIt);
            localResidual().eval(elemCtx);

            int numScv = elemCtx.numScv();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx) {
                int globalI = vertexMapper().map(*elemIt, scvIdx, dim);
                dest[globalI] += localResidual().residual(scvIdx);
            }
        };

        // calculate the square norm of the residual
        Scalar result2 = dest.two_norm2();
        result2 = gridView().comm().sum(result2);

        // add up the residuals on the process borders
        if (gridView().comm().size() > 1) {
            VertexHandleSum<EqVector, GlobalEqVector, VertexMapper>
                sumHandle(dest, vertexMapper());
            gridView().communicate(sumHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param dest Stores the result
     */
    void globalStorage(EqVector &dest)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = gridView_().template begin<0>();
        const ElementIterator elemEndIt = gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);
            localResidual().evalStorage(elemCtx, /*timeIdx=*/0);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                dest += localResidual().storageTerm()[i];
        };

        dest = gridView_().comm().sum(dest);
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
     */
    const SolutionVector &solution(int timeIdx) const
    { return solution_[timeIdx]; }

    /*!
     * \brief Reference to the solution at a given history index as a block vector.
     */
    SolutionVector &solution(int timeIdx)
    { return solution_[timeIdx]; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianAssembler()
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    LocalJacobian &localJacobian()
    { return localJacobian_; }
    /*!
     * \copydoc localJacobian()
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }

    /*!
     * \brief Returns the local residual function.
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual()
     */
    const LocalResidual &localResidual() const
    { return localJacobian().localResidual(); }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param vertIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int vertIdx, int pvIdx) const
    {
        Scalar absPv = std::abs(this->solution(/*timeIdx=*/1)[vertIdx][pvIdx]);
        return 1.0/std::max(absPv, 1.0);
    }

    /*!
     * \brief Returns the relative error between two vectors of
     *        primary variables.
     *
     * \param vertexIdx The global index of the control volume's
     *                  associated vertex
     * \param pv1 The first vector of primary variables
     * \param pv2 The second vector of primary variables
     *
     * \todo The vertexIdx argument is pretty hacky. it is required by
     *       models with pseudo primary variables (i.e. the primary
     *       variable switching models). the clean solution would be
     *       to access the pseudo primary variables from the primary
     *       variables.
     */
    Scalar relativeErrorVertex(int vertexIdx,
                               const PrimaryVariables &pv1,
                               const PrimaryVariables &pv2)
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
     * \param controller The controller which specifies the behaviour
     *                   of the non-linear solver
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < solution(/*timeIdx=*/0).size(); ++i)
            Valgrind::CheckDefined(solution(/*timeIdx=*/0)[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        bool converged = solver.execute(controller);
        if (converged) {
            asImp_().updateSuccessful();
        }
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
    {
        updateBoundaryTypes_();
    }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    { };

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
    };

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
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << solution_[/*timeIdx=*/0][vertIdx][eqIdx] << " ";
        }
    };

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param vert The DUNE Codim<dim> entity which's data should be
     *             deserialized
     */
    void deserializeEntity(std::istream &instream,
                           const Vertex &vert)
    {
        int vertIdx = dofMapper().map(vert);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << vertIdx);
            instream >> solution_[/*timeIdx=*/0][vertIdx][eqIdx];
        }
    };

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
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); };

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); };

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
     * \brief Return a pointer to the BoundaryTypes for a given global
     *        vertex index or 0 if the vertex is not on the boundary.
     */
    const BoundaryTypes *boundaryTypes(int globalIdx) const
    {
        static BoundaryTypes dummy;
        int bvertIdx = boundaryVertexIndex_[globalIdx];
        if (bvertIdx < 0)
            return &dummy;
        return &boundaryTypes_[bvertIdx];
    }

    /*!
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "box"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        oss << pvIdx;
        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
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
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    { };

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
                                 const GlobalEqVector &deltaU)
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

        asImp_().addOutputVtkFields(u, writer);
    }

    /*!
     * \brief Add the quantities of a time step which ought to be written to disk.
     *
     * This should be overwritten by the actual model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param sol The global vector of primary variable values.
     * \param writer The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
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
            elemCtx.updateAll(*elemIt);

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
    { return GET_PARAM(TypeTag, bool, EnableHints); }

    void registerVtkModules_()
    {
        // add the VTK output modules available on all model
        auto *vtkMod = new Dumux::BoxVtkPrimaryVarsModule<TypeTag>(this->problem_());
        this->vtkOutputModules_.push_back(vtkMod);
    };

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
    void updateBoundaryTypes_()
    {
        // resize the vectors
        boundaryVertexIndex_.resize(numDofs());
        std::fill(boundaryVertexIndex_.begin(),
                  boundaryVertexIndex_.end(),
                  -1);

        int numBoundaryVertices = 0;
        BoxBoundaryContext<TypeTag> boundaryCtx(problem_());

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
            const ReferenceElement &refElem = ReferenceElements::general(geoType);

            boundaryCtx.update(*elemIt);

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
                int numFaceVerts = refElem.size(faceIdx, 1, dim);
                for (int faceVertIdx = 0;
                     faceVertIdx < numFaceVerts;
                     ++faceVertIdx)
                {
                    // find the local element index of the face's
                    // vertex
                    int scvIdx = refElem.subEntity(/*entityIdx=*/faceIdx,
                                                   /*entityCodim=*/1,
                                                   /*subEntityIdx=*/faceVertIdx,
                                                   /*subEntityCodim=*/dim);
                    int globalIdx = vertexMapper().map(*elemIt, scvIdx, /*codim=*/dim);
                    if (boundaryVertexIndex_[globalIdx] >= 0)
                        continue; // vertex has already been visited

                    // add a BoundaryTypes object
                    if (boundaryTypes_.size() <= numBoundaryVertices)
                        boundaryTypes_.resize(numBoundaryVertices + 1);
                    BoundaryTypes &bTypes = boundaryTypes_[numBoundaryVertices];
                    boundaryVertexIndex_[globalIdx] = numBoundaryVertices;
                    ++numBoundaryVertices;

                    problem_().boundaryTypes(bTypes, boundaryCtx, scvIdx, /*timeIdx=*/0);
                } // loop over intersection's vertices
            } // loop over intersections
        } // loop over elements
    };

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

        // add up the primary variables and the volumes of the boxes
        // which cross process borders
        if (gridView().comm().size() > 1) {
            VertexHandleSum<Dune::FieldVector<Scalar, 1>,
                Dune::BlockVector<Dune::FieldVector<Scalar, 1> >,
                VertexMapper> sumVolumeHandle(boxVolume_, vertexMapper());
            gridView().communicate(sumVolumeHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }
    }

    // the hint cache for the previous and the current volume
    // variables
    mutable std::vector<bool> hintsUsable_[historySize];
    mutable std::vector<VolumeVariables> hints_[historySize];

    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; };

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
    SolutionVector solution_[historySize];

    // all the index of the BoundaryTypes object for a vertex
    std::vector<int> boundaryVertexIndex_;
    std::vector<BoundaryTypes> boundaryTypes_;

    std::list<BoxVtkOutputModule<TypeTag>*> vtkOutputModules_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;
};
}

#endif
