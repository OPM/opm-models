// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Caculates the jacobian of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_JACOBIAN_HH
#define DUMUX_BOX_JACOBIAN_HH

#include <dune/common/exceptions.hh>
#include <dune/disc/operators/localstiffness.hh>

#include <dumux/auxiliary/valgrind.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <boost/format.hpp>

#include "boxproperties.hh"


namespace Dune
{
/*!
 * \ingroup BoxScheme
 * \brief Element-wise caculation of the jacobian matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation = typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) >
class BoxJacobian : public Dune::LocalStiffness<typename GET_PROP_TYPE(TypeTag, PTAG(GridView)),
                                                typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                                GET_PROP_VALUE(TypeTag, PTAG(NumEq)) >
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))    Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef BoxJacobian<TypeTag, Implementation> ThisType;
    enum {
        numEq     = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        dim       = GridView::dimension,
        dimWorld  = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype                CoordScalar;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity               Element;
    typedef typename GridView::template Codim<0>::Iterator             ElementIterator;
    typedef typename Element::EntityPointer                            ElementPointer;

    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;
    typedef typename RefElemProp::ReferenceElement              ReferenceElement;

    typedef typename GridView::IntersectionIterator                    IntersectionIterator;
    typedef typename Element::Geometry                                 Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry))   FVElementGeometry;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::DofEntityMapper         DofEntityMapper;
    typedef typename SolutionTypes::VertexMapper            VertexMapper;
    typedef typename SolutionTypes::ElementMapper           ElementMapper;
    typedef typename SolutionTypes::Solution                Solution;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;
    typedef typename SolutionTypes::JacobianAssembler       JacobianAssembler;


    typedef typename GET_PROP(TypeTag, PTAG(VertexData))::type VertexData;
    typedef typename std::vector<VertexData>                   VertexDataArray;
    
public:
    BoxJacobian(Problem &problem)
        : problem_(problem),
          gridView_(problem.gridView()),
          curElementPtr_(* ++gridView_.template begin<0>())
    {
    }

    /*!
     * \brief Compute the global residual right hand side
     *        of an equation we would like to have zero.
     */
    void evalGlobalResidual(SolutionFunction &residual)
    {
        *residual = 0;
        SolutionOnElement tmpSol, tmpSolOld;
        SolutionOnElement localResid;
        localResid.resize(12);
        ElementIterator elemIt = gridView_.template begin<0>();
        const ElementIterator elemEndIt = gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            this->setCurrentElement(*elemIt);
            this->restrictToElement(tmpSol, this->problem().model().curSolFunction());
            this->restrictToElement(tmpSolOld, this->problem().model().prevSolFunction());
            this->asImp_().setCurrentSolution(tmpSol);
            this->asImp_().setPreviousSolution(tmpSolOld);

            asImp_().evalLocalResidual(localResid, true); // with boundary

            for (int i = 0; i < elemIt->template count<dim>(); ++i) {
                int globalI = this->problem().model().vertexMapper().map(*elemIt, i, dim);
                (*residual)[globalI] += localResid[i];
            }
        };
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element, given a local solution 'localU'.
     */
    void assemble(const Element &element, const SolutionOnElement& localU, int orderOfShapeFns = 1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

#if HAVE_VALGRIND
        for (size_t i = 0; i < localU.size(); ++i)
            Valgrind::CheckDefined(localU[i]);
#endif // HAVE_VALGRIND

        SolutionOnElement mutableLocalU(localU);
        asImp_().assemble_(element, mutableLocalU);
    }

    /*!
     * \brief Assemble the linear system of equations for the
     *        verts of a element.
     */
    void assemble(const Element &element, int orderOfShapeFns = 1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        int numVertices = curElementGeom_.numVertices;
        SolutionOnElement localU(numVertices);
        restrictToElement(localU, problem_.model().curSolFunction());
        asImp_().assemble_(element, localU);
    }
  
    /*!
     * \brief Express the boundary conditions for a element in terms
     *        of a linear equation system.
     */
    void assembleBoundaryCondition(const Element &element, int orderOfShapeFns=1)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        // reset the right hand side and the boundary
        // condition type vector
        asImp_().resetRhs_();

        // set the boundary types
        asImp_().updateBoundaryTypes_();
        
        // apply the neumann conditions
        asImp_().applyBoundaryCondition_();
    };

    void updateBoundaryTypes(const Element &element)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);
        asImp_().updateBoundaryTypes_();
    }

    void applyBoundaryCondition_()
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        // evaluate boundary conditions for all intersections of
        // the current element
        IntersectionIterator isIt = curElement_().ileafbegin();
        const IntersectionIterator &endIt = curElement_().ileafend();
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on boundaries.
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all verts of the
            // current face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);

                // handle boundary conditions for a single
                // sub-control volume face
                applyBoundaryCondition_(isIt,
                                        elemVertIdx,
                                        boundaryFaceIdx);
            }
        }
    }

    // handle boundary conditions for a single
    // sub-control volume face
    void applyBoundaryCondition_(const IntersectionIterator &isIt,
                                 int scvIdx,
                                 int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundaries
        PrimaryVarVector values(0.0);
        bool wasEvaluated = false;

        // loop over all primary variables to deal with mixed
        // boundary conditions
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (this->bctype[scvIdx][eqIdx]
                != BoundaryConditions::neumann)
            {
                // dirichlet boundary conditions are treated as
                // zeros on the right hand side
                this->b[scvIdx][eqIdx] = 0;
                continue;
            }

            // if we are here we've got a neumann boundary condition

            // make sure that we evaluate call the problem's
            // neumann() method exactly once if
            if (!wasEvaluated)
            {
                // make sure that we only evaluate
                // the neumann fluxes once
                wasEvaluated = true;

                problem_.neumann(values,
                                 curElement_(),
                                 curElementGeom_,
                                 *isIt,
                                 scvIdx,
                                 boundaryFaceIdx);
                Valgrind::CheckDefined(values);
                // TODO (?): multiple integration
                // points
                values *= curElementGeom_.boundaryFace[boundaryFaceIdx].area;
            }

            Valgrind::CheckDefined(values[eqIdx]);

            this->b[scvIdx][eqIdx] += values[eqIdx];
            Valgrind::CheckDefined(this->b[scvIdx][eqIdx]);
        }
    }

    /*!
     * \brief Compute the local residual, i.e. the right hand side
     *        of an equation we would like to have zero.
     */
    void evalLocalResidual(SolutionOnElement &residual,
                           bool withBoundary = true)
    {
        // reset residual
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int j = 0; j < numEq; ++j)
                residual[i][j] = Scalar(0);
        }
        
        if (withBoundary)
            this->asImp_().assembleBoundaryCondition(this->curElement_());

        this->evalFluxes_(residual);

        // evaluate the local rate
        for (int i=0; i < curElementGeom_.numVertices; i++)
        {
            PrimaryVarVector massContrib(0), tmp(0);

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            this->asImp_().computeStorage(massContrib, i, false);
            this->asImp_().computeStorage(tmp, i, true);

            massContrib -= tmp;
            massContrib *= curElementGeom_.subContVol[i].volume/problem_.timeStepSize();
            
            for (int j = 0; j < numEq; ++j)
                residual[i][j] += massContrib[j];

            // subtract the source term from the local rate
            PrimaryVarVector source;
            this->asImp_().computeSource(source, i);
            source *= curElementGeom_.subContVol[i].volume;

            for (int j = 0; j < numEq; ++j) {
                residual[i][j] -= source[j];
                
                // make sure that only defined quantities where used
                // to calculate the residual.
                Valgrind::CheckDefined(residual[i][j]);
            }
        }

        if (withBoundary) {
            for (int i = 0; i < curElementGeom_.numVertices; i++) {
                for (int j = 0; j < numEq; j++) {
                    if (this->bctype[i][j] == BoundaryConditions::dirichlet)
                        residual[i][j] = 0;
                    else
                        residual[i][j] += this->b[i][j];
                }
            }
        }

#if HAVE_VALGRIND
        for (int i=0; i < curElementGeom_.numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Set the current grid element.
     */
    void setCurrentElement(const Element &element)
    {
        if (curElementPtr_ != element) {
            // update the FV element geometry if the current element
            // is to be changed
            curElementPtr_ = element;
            curElementGeom_.update(curElement_());

            // tell LocalStiffness (-> parent class) the current
            // number of degrees of freedom
            this->setcurrentsize(curElementGeom_.numVertices);
	}
    };


    /*!
     * \brief Restrict the global function 'globalFn' to the vertices
     *        of the current element, save the result to 'dest'.
     */
    void restrictToElement(SolutionOnElement &dest,
                           const SolutionFunction &globalFn) const
    {
        const DofEntityMapper &dofMapper = this->problem_.model().dofEntityMapper();
        // we assert that the i-th shape function is
        // associated to the i-th vert of the element.
        int n = curElement_().template count<dim>();
        dest.resize(n);
        for (int i = 0; i < n; i++) {
            dest[i] = (*globalFn)[dofMapper.map(curElement_(), i, dim)];
        }
    }

    /*!
     * \brief Initialize the static data of all elements. The current
     *        solution is not yet valid when this method is
     *        called. (updateStaticData() is called as soon the
     *        current solution is valid.)
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void initStaticData()
    { }

    /*!
     * \brief Update the static data of all elements with the current solution
     *
     * This method should be overwritten by the child class if
     * necessary.
     */
    void updateStaticData(SolutionFunction &curSol, SolutionFunction &oldSol)
    { };

    /*!
     * \brief Update the model specific vertex data of a whole
     *        element.
     */
    void updateElementData_(VertexDataArray &dest, const SolutionOnElement &sol, bool isOldSol)
    {
        dest.resize(sol.size());

#ifdef ENABLE_VALGRIND
        for (int i = 0; i < dest.size(); ++i)
            Valgrind::SetUndefined(dest[i]);
#endif // ENABLE_VALGRIND

        int numVertices = this->curElement_().template count<dim>();
        for (int vertIdx = 0; vertIdx < numVertices; vertIdx++) {
            dest[vertIdx].update(sol[vertIdx],
                                 curElement_(),
                                 curFvElementGeometry(),
                                 vertIdx,
                                 problem(), 
                                 isOldSol);
        }
    }


    /*!
     * \brief This returns the finite volume geometric information of the current
     *        element. <b>Use this method with great care!</b>
     */
    const FVElementGeometry &curFvElementGeometry() const
    { return curElementGeom_; }

    /*!
     * \brief Set current local solution
     */
    void setCurrentSolution(const SolutionOnElement &sol)
    {
        curElemDat_.resize(sol.size());
        asImp_().updateElementData_(curElemDat_, sol, false);
    }

    /*!
     * \brief Set local solution of the last time step
     */
    void setPreviousSolution(const SolutionOnElement &sol)
    {
        prevElemDat_.resize(sol.size());
        asImp_().updateElementData_(prevElemDat_, sol, true);
    }

    /*!
     * \brief Vary a single component of a single vert of the
     *        local solution for the current element.
     *
     * This method is a optimization, since if varying a single
     * component at a degree of freedom not the whole element cache
     * needs to be recalculated. (Updating the element cache is very
     * expensive since material laws need to be evaluated.)
     */
    void deflectCurrentSolution(SolutionOnElement &curSol, 
                                int vertIdx,
                                int eqIdx,
                                Scalar value)
    {
        // stash away the orignial vertex data so that we do not need
        // to re calculate them if restoreCurSolution() is called.
        curVertexDataStash_ = curElemDat_[vertIdx];

        // recalculate the vertex data for the box which should be
        // changed
        curSol[vertIdx][eqIdx] = value;
        
        Valgrind::SetUndefined(curElemDat_[vertIdx]);
        curElemDat_[vertIdx].update(curSol[vertIdx], 
                                    curElement_(),
                                    curFvElementGeometry(),
                                    vertIdx, 
                                    problem(),
                                    false);
        Valgrind::CheckDefined(curElemDat_[vertIdx]);
    }

    /*!
     * \brief Restore the local jacobian to the state before
     *        deflectCurSolution() was called.
     *
     * This only works if deflectSolution was only called with
     * (vert, component) as arguments.
     */
    void restoreCurrentSolution(SolutionOnElement &curSol, 
                                int vertIdx, 
                                int eqIdx,
                                Scalar origValue)
    {
        curSol[vertIdx][eqIdx] = origValue;
        curElemDat_[vertIdx] = curVertexDataStash_;
    };

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem()
    { return problem_; };

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem() const
    { return problem_; };

    /*!
     * \brief Returns a reference to the model.
     */
    Model &model()
    { return problem_.model(); };

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model() const
    { return problem_.model(); };

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper() const
    { return model().vertexMapper(); };

    /*!
     * \brief Returns a reference to the element mapper.
     */
    const ElementMapper &elementMapper() const
    { return model().elementMapper(); };

protected:
    const Element &curElement_() const
    { return *curElementPtr_; }

    // The problem we would like to solve
    Problem &problem_;

    // The grid view we are dealing with
    const GridView gridView_;

    ElementPointer    curElementPtr_;
    FVElementGeometry curElementGeom_;

    // current and previous element data. (this is model specific.)
    VertexDataArray  curElemDat_;
    VertexDataArray  prevElemDat_;
    
    // temporary variable to store the variable vertex data
    VertexData   curVertexDataStash_;

protected:
    void updateBoundaryTypes_()
    {
        Dune::GeometryType      geoType = curElement_().geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        int numVerts = curElement_().template count<dim>();
        for (int i = 0; i < numVerts; ++i)
            for (int comp=0; comp < numEq; comp++)
                this->bctype[i][comp] = BoundaryConditions::neumann;

        // evaluate boundary conditions
        IntersectionIterator isIt = curElement_().ileafbegin();
        const IntersectionIterator &endIt = curElement_().ileafend();
        for (; isIt != endIt; ++isIt)
        {
            // Ignore non- boundary faces.
            if (!isIt->boundary())
                continue;

            // Set the bctype for all vertices of the face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 faceVertIdx++)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);
                int boundaryFaceIdx =
                    curElementGeom_.boundaryFaceIndex(faceIdx,
                                                      faceVertIdx);

                // set the boundary types
                BoundaryTypeVector tmp;
                tmp.reset();
                problem_.boundaryTypes(tmp,
                                       curElement_(),
                                       curFvElementGeometry(),
                                       *isIt,
                                       elemVertIdx,
                                       boundaryFaceIdx);
                tmp.checkWellPosed();
                Valgrind::CheckDefined(tmp);

                // copy boundary type to the bctype array.
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    // make sure that dirichlet boundaries have
                    // priority over neumann ones
                    if (this->bctype[elemVertIdx][eqIdx]
                        == BoundaryConditions::dirichlet)
                    {
                        continue;
                    }
                    
                    if (tmp.isDirichlet(eqIdx))
                        this->bctype[elemVertIdx][eqIdx] = BoundaryConditions::dirichlet;
                    else
                        this->bctype[elemVertIdx][eqIdx] = BoundaryConditions::neumann;
                    Valgrind::CheckDefined(this->bctype[elemVertIdx][eqIdx]);
                }
            }
        }
    };

    void assemble_(const Element &element, SolutionOnElement& localU)
    {
        // set the current grid element
        asImp_().setCurrentElement(element);

        // reset the right hand side and the bctype array
        resetRhs_();

        int numVertices = curElementGeom_.numVertices;

        // restrict the previous global solution to the current element
        SolutionOnElement localOldU(numVertices);
        restrictToElement(localOldU, problem_.model().prevSolFunction());

        this->asImp_().setCurrentSolution(localU);
        this->asImp_().setPreviousSolution(localOldU);

        // approximate the local stiffness matrix numerically
        // TODO: do this analytically if possible
        SolutionOnElement partialStiffness(numVertices);
        for (int j = 0; j < numVertices; j++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().assemblePartialStiffness_(partialStiffness, 
                                                   localU,
                                                   j,
                                                   pvIdx);
                
                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalStiffness_(j,
                                      pvIdx,
                                      partialStiffness);
            }
        }

#if !HAVE_DUNE_PDELAB
        // calculate the right hand side
        SolutionOnElement residU(numVertices);
        asImp_().evalLocalResidual(residU, true); // include boundary conditions for the residual

        for (int i=0; i < numVertices; i++) {
            for (int eqIdx=0; eqIdx < numEq; eqIdx++) {
                // TODO: in most cases this is not really a boundary
                // vertex but an interior vertex, so
                // Dune::BoundaryConditions::neumann is misleading...
                if (this->bctype[i][eqIdx] == Dune::BoundaryConditions::neumann)
                    this->b[i][eqIdx] = residU[i][eqIdx];
            }
        }
#endif
    };
    
protected:
    void evalFluxes_(SolutionOnElement &residual)
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int k = 0; k < curElementGeom_.numEdges; k++)
        {
            int i = curElementGeom_.subContVolFace[k].i;
            int j = curElementGeom_.subContVolFace[k].j;

            PrimaryVarVector flux;
            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, k);
            Valgrind::CheckDefined(flux);

            // subtract fluxes from the local mass rates of the
            // respective sub control volume adjacent to the face. We
            // ignore dirichlet cells because for them, the mass
            // change inside the cell is not equal to the flux out of
            // the cell.
            for (int eq = 0; eq < numEq; ++ eq) {
                if (this->bctype[i][eq] == BoundaryConditions::neumann) 
                    residual[i][eq] -= flux[eq];
                if (this->bctype[j][eq] == BoundaryConditions::neumann) 
                    residual[j][eq] += flux[eq];
            }
        }
    }


    /*!
     * \brief Update the stiffness matrix for all equations on all
     *        vertices of the current element with the partial
     *        derivative to the pvIdx's primary variable at the
     *        vertexIdx-th vertex.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than central differences ought to be used.
     */
    void assemblePartialStiffness_(SolutionOnElement &dest,
                                   SolutionOnElement &elemSol,
                                   int vertexIdx,
                                   int pvIdx)
    {
        Scalar eps = asImp_().numericEpsilon_(elemSol, vertexIdx, pvIdx);
        Scalar uJ = elemSol[vertexIdx][pvIdx];
        
        // vary the pvIdx-th primary variable at the element's j-th
        // vertex and calculate the residual, don't include the
        // boundary conditions
        asImp_().deflectCurrentSolution(elemSol, vertexIdx, pvIdx, uJ + eps);
        asImp_().evalLocalResidual(dest, false);
        asImp_().restoreCurrentSolution(elemSol, vertexIdx, pvIdx, uJ);

        asImp_().deflectCurrentSolution(elemSol, vertexIdx, pvIdx, uJ - eps);
        SolutionOnElement tmp(curElementGeom_.numVertices);
        asImp_().evalLocalResidual(tmp, false);
        asImp_().restoreCurrentSolution(elemSol, vertexIdx, pvIdx, uJ);
        
        // central differences
        dest -= tmp;
        dest /= 2*eps;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < dest.size(); ++i)
            Valgrind::CheckDefined(dest[i]);
#endif
    }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemSol    The current solution on the element
     * \param vertexIdx  The local index of the element's vertex for 
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon_(const SolutionOnElement &elemSol,
                           int vertIdx,
                           int pvIdx) const
    {
        return std::max(fabs(1e-7*elemSol[vertIdx][pvIdx]), 1e-7);
        //return std::max(std::abs(1e-9*elemSol[vertIdx][pvIdx]), Scalar(1e-12));
    }

    /*!
     * \brief Updates the current local stiffness matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at vertex j .
     */
    void updateLocalStiffness_(int j,
                               int pvIdx,
                               const SolutionOnElement &stiffness)
    {
        for (int i = 0; i < curElementGeom_.numVertices; i++) {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                // A[i][j][eqIdx][pvIdx] is the approximate rate of
                // change of the residual of equation 'eqIdx' at
                // vertex 'i' depending on the primary variable
                // 'pvIdx' at vertex 'j'.
                this->A[i][j][eqIdx][pvIdx] = stiffness[i][eqIdx];
                Valgrind::CheckDefined(this->A[i][j][eqIdx][pvIdx]);
            }
        }
    }

    // reset the right hand side
    void resetRhs_()
    {
        int numVertices = this->currentsize();
        for (int i=0; i < numVertices; i++) {
            for (int j = 0; j < numEq; ++j) {
                this->bctype[i][j] = Dune::BoundaryConditions::neumann;
                this->b[i][j] = 0;
            }
        }
    };

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
}

#endif
