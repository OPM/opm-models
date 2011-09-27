// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
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
 * \brief Represents all quantities which are relevant for an element
 */
#ifndef DUMUX_BOX_ELEMENT_CONTEXT_HH
#define DUMUX_BOX_ELEMENT_CONTEXT_HH

#include "boxproperties.hh"


namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class BoxElementContext
{
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    // the history size of the time discretization in number of steps
    enum { timeDiscHistorySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize) };

    struct ScvStore_ {
        VolumeVariables volVars[timeDiscHistorySize];
        const VolumeVariables *hint[timeDiscHistorySize];
    };
    typedef std::vector<ScvStore_> ScvVarsVector;
    typedef std::vector<FluxVariables> ScvfVarsVector;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit BoxElementContext(const Problem &problem)
        : gridView_(problem.gridView())
    {
        // remember the problem object
        problemPtr_ = &problem;
        modelPtr_ = &problem.model();
    }

    /*!
     * \brief Construct the volume variables of an element from scratch.
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     */
    void updateAll(const Element &elem)
    {
        updateFVElemGeom(elem);
        updateBoundaryTypes();
        updateAllScvVars();
        updateAllScvfVars();
    }

    void updateFVElemGeom(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the finite element geometry
        fvElemGeom_.update(gridView_, elem);

        // resize the SCV and the SCVF arrays
        if (fvElemGeom_.numVertices > scvVars_.size()) {
            scvVars_.resize(fvElemGeom_.numVertices);
            boundaryTypes_.resize(fvElemGeom_.numVertices);
        }
        if (fvElemGeom_.numEdges > scvfVars_.size())
            scvfVars_.resize(fvElemGeom_.numEdges);
    }

    void updateBoundaryTypes()
    {
        // fill the boundary types stuff
        hasNeumann_ = false;
        hasDirichlet_ = false;
        for (int scvIdx = 0; scvIdx < numScv(); ++scvIdx) {
            int globalIdx = globalIndex(scvIdx);
            boundaryTypes_[scvIdx] = model().boundaryTypes(globalIdx);
            hasDirichlet_ = hasDirichlet_ || boundaryTypes_[scvIdx]->hasDirichlet();
            hasNeumann_ = hasNeumann_ || boundaryTypes_[scvIdx]->hasNeumann();
        }

        // set the default evaluation points
        scvIdxSaved_ = -1;
        scvfVarsEval_ = &scvfVars_;
    }

    void updateAllScvVars()
    {
        for (int historyIdx = 0; historyIdx < timeDiscHistorySize; ++ historyIdx)
            updateScvVars(historyIdx);
    };

    void updateScvVars(int historyIdx)
    {
        // update the volume variables for the whole history
        const VertexMapper &vertexMapper = problem().vertexMapper();
        const SolutionVector &globalSol = model().solution(historyIdx);

        int nScv = numScv();
        for (int scvIdx = 0; scvIdx < nScv; scvIdx++) {
            int globalIdx = vertexMapper.map(element(), scvIdx, dim);
            const PrimaryVariables &scvSol = globalSol[globalIdx];

            scvVars_[scvIdx].hint[historyIdx] = model().hint(globalIdx, historyIdx);
            updateScvVars(scvSol, scvIdx, historyIdx);
        }
    }

    void updateScvVars(const PrimaryVariables &priVars, int scvIdx, int historyIdx)
    {
        scvVars_[scvIdx].volVars[historyIdx].update(priVars,
                                                    /*context=*/*this,
                                                    scvIdx,
                                                    historyIdx);
    }

    void updateAllScvfVars()
    {
        for (int scvfIdx = 0; scvfIdx < numScvf(); scvfIdx++) {
            scvfVars_[scvfIdx].update(/*context=*/ *this,
                                      /*localIndex=*/scvfIdx);
        }
    }

    /*!
     * \brief Construct the volume variables for all of vertices of an
     *        element given a solution vector computed by PDELab.
     *
     * \tparam ElemSolVectorType The container type which stores the
     *                           primary variables of the element
     *                           using _local_ indices
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     * \param fvElemGeom The finite volume geometry of the element
     * \param elementSolVector The local solution for the element using PDELab ordering
     */
    template<typename ElemSolVectorType>
    void updatePDELab(const Element &element,
                      const ElemSolVectorType &elementSolVector)
    {
        updateFVElemGeom(element);
        updateBoundaryTypes(element);

        // update the current time step's volume variables
        PrimaryVariables scvSol;
        for (int scvIdx = 0; scvIdx < numScv(); scvIdx++)
        {
            // reorder the solution
            for (int eqnIdx = 0; eqnIdx < numEq; eqnIdx++)
                scvSol[eqnIdx] = elementSolVector[scvIdx + eqnIdx*numScv()];

            // update the volume variables for the newest history index
            updateScvVars_(scvSol, /*historyIdx=*/0, scvIdx);
        }
    }

    /*!
     * \brief Return a reference to the problem.
     */
    const Problem &problem() const
    { return *problemPtr_; }

    /*!
     * \brief Return a reference to the model.
     */
    const Model &model() const
    { return *modelPtr_; }

    /*!
     * \brief Return a reference to the grid view.
     */
    const GridView &gridView() const
    { return gridView_; }

    /*!
     * \brief Return the current element.
     */
    const Element &element() const
    { return *elemPtr_; }

    /*!
     * \brief Return the number of sub-control volumes of the current element.
     */
    int numScv() const
    { return fvElemGeom_.numVertices; }

    /*!
     * \brief Return the number of sub-control volume faces of the current element.
     */
    int numScvf() const
    { return fvElemGeom_.numEdges; }

    /*!
     * \brief Return the current finite element geometry.
     */
    const FVElementGeometry &fvElemGeom() const
    { return fvElemGeom_; }

    /*!
     * \brief Return the position of a local entities in global coordinates
     */
    const GlobalPosition &pos(int scvIdx) const
    { return fvElemGeom_.subContVol[scvIdx].global; }

    /*!
     * \brief Return the global index for a sub-control volume
     */
    int globalIndex(int scvIdx) const
    { return model().vertexMapper().map(element(), scvIdx, dim); }

    /*!
     * \brief Returns whether the current element is on the domain's
     *        boundary.
     */
    bool onBoundary() const
    { return hasNeumann_ || hasDirichlet_; };

    /*!
     * \brief Returns whether the current element has a Neumann boundary segment.
     */
    bool hasNeumann() const
    { return hasNeumann_; };

    /*!
     * \brief Returns whether the current element has a Dirichlet vertex
     */
    bool hasDirichlet() const
    { return hasDirichlet_; };

    /*!
     * \brief Returns the boundary types for a given vertex
     */
    const BoundaryTypes &boundaryTypes(int scvIdx) const
    {
        return *boundaryTypes_[scvIdx];
    }

    /*!
     * \brief Save the current flux variables and use them as the
     *        evaluation point.
     */
    void saveScvfVars()
    {
        scvfVarsSaved_ = scvfVars_;

        // change evaluation point
        scvfVarsEval_ = &scvfVarsSaved_;
    }

    /*!
     * \brief Restore current flux variables from the saved ones.
     */
    void restoreScvfVars()
    {
        //scvfVarsSaved_ = scvfVars_; // not needed

        // change evaluation point
        scvfVarsEval_ = &scvfVars_;
    }

    /*!
     * \brief Return a reference to the volume variables of a
     *        sub-control volume at a given time.
     *
     * If the time step index is not given, return the volume
     * variables for the current time. 
     *
     * \param scvIdx The local index of the sub-control volume for
     *               which the volume variables are requested
     * \param historyIdx The index of the time step for which the 
     *                    volume variables are requested. 0 means 
     *                    current time step, 1 previous time step,
     *                    2 next-to-previous, etc.
     */
    const VolumeVariables &volVars(int scvIdx, int historyIdx = 0) const
    { return scvVars_[scvIdx].volVars[historyIdx]; }

    const VolumeVariables *hint(int scvIdx, int historyIdx = 0) const
    { return scvVars_[scvIdx].hint[historyIdx]; }
    /*!
     * \copydoc volVars()
     */
    VolumeVariables &volVars(int scvIdx, int historyIdx = 0)
    { return scvVars_[scvIdx].volVars[historyIdx]; }

    PrimaryVariables &primaryVars(int scvIdx, int historyIdx = 0)
    { return scvVars_[scvIdx].volVars[historyIdx].primaryVars(); }
    const PrimaryVariables &primaryVars(int scvIdx, int historyIdx = 0) const
    { return scvVars_[scvIdx].volVars[historyIdx].primaryVars(); }

    /*!
     * \brief Returns the volume variables at the evaluation point.
     */
    void saveScvVars(int scvIdx)
    { 
        scvIdxSaved_ = scvIdx;
        scvVarsSaved_ = scvVars_[scvIdx].volVars[/*historyIdx=*/0];
    }

    /*!
     * \brief Restores the volume variables at the evaluation point.
     */
    void restoreScvVars(int scvIdx)
    { 
        scvIdxSaved_ = -1;
        scvVars_[scvIdx].volVars[/*historyIdx=*/0] = scvVarsSaved_;
    }

    /*!
     * \brief Return a reference to the flux variables of a
     *        sub-control volume face.
     *
     * \param scvfIdx The local index of the sub-control volume face for
     *               which the flux variables are requested
     */
    const FluxVariables &fluxVars(int scvIdx) const
    { return scvfVars_[scvIdx]; }

    /*!
     * \brief Return a reference to the flux variables of a
     *        sub-control volume face for the evaluation point.
     *
     * \param scvfIdx The local index of the sub-control volume face for
     *               which the flux variables are requested
     */
    const FluxVariables &evalPointFluxVars(int scvfIdx) const
    { 
        return (*scvfVarsEval_)[scvfIdx];
    }

    /*!
     * \brief Returns the volume variables for history index 0 at the
     *        evaluation point.
     */
    const VolumeVariables &evalPointVolVars(int scvIdx) const
    { 
        if (scvIdxSaved_ == scvIdx)
            return scvVarsSaved_;
        return volVars(scvIdx, /*historyIdx=*/0);
    }

protected:
    ScvVarsVector scvVars_;  

    int scvIdxSaved_;
    VolumeVariables scvVarsSaved_;

    ScvfVarsVector scvfVars_;
    ScvfVarsVector scvfVarsSaved_;

    ScvfVarsVector *scvfVarsEval_;

    const Problem *problemPtr_;
    const Model *modelPtr_;
    const Element *elemPtr_;
    const GridView gridView_;
    FVElementGeometry fvElemGeom_;
    bool hasNeumann_;
    bool hasDirichlet_;
    std::vector<const BoundaryTypes*> boundaryTypes_;
};

} // namespace Dumux

#endif
