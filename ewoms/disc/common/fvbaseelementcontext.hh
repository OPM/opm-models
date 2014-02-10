/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \copydoc Ewoms::FvBaseElementContext
 */
#ifndef EWOMS_FV_BASE_ELEMENT_CONTEXT_HH
#define EWOMS_FV_BASE_ELEMENT_CONTEXT_HH

#include "fvbaseproperties.hh"

#include <dune/common/fvector.hh>

#include <vector>

namespace Ewoms {

/*!
 * \ingroup Discretization
 *
 * \brief This class stores an array of VolumeVariables objects, one
 *        volume variables object for each of the element's vertices
 */
template<class TypeTag>
class FvBaseElementContext
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    // the history size of the time discretization in number of steps
    enum { timeDiscHistorySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize) };

    struct DofStore_ {
        VolumeVariables volVars[timeDiscHistorySize];
        PrimaryVariables priVars[timeDiscHistorySize];
        const VolumeVariables *hint[timeDiscHistorySize];
    };
    typedef std::vector<DofStore_> DofVarsVector;
    typedef std::vector<FluxVariables> FluxVarsVector;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, GradientCalculator) GradientCalculator;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    static const int dim = GridView::dimension;
    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static const int requireScvCenterGradients =
        GET_PROP_VALUE(TypeTag, RequireScvCenterGradients);

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit FvBaseElementContext(const Problem &problem)
        : gridView_(problem.gridView())
        , stencil_(gridView_)
    {
        // remember the problem object
        problemPtr_ = &problem;
        modelPtr_ = &problem.model();
    }

    /*!
     * \brief Construct all volume and flux variables of an element
     *        from scratch.
     *
     * \param elem The DUNE Codim<0> entity for which the volume
     *             variables ought to be calculated
     */
    void updateAll(const Element &elem)
    {
        updateStencil(elem);
        updateAllVolVars();
        updateAllFluxVars();
    }

    /*!
     * \brief Compute the finite volume geometry for an element.
     *
     * \param elem The grid element for which the finite volume
     *             geometry ought to be computed.
     */
    void updateStencil(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the stencil. the center gradients are kind of
        // expensive to calculate and most models don't need them, so
        // that we only do this if the model needs them
        stencil_.update(elem);
        if (requireScvCenterGradients)
            stencil_.updateCenterGradients();

        // resize the arrays containing the flux and the volume
        // variables
        volVars_.resize(stencil_.numDof());
        fluxVars_.resize(stencil_.numInteriorFaces());
    }

    /*!
     * \brief Update the topological part of the stencil, but nothing else.
     *
     * \param elem The grid element for which the finite volume
     *             geometry ought to be computed.
     */
    void updateStencilTopology(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the finite element geometry
        stencil_.updateTopology(elem);
    }

    /*!
     * \brief Compute the volume variables of all sub-control volumes
     *        of the current element for all time indices.
     */
    void updateAllVolVars()
    {
        for (int timeIdx = 0; timeIdx < timeDiscHistorySize; ++ timeIdx)
            updateVolVars(timeIdx);
        dofIdxSaved_ = -1;
    }

    /*!
     * \brief Compute the volume variables of all sub-control volumes
     *        of the current element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateVolVars(int timeIdx)
    {
        // update the volume variables for the whole history
        const SolutionVector &globalSol = model().solution(timeIdx);

        // update the non-gradient quantities
        int nDof = numDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < nDof; dofIdx++) {
            int globalIdx = globalSpaceIndex(dofIdx, timeIdx);
            const PrimaryVariables &volSol = globalSol[globalIdx];

            volVars_[dofIdx].hint[timeIdx] = model().hint(globalIdx, timeIdx);
            updateSingleVolVars_(volSol, dofIdx, timeIdx);
        }

        // update gradients
        for (int dofIdx = 0; dofIdx < nDof; dofIdx++) {
            volVars_[dofIdx].volVars[timeIdx].updateScvGradients(/*context=*/*this,
                                                                 dofIdx,
                                                                 timeIdx);
        }
    }

    /*!
     * \brief Compute the volume variables of a single sub-control
     *        volume of the current element for a single time index.
     *
     * \param priVars The PrimaryVariables which should be used to
     *                calculate the volume variables.
     * \param dofIdx The local index in the current element of the
     *               sub-control volume which should be updated.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateVolVars(const PrimaryVariables &priVars, int dofIdx, int timeIdx)
    {
        updateSingleVolVars_(priVars, dofIdx, timeIdx);

        // update gradients inside a sub control volume
        int nDof = numDof(/*timeIdx=*/0);
        for (int dofIdx = 0; dofIdx < nDof; dofIdx++) {
            volVars_[dofIdx].volVars[timeIdx].updateScvGradients(/*context=*/*this,
                                                                 dofIdx,
                                                                 timeIdx);
        }
    }

    /*!
     * \brief Compute the flux variables of all sub-control volume
     *        faces of the current element for all time indices.
     */
    void updateAllFluxVars()
    {
        updateFluxVars(/*timeIdx=*/0);
    }

    /*!
     * \brief Compute the flux variables of all sub-control volume
     *        faces of the current element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateFluxVars(int timeIdx)
    {
        fluxVarsEval_ = &fluxVars_;

        gradientCalculator_.prepare(/*context=*/*this, timeIdx);

        for (int fluxIdx = 0; fluxIdx < numInteriorFaces(timeIdx); fluxIdx++) {
            fluxVars_[fluxIdx].update(/*context=*/ *this,
                                      /*localIndex=*/fluxIdx,
                                      timeIdx);
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
    int numDof(int timeIdx) const
    { return stencil(timeIdx).numDof(); }

    /*!
     * \brief Return the number of primary degrees of freedom of the current element.
     */
    int numPrimaryDof(int timeIdx) const
    { return stencil(timeIdx).numPrimaryDof(); }

    /*!
     * \brief Return the number of non-boundary faces which need to be
     *        considered for the flux apporixmation.
     */
    int numInteriorFaces(int timeIdx) const
    { return stencil(timeIdx).numInteriorFaces(); }

    /*!
     * \brief Return the number of boundary faces which need to be
     *        considered for the flux apporixmation.
     */
    int numBoundaryFaces(int timeIdx) const
    { return stencil(timeIdx).numBoundaryFaces(); }

    /*!
     * \brief Return the current finite element geometry.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const Stencil &stencil(int timeIdx) const
    { return stencil_; }

    /*!
     * \brief Return the position of a local entities in global coordinates
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const GlobalPosition &pos(int dofIdx, int timeIdx) const
    { return stencil_.subControlVolume(dofIdx).globalPos(); }

    /*!
     * \brief Return the global spatial index for a sub-control volume
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    int globalSpaceIndex(int dofIdx, int timeIdx) const
    { return stencil(timeIdx).globalSpaceIndex(dofIdx); }

    /*!
     * \brief Returns whether the current element is on the domain's
     *        boundary.
     */
    bool onBoundary() const
    { return element().hasBoundaryIntersections(); }

    /*!
     * \brief Save the current flux variables and use them as the
     *        evaluation point.
     */
    void saveFluxVars()
    {
        fluxVarsSaved_ = fluxVars_;

        // change evaluation point
        fluxVarsEval_ = &fluxVarsSaved_;
    }

    /*!
     * \brief Restore current flux variables from the saved ones.
     */
    void restoreFluxVars()
    {
        //fluxVarsSaved_ = fluxVars_; // not needed

        // change evaluation point
        fluxVarsEval_ = &fluxVars_;
    }

    /*!
     * \brief Return a reference to the volume variables of a
     *        sub-control volume at a given time.
     *
     * If the time step index is not given, return the volume
     * variables for the current time.
     *
     * \param dofIdx The local index of the sub-control volume for
     *               which the volume variables are requested
     * \param timeIdx The index of the time step for which the
     *                    volume variables are requested. 0 means
     *                    current time step, 1 previous time step,
     *                    2 next-to-previous, etc.
     */
    const VolumeVariables &volVars(int dofIdx, int timeIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));
        return volVars_[dofIdx].volVars[timeIdx];
    }

    /*!
     * \brief Return the hint for a given local index.
     *
     * \sa Discretization::hint(int, int)
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const VolumeVariables *hint(int dofIdx, int timeIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));
        return volVars_[dofIdx].hint[timeIdx];
    }
    /*!
     * \copydoc volVars()
     */
    VolumeVariables &volVars(int dofIdx, int timeIdx)
    {
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));
        return volVars_[dofIdx].volVars[timeIdx];
    }

    /*!
     * \brief Return the primary variables for a given local index.
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    PrimaryVariables &primaryVars(int dofIdx, int timeIdx)
    {
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));
        return volVars_[dofIdx].priVars[timeIdx];
    }
    /*!
     * \copydoc primaryVars()
     */
    const PrimaryVariables &primaryVars(int dofIdx, int timeIdx) const
    {
        assert(0 <= dofIdx && dofIdx < numDof(timeIdx));
        return volVars_[dofIdx].priVars[timeIdx];
    }

    /*!
     * \brief Returns the volume variables at the evaluation point.
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     */
    void saveVolVars(int dofIdx)
    {
        assert(0 <= dofIdx && dofIdx < numDof(/*timeIdx=*/0));

        dofIdxSaved_ = dofIdx;
        volVarsSaved_ = volVars_[dofIdx].volVars[/*timeIdx=*/0];
        priVarsSaved_ = volVars_[dofIdx].priVars[/*timeIdx=*/0];
    }

    /*!
     * \brief Restores the volume variables at the evaluation point.
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     */
    void restoreVolVars(int dofIdx)
    {
        dofIdxSaved_ = -1;
        volVars_[dofIdx].priVars[/*timeIdx=*/0] = priVarsSaved_;
        volVars_[dofIdx].volVars[/*timeIdx=*/0] = volVarsSaved_;
    }

    /*!
     * \brief Return a reference to the gradient calculation class of
     *        the chosen spatial discretization.
     */
    const GradientCalculator& gradientCalculator() const
    { return gradientCalculator_; }

    /*!
     * \brief Return a reference to the flux variables of a
     *        sub-control volume face.
     *
     * \param fluxIdx The local index of the sub-control volume face for
     *               which the flux variables are requested
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const FluxVariables &fluxVars(int fluxIdx, int timeIdx) const
    { return fluxVars_[fluxIdx]; }

    /*!
     * \brief Return a reference to the flux variables of a
     *        sub-control volume face for the evaluation point.
     *
     * \param fluxIdx The local index of the sub-control volume face for
     *               which the flux variables are requested
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const FluxVariables &evalPointFluxVars(int fluxIdx, int timeIdx) const
    {
        if (timeIdx != 0)
            return fluxVars(fluxIdx, timeIdx);
        return (*fluxVarsEval_)[fluxIdx];
    }

    /*!
     * \brief Returns the volume variables for history index 0 at the
     *        evaluation point.
     *
     * \param dofIdx The local index of the sub-control-volume index
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const VolumeVariables &evalPointVolVars(int dofIdx, int timeIdx) const
    {
        if (timeIdx != 0)
            return volVars(dofIdx, timeIdx);
        if (dofIdxSaved_ == dofIdx)
            return volVarsSaved_;
        return volVars(dofIdx, /*timeIdx=*/0);
    }

protected:
    void updateSingleVolVars_(const PrimaryVariables &priVars, int dofIdx, int timeIdx)
    {
        volVars_[dofIdx].priVars[timeIdx] = priVars;
        volVars_[dofIdx].volVars[timeIdx].update(/*context=*/*this, dofIdx, timeIdx);
    }

    DofVarsVector volVars_;

    int dofIdxSaved_;
    VolumeVariables volVarsSaved_;
    PrimaryVariables priVarsSaved_;

    GradientCalculator gradientCalculator_;

    FluxVarsVector fluxVars_;
    FluxVarsVector fluxVarsSaved_;

    FluxVarsVector *fluxVarsEval_;

    const Problem *problemPtr_;
    const Model *modelPtr_;
    const Element *elemPtr_;
    const GridView gridView_;
    Stencil stencil_;
};

} // namespace Ewoms

#endif
