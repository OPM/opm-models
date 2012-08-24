// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \ingroup BoxModel
 *
 * \brief Represents all quantities which available on boundary segments
 */
#ifndef DUMUX_BOX_BOUNDARY_CONTEXT_HH
#define DUMUX_BOX_BOUNDARY_CONTEXT_HH

#include "boxproperties.hh"
#include "boxelementcontext.hh"

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief Represents all quantities which available on boundary segments
 */
template<class TypeTag>
class BoxBoundaryContext
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    enum { dim = GridView::dimension };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit BoxBoundaryContext(const ElementContext &elemCtx)
        : elemCtx_(elemCtx)
        , intersectionIt_(gridView().ibegin(element()))
    { }

    /*!
     * \copydoc Dumux::BoxElementContext::problem()
     */
    const Problem &problem() const
    { return elemCtx_.problem(); }

    /*!
     * \copydoc Dumux::BoxElementContext::model()
     */
    const Model &model() const
    { return elemCtx_.model(); }

    /*!
     * \copydoc Dumux::BoxElementContext::gridView()
     */
    const GridView &gridView() const
    { return elemCtx_.gridView(); }

    /*!
     * \copydoc Dumux::BoxElementContext::element()
     */
    const Element &element() const
    { return elemCtx_.element(); }

    /*!
     * \brief Returns a reference to the element context object.
     */
    const ElementContext &elemContext() const
    { return elemCtx_; }

    /*!
     * \copydoc Dumux::BoxElementContext::numScv()
     */
    int numScv() const
    { return elemCtx_.numScv(); }

    /*!
     * \copydoc Dumux::BoxElementContext::numScvf()
     */
    int numScvf() const
    { return elemCtx_.numScvf(); }

    /*!
     * \copydoc Dumux::BoxElementContext::fvElemGeom()
     */
    const FVElementGeometry &fvElemGeom(int timeIdx) const
    { return elemCtx_.fvElemGeom(timeIdx); }

    /*!
     * \brief Return the position of a local entity in global coordinates.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const GlobalPosition &pos(int boundaryFaceIdx, int timeIdx) const
    { return fvElemGeom(timeIdx).boundaryFace[boundaryFaceIdx].ipGlobal; }

    /*!
     * \brief Return the position of a control volume's center in global coordinates.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const GlobalPosition &cvCenter(int boundaryFaceIdx, int timeIdx) const
    { return fvElemGeom(timeIdx).subContVol[insideScvIndex(boundaryFaceIdx, timeIdx)].global; }

    /*!
     * \brief Return the local sub-control volume index of the
     *        interior of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    short insideScvIndex(int boundaryFaceIdx, int timeIdx) const
    { return fvElemGeom(timeIdx).subContVolIndex(boundaryFaceIdx); }

    /*!
     * \brief Return the global space index of the sub-control volume
     *        at the interior of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    int globalSpaceIndex(int boundaryFaceIdx, int timeIdx) const
    { return elemCtx_.globalSpaceIndex(insideScvIndex(boundaryFaceIdx, timeIdx), timeIdx); }

    /*!
     * \brief Return the volume variables for the finite volume in the
     *        interiour of a boundary segment
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const VolumeVariables &volVars(int boundaryFaceIdx, int timeIdx) const
    {
        short insideScvIdx = this->insideScvIndex(boundaryFaceIdx, timeIdx);
        return elemCtx_.volVars(insideScvIdx, timeIdx);
    }

    /*!
     * \brief Return the flux variables for a given boundary face.
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     * \param timeIdx The index of the solution used by the time discretization
     */
    const FluxVariables &fluxVars(int boundaryFaceIdx, int timeIdx) const
    { return elemCtx_.boundaryFluxVars(boundaryFaceIdx, timeIdx); }

    /*!
     * \brief Return the intersection for the neumann segment
     *
     * TODO/HACK: The intersection should take a local index as an
     * argument. since that's not supported efficiently by the DUNE
     * grid interface, we just ignore the index argument here!
     *
     * \param boundaryFaceIdx The local index of the boundary segment
     */
    const Intersection &intersection(int boundaryFaceIdx) const
    { return *intersectionIt_; }

    /*!
     * \brief Return the intersection for the neumann segment
     *
     * TODO/HACK: the intersection iterator can basically be
     * considered as an index which is manipulated externally, but
     * context classes should not store any indices. it is done this
     * way for performance reasons
     */
    IntersectionIterator &intersectionIt()
    { return intersectionIt_; }

protected:
    const ElementContext &elemCtx_;
    IntersectionIterator intersectionIt_;
};

} // namespace Dumux

#endif
