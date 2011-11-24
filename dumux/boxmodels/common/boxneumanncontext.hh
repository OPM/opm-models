// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Represents all quantities which available on Neumann boundary segments
 */
#ifndef DUMUX_BOX_NEUMANN_CONTEXT_HH
#define DUMUX_BOX_NEUMANN_CONTEXT_HH

#include "boxproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief Represents all quantities which available on Neumann boundary segments
 */
template<class TypeTag>
class BoxNeumannContext
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit BoxNeumannContext(const ElementContext &elemCtx)
        : elemCtx_(elemCtx)
        , intersectionIt_(gridView().ibegin(element()))
    { }

    /*!
     * \brief Return a reference to the problem.
     */
    const Problem &problem() const
    { return elemCtx_.problem(); }

    /*!
     * \brief Return a reference to the model.
     */
    const Model &model() const
    { return elemCtx_.model(); }

    /*!
     * \brief Return a reference to the grid view.
     */
    const GridView &gridView() const
    { return elemCtx_.gridView(); }

    /*!
     * \brief Return the current element.
     */
    const Element &element() const
    { return elemCtx_.element(); }

    /*!
     * \brief Returns a reference to the element variables.
     */
    const ElementContext &elemCtx() const
    { return elemCtx_; }

    /*!
     * \brief Return the number of sub-control volumes of the current element.
     */
    int numScv() const
    { return elemCtx_.numScv(); }

    /*!
     * \brief Return the number of sub-control volume faces of the current element.
     */
    int numScvf() const
    { return elemCtx_.numScvf(); }

    /*!
     * \brief Return the current finite element geometry.
     */
    const FVElementGeometry &fvElemGeom(int timeIdx) const
    { return elemCtx_.fvElemGeom(timeIdx); };

    /*!
     * \brief Return the position of a local entities in global coordinates
     */
    const GlobalPosition &pos(int boundaryFaceIdx, int timeIdx) const
    { return fvElemGeom(timeIdx).boundaryFace[boundaryFaceIdx].ipGlobal; }

    /*!
     * \brief Return the local sub-control volume index of the interior of a neumann segment
     */
    short scvIdx(int boundaryFaceIdx, int timeIdx) const
    { return fvElemGeom(timeIdx).subContVolIndex(boundaryFaceIdx); }

    /*!
     * \brief Return the volume variables for the finite volume in the
     *        interiour of the current Neumann segment
     */
    const VolumeVariables &volVars(int boundaryFaceIdx, int timeIdx) const
    {
        short insideScvIdx = scvIdx(boundaryFaceIdx, timeIdx);
        return elemCtx_.volVars(insideScvIdx, timeIdx);
    }

    /*!
     * \brief Return the intersection for the neumann segment
     *
     * TODO/HACK: The intersection should take a local index as an
     * argument. since that's not supported efficiently by the DUNE
     * grid interface, we just ignore the index argument here!
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
