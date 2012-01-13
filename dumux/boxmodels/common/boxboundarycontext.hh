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
 * \brief Represents all quantities which available on boundary elements
 */
#ifndef DUMUX_BOX_BOUNDARY_CONTEXT_HH
#define DUMUX_BOX_BOUNDARY_CONTEXT_HH

#include "boxproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxModel
 *
 * \brief Represents all quantities which available on boundary elements
 */
template<class TypeTag>
class BoxBoundaryContext
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit BoxBoundaryContext(const Problem &problem)
        : gridView_(problem.gridView())
    {
        // remember the problem object
        problemPtr_ = &problem;
        modelPtr_ = &problem.model();
    }

    /*!
     * \brief Set the boundary context
     *
     * \param problem The problem which needs to be simulated.
     * \param element The DUNE Codim<0> entity for which the volume variables ought to be calculated
     */
    void update(const Element &elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        fvElemGeomUpToDate_ = false;
    };

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
    { return element().template count<dim>(); }

    /*!
     * \brief Return the number of sub-control volume faces of the current element.
     */
    int numScvf() const
    { return element().template count<dim - 1>(); }

    /*!
     * \brief Return the current finite element geometry.
     */
    const FVElementGeometry &fvElemGeom(int timeIdx) const
    {
        if (!fvElemGeomUpToDate_) {
            fvElemGeomUpToDate_ = true;
            fvElemGeom_.update(gridView_, element());
        }
        return fvElemGeom_;
    }

    /*!
     * \brief Return the position of a local entities in global coordinates
     */
    const GlobalPosition pos(int scvIdx, int timeIdx) const
    { return element().geometry().corner(scvIdx); }

protected:
    const Problem *problemPtr_;
    const Model *modelPtr_;
    const Element *elemPtr_;
    const GridView gridView_;

    bool fvElemGeomUpToDate_;
    FVElementGeometry fvElemGeom_;
};

} // namespace Dumux

#endif
