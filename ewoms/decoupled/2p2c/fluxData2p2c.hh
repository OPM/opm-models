// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
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
 * \copydoc Ewoms::FluxData2P2C
 */
#ifndef EWOMS_FLUXDATA2P2C_HH
#define EWOMS_FLUXDATA2P2C_HH

#include <ewoms/decoupled/common/decoupledproperties.hh>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

namespace Ewoms
{

/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FluxData2P2C
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };


    enum
    {
        numEquations = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typename Dune::BlockVector<typename Dune::FieldVector<bool, numEquations> > isUpwindCell_;

public:

    //! Constructor
    FluxData2P2C()
    {
        isUpwindCell_.resize(2 * dim);
        for (int i = 0; i<2*dim; i++)
        {
            isUpwindCell_[i] = false;
        }
    }
    //! resizes the upwind vector for the case of hanging nodes
    void resize(int size)
    {
        isUpwindCell_.resize(size);
    }
    //! returns the size of the upwind vector which equals number of faces
    int size()
    {
        return isUpwindCell_.size();
    }
    //! functions returning upwind information
    /* \param indexInInside The local inside index of the intersection
     * \param equationIdx The equation index
     */
    const bool& isUpwindCell(int indexInInside, int equationIdx) const
    {
        return isUpwindCell_[indexInInside][equationIdx];
    }
    //! Sets the upwind information
    /* \param indexInInside The local inside index of the intersection
     * \param equationIdx The equation index
     * \return value set true or false
     */
    void setUpwindCell(int indexInInside, int equationIdx, bool value)
    {
        isUpwindCell_[indexInInside][equationIdx] = value;
    }

    //! Console output for the FluxData
    void outputFluxData()
    {
        for(int banana=0; banana<isUpwindCell_.size(); banana++)
            printvector(std::cout, isUpwindCell_, "upwindInformation", "row", 3);
    }
};
}

#endif
