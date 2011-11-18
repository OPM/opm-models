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
 * \ingroup BoxModels
 * \brief The base class for the problems of box models which assume a
 *        porous medium.
 */
#ifndef DUMUX_BOX_POROUS_PROBLEM_HH
#define DUMUX_BOX_POROUS_PROBLEM_HH

#include <dumux/boxmodels/common/boxproblem.hh>
#include <dumux/common/math.hh>

namespace Dumux {
/*!
 * \ingroup BoxModels
 */

/*!
 * \brief The base class for the problems of box models which assume a
 *        porous medium.
 */
template<class TypeTag>
class BoxPorousProblem : public BoxProblem<TypeTag>
{
    typedef Dumux::BoxProblem<TypeTag> ParentType;
    
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

    enum {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    BoxPorousProblem(TimeManager &timeManager, const GridView &gv)
        : ParentType(timeManager, gv)
    { }

    /*!
     * \brief Averages the intrinsic permeability (Scalar).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(Tensor &result,
               Scalar K1,
               Scalar K2) const
    {
        const Scalar K = Dumux::harmonicMean(K1, K2);
        for (int i = 0; i < dimWorld; ++i) {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Averages the intrinsic permeability (Tensor).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(Tensor &result,
               const Tensor &K1,
               const Tensor &K2) const
    {
        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = harmonicMean(K1[i][j], K2[i][j]);
    }

    template <class Context>
    const Tensor intrinsicPermeability(const Context &context,
                                       int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::intrinsicPermeability()");
    }

    template <class Context>
    Scalar porosity(const Context &context,
                    int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::porosity()");
    }
    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context,
                             int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::heatCapacitySolid()");
    }

    /*!
     * \brief Returns the thermal conductivity [W / (K m)] of the solid
     *        phase disregarding the pores in a sub-control volume.
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::heatConductionParams()");
    }

    /*!
     * \brief Define the tortuosity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar tortuosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::tortuosity()");
    }

    /*!
     * \brief Define the dispersivity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar dispersivity(const Context &context,
                        int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::dispersivity()");
    }


private:
    const Tensor &toTensor_(const Tensor &val) const
    { return val; };

    Tensor toTensor_(Scalar val) const
    {
        Tensor ret(0.0);
        for (int i = 0; i < Tensor::rows; ++i)
            ret[i][i] = val;
        return ret;
    };
};

} // namespace Dumux

#endif
