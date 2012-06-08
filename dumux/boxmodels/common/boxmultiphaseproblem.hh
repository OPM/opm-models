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
 * \ingroup BoxModel
 *
 * \brief The base class for the problems of box models which deal
 *        with a multi-phase flow through a porous medium.
 */
#ifndef DUMUX_BOX_MULTI_PHASE_PROBLEM_HH
#define DUMUX_BOX_MULTI_PHASE_PROBLEM_HH

#include <dumux/boxmodels/common/boxporousproblem.hh>
#include <dumux/common/math.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(MaterialLawParams);
}

/*!
 * \ingroup BoxModel
 *
 * \brief The base class for the problems of box models which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class BoxMultiPhaseProblem : public BoxPorousProblem<TypeTag>
{
    typedef Dumux::BoxPorousProblem<TypeTag> ParentType;
    
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    BoxMultiPhaseProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dimWorld-1]  = -9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * \param context Container for the volume variables, element,
     *                fvElementGeometry, etc
     * \param spaceIdx The local index of the sub control volume inside
     *                 the element
     * \param timeIdx The index used by the time discretization.
     */     
    template <class Context>
    const MaterialLawParams& 
    materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::materialLawParams()");
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * \param context Container for the volume variables, element,
     *                fvElementGeometry, etc
     * \param spaceIdx The local index of the sub control volume inside
     *                 the element
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }


    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * \param context Container for the volume variables, element,
     *                fvElementGeometry, etc
     * \param spaceIdx The local index of the sub control volume inside
     *                 the element
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimVector &gravity(const Context &context,
                          int spaceIdx, int timeIdx) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>EnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const DimVector &gravity() const
    { return gravity_; }

    // \}

protected:
    DimVector gravity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Dumux

#endif
