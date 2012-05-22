// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Base class for all problems which use the single-phase,
 *        two-component box model
 */
#ifndef DUMUX_1P2C_PROBLEM_HH
#define DUMUX_1P2C_PROBLEM_HH

#include <dumux/boxmodels/common/boxporousproblem.hh>
#include <dune/common/fvector.hh>
#include "1p2cproperties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxBaseProblems
 * \brief Base class for all problems which use the single-phase, two-component box model
 *
 */
template<class TypeTag>
class OnePTwoCBoxProblem : public BoxPorousProblem<TypeTag>
{
    typedef BoxPorousProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    OnePTwoCBoxProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          gravity_(0)
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
     * \brief Return the parameters for the material law.
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        static MaterialLawParams matParams;
        return matParams;
    }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * \param context Container for the volume variables, element,
     *                fvElementGeometry, etc
     * \param spaceIdx The local index of the sub control volume inside
     *                 the element
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

    /*!
     * \brief Specify whether two-point gradients should be used
     *        instead of finite-element ones.
     */
    template <class Context>
    bool useTwoPointGradient(const Context &context,
                             int spaceIdx,
                             int timeIdx) const
    { return false; }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    DimVector gravity_;
};

}

#endif
