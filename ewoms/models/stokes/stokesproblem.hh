// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
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
 * \copydoc Ewoms::StokesProblem
 */
#ifndef EWOMS_STOKES_PROBLEM_HH
#define EWOMS_STOKES_PROBLEM_HH

#include "stokesproperties.hh"

#include <ewoms/disc/vcfv/vcfvproblem.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup VCFVStokesProblems
 * \brief Base class for all problems which use the Stokes VCVF discretization.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 */
template<class TypeTag>
class StokesProblem : public VcfvProblem<TypeTag>
{
    typedef VcfvProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \copydoc VcfvProblem::VcfvProblem(TimeManager &, const GridView &)
     */
    StokesProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , gravity_(0)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \brief Register all run-time parameters for the problem.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableGravity, "Use the gravity correction for the pressure gradients.");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature at a spatial and temporal
     *        position within the domain.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { OPM_THROW(std::logic_error, "Not implemented:Problem does not provide a temperature() method"); }

    /*!
     * \brief Returns the heat capacity [J/(K m^3)] of the solid phase
     *        with no pores in the sub-control volume.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context,
                             int spaceIdx, int timeIdx) const
    { return 0; }

    /*!
     * \brief Returns the parameter object for the heat conductivity law in
     *        a sub-control volume.
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        static const HeatConductionLawParams dummy;
        return dummy;
    }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else
     * \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    const DimVector &gravity(const Context &context, int spaceIdx, int timeIdx) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const DimVector &gravity() const
    { return gravity_; }

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

} // namespace Ewoms

#endif
