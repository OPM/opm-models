// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::VcfvMultiPhaseProblem
 */
#ifndef EWOMS_VCFV_MULTI_PHASE_PROBLEM_HH
#define EWOMS_VCFV_MULTI_PHASE_PROBLEM_HH

#include <ewoms/disc/vcfv/vcfvproblem.hh>
#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <ewoms/common/math.hh>

#include <opm/material/fluidmatrixinteractions/mp/NullMaterialLaw.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(HeatConductionLawParams);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(VelocityModule);
}

/*!
 * \ingroup VcfvModel
 *
 * \brief The base class for the problems of VCVF discretizations which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class VcfvMultiPhaseProblem
    : public VcfvProblem<TypeTag>
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityBaseProblem
{
    typedef Ewoms::VcfvProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;
    typedef typename Opm::NullMaterialLaw<GET_PROP_VALUE(TypeTag, NumPhases),
                                          typename GET_PROP_TYPE(TypeTag, Scalar)>::Params MaterialLawParams;
    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc VcfvProblem::VcfvProblem(TimeManager &, const GridView &)
     */
    VcfvMultiPhaseProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    { init_(); }

    /*!
     * \brief Register all run-time parameters for the problem and the model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableGravity, "Use the gravity correction for the pressure gradients.");
    };

    /*!
     * \brief Averages the intrinsic permeability Tensor.
     *
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first node
     * \param K2 intrinsic permeability of the second node
     */
    void meanK(DimMatrix &result,
               const DimMatrix &K1,
               const DimMatrix &K2) const
    {
        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = harmonicMean(K1[i][j], K2[i][j]);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$ at a given position
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context,
                                           int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::intrinsicPermeability()");
    }

    /*!
     * \brief Returns the porosity [] of the porous medium for a given
     *        control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
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
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context,
                             int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::heatCapacitySolid()");
    }

    /*!
     * \brief Returns the parameter object for the heat conductivity law in
     *        a sub-control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::heatConductionParams()");
    }

    /*!
     * \brief Define the tortuosity.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar tortuosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::tortuosity()");
    }

    /*!
     * \brief Define the dispersivity.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    Scalar dispersivity(const Context &context,
                        int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::dispersivity()");
    }

    /*!
     * \brief Returns the material law parameters \f$\mathrm{[K]}\f$ within a control volume.
     *
     * If you get a compiler error at this method, you set the
     * MaterialLaw property to something different than
     * Opm::NullMaterialLaw. In this case, you have to overload the
     * matererialLaw() method in the derived class!
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */
    template <class Context>
    const MaterialLawParams &
    materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        static MaterialLawParams dummy;
        return dummy;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
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
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
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
    /*!
     * \brief Converts a Scalar value to an isotropic Tensor
     *
     * This is convenient e.g. for specifying intrinsic permebilities:
     * \code{.cpp}
     * auto permTensor = this->toDimMatrix_(1e-12);
     * \endcode
     *
     * \param val The scalar value which should be expressed as a tensor
     */
    DimMatrix toDimMatrix_(Scalar val) const
    {
        DimMatrix ret(0.0);
        for (int i = 0; i < DimMatrix::rows; ++i)
            ret[i][i] = val;
        return ret;
    }

    DimVector gravity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void init_()
    {
        gravity_ = 0.0;
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dimWorld-1]  = -9.81;
    }
};

} // namespace Ewoms

#endif
