/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \brief This method contains all callback classes for quantities
 *        that are required by some flux variables
 */
#ifndef EWOMS_QUANTITY_CALLBACKS_HH
#define EWOMS_QUANTITY_CALLBACKS_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief Callback class for temperature.
 */
template <class TypeTag>
class TemperatureCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    TemperatureCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Return the temperature given the index of a
     *        degree of freedom within an element context.
     *
     * In this context, we assume that thermal equilibrium applies,
     * i.e. that the temperature of all phases is equal.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().temperature(/*phaseIdx=*/0); }

private:
    const ElementContext& elemCtx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a phase pressure.
 */
template <class TypeTag>
class PressureCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    PressureCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the fluid phase for which the pressure
     *        should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Return the pressure of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().pressure(phaseIdx_); }

private:
    const ElementContext& elemCtx_;
    short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a phase pressure.
 */
template <class TypeTag, class FluidState>
class BoundaryPressureCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    BoundaryPressureCallback(const ElementContext& elemCtx, const FluidState& boundaryFs)
        : elemCtx_(elemCtx)
        , boundaryFs_(boundaryFs)
    {}

    /*!
     * \brief Set the index of the fluid phase for which the pressure
     *        should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Return the pressure of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().pressure(phaseIdx_); }

    Scalar boundaryValue() const
    { return boundaryFs_.pressure(phaseIdx_); }

private:
    const ElementContext& elemCtx_;
    const FluidState& boundaryFs_;
    short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the density of a phase.
 */
template <class TypeTag>
class DensityCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    DensityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the fluid phase for which the density
     *        should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Return the density of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().density(phaseIdx_); }

private:
    const ElementContext& elemCtx_;
    short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the molar density of a phase.
 */
template <class TypeTag>
class MolarDensityCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    MolarDensityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the fluid phase for which the molar
     *        density should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Return the molar density of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().molarDensity(phaseIdx_); }

private:
    const ElementContext& elemCtx_;
    short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the viscosity of a phase.
 */
template <class TypeTag>
class ViscosityCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    ViscosityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the fluid phase for which the viscosity
     *        should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Return the viscosity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().viscosity(phaseIdx_); }

private:
    const ElementContext& elemCtx_;
    short phaseIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the velocity of a phase at the center of a DOF.
 */
template <class TypeTag>
class VelocityCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    VelocityCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Return the velocity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    const DimVector &operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).velocityCenter(); }

private:
    const ElementContext& elemCtx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for the velocity of a phase at the center of a DOF.
 */
template <class TypeTag>
class VelocityComponentCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    VelocityComponentCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the component of the velocity
     *        which should be returned.
     */
    void setDimIndex(short dimIdx)
    { dimIdx_ = dimIdx; }

    /*!
     * \brief Return the velocity of a phase given the index of a
     *        degree of freedom within an element context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).velocityCenter()[dimIdx_]; }

private:
    const ElementContext& elemCtx_;
    short dimIdx_;
};

/*!
 * \ingroup Discretization
 *
 * \brief Callback class for a mole fraction of a component in a phase.
 */
template <class TypeTag>
class MoleFractionCallback
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    MoleFractionCallback(const ElementContext& elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \brief Set the index of the fluid phase for which a mole
     *        fraction should be returned.
     */
    void setPhaseIndex(short phaseIdx)
    { phaseIdx_ = phaseIdx; }

    /*!
     * \brief Set the index of the component for which the mole
     *        fraction should be returned.
     */
    void setComponentIndex(short compIdx)
    { compIdx_ = compIdx; }

    /*!
     * \brief Return the mole fraction of a component in a phase given
     *        the index of a degree of freedom within an element
     *        context.
     */
    Scalar operator()(int dofIdx) const
    { return elemCtx_.volVars(dofIdx, /*timeIdx=*/0).fluidState().moleFraction(phaseIdx_, compIdx_); }

private:
    const ElementContext& elemCtx_;
    short phaseIdx_;
    short compIdx_;
};

} // namespace Ewoms

#endif
