// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::VcfvVolumeVariables
 */
#ifndef EWOMS_VCFV_VOLUME_VARIABLES_HH
#define EWOMS_VCFV_VOLUME_VARIABLES_HH

#include "vcfvproperties.hh"

#include <opm/material/Valgrind.hpp>

namespace Ewoms {

/*!
 * \ingroup VcfvModel
 *
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
template <class TypeTag>
class VcfvVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    // default constructor
    VcfvVolumeVariables()
    { evalPoint_ = 0; }

    // copy constructor
    VcfvVolumeVariables(const VcfvVolumeVariables &v)
    {
        evalPoint_ = 0;
        extrusionFactor_ = v.extrusionFactor_;
    }

    /*!
     * \brief Assignment operator
     */
    VcfvVolumeVariables &operator=(const VcfvVolumeVariables &v)
    {
        evalPoint_ = 0;
        extrusionFactor_ = v.extrusionFactor_;

        return *this;
    }

    /*!
     * \brief Register all run-time parameters for the volume variables.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Sets the evaluation point used by the local jacobian.
     *
     * The evaluation point is only used by semi-smooth models.
     *
     * \param ep A pointer to the VolumeVariables object of the evaluation point
     */
    void setEvalPoint(const Implementation *ep)
    {
        evalPoint_ = ep;
        Valgrind::CheckDefined(evalPoint_);
    }

    /*!
     * \brief Returns the evaluation point used by the local jacobian.
     *
     * The evaluation point is only used by semi-smooth models.
     */
    const Implementation &evalPoint() const
    { return (evalPoint_ == 0)?asImp_():*evalPoint_; }

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param elemCtx The execution context from which the method is called.
     * \param scvIdx The index of the sub-control volume for which the
     *               volume variables should be calculated.
     * \param timeIdx The index for the time discretization for which
     *                the volume variables should be calculated
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        extrusionFactor_ = elemCtx.problem().extrusionFactor(elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \brief Update all gradients for a given control volume.
     *
     * \param elemCtx The execution context from which the method is called.
     * \param scvIdx The index of the sub-control volume for which the
     *               volume variables should be calculated.
     * \param timeIdx The index for the time discretization for which
     *                the volume variables should be calculated
     */
    void updateScvGradients(const ElementContext &elemCtx,
                            int scvIdx,
                            int timeIdx)
    { }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(evalPoint_);
        if (evalPoint_ && evalPoint_ != this)
            evalPoint_->checkDefined();
#endif
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    // the evaluation point of the local jacobian
    const Implementation *evalPoint_;

    Scalar extrusionFactor_;
};

} // end namepace

#endif
