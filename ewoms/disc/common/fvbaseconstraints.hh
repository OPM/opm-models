// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \copydoc Ewoms::FvBaseConstraints
 */
#ifndef EWOMS_FV_BASE_CONSTRAINTS_HH
#define EWOMS_FV_BASE_CONSTRAINTS_HH

#include <opm/core/utility/PropertySystem.hpp>
#include <opm/material/Valgrind.hpp>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
}}

namespace Ewoms {
/*!
 * \brief Class to specify constraints in a ECFV discretization.
 */
template <class TypeTag>
class FvBaseConstraints : public GET_PROP_TYPE(TypeTag, PrimaryVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    FvBaseConstraints()
    { reset(); }

//! \cond SKIP
    /*!
     * \brief Use the assignment operators from Dune::FieldVector
     */
    using ParentType::operator=;
//! \endcond

    /*!
     * \brief Reset the constraints types.
     *
     * If nothing is changed after this method has been called, no
     * DOFs will be constraint.
     */
    void reset()
    {
        for (int i=0; i < numEq; ++i) {
            constraintInfo_[i].isConstraint = 0;

            eq2pvIdx_[i] = i;
            pv2eqIdx_[i] = i;
        };
    }

    /*!
     * \brief Returns true if constraints for any equation have been
     *        specified.
     */
    bool isConstraint() const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (constraintInfo_[eqIdx].isConstraint)
                return true;
        return false;
    }

    /*!
     * \brief Returns true if constraints for given equation have been
     *        specified.
     *
     * \param eqIdx The index of the equation
     */
    bool isConstraint(int eqIdx) const
    { return constraintInfo_[eqIdx].isConstraint; }

    /*!
     * \brief Set all to be constraint
     *
     * (with equation index == primary variable index)
     */
    void setAllConstraint()
    {
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            constraintInfo_[eqIdx].isConstraint = 1;

            eq2pvIdx_[eqIdx] = eqIdx;
            pv2eqIdx_[eqIdx] = eqIdx;
        }
        //Valgrind::SetDefined(constraintInfo_[eqIdx]);
    }

    /*!
     * \brief Set a constraint for single equation
     *
     * \param eqIdx The index of the equation which should used to set
     *              the constraint
     * \param pvIdx The index of the primary variable which should be set.
     * \param value The value of the constraint DOF
     */
    void setConstraint(int eqIdx, int pvIdx, Scalar value)
    {
        setConstraint(eqIdx, pvIdx);
        ParentType::operator[](pvIdx) = value;

        //Valgrind::SetDefined(constraintInfo_);
    }

    /*!
     * \brief Set a constraint for single equation
     *
     * \param eqIdx The index of the equation which should used to set
     *              the constraint
     * \param pvIdx The index of the primary variable which should be set.
     */
    void setConstraint(int eqIdx, int pvIdx)
    {
        constraintInfo_[eqIdx].isConstraint = 1;

        // update the equation <-> primary variable mapping
        eq2pvIdx_[eqIdx] = pvIdx;
        pv2eqIdx_[pvIdx] = eqIdx;

        //Valgrind::SetDefined(constraintInfo_[eqIdx]);
    }

    /*!
     * \brief Returns the index of the equation which should be used
     *        to constrain the pvIdx's primary variable.
     *
     * \param pvIdx The index of the primary variable which is be set
     *              by the constraint.
     */
    unsigned pvToEqIndex(unsigned pvIdx) const
    { return pv2eqIdx_[pvIdx]; }

    /*!
     * \brief Returns the index of the primary variable for which a
     *        given equation should be used as a constraint.
     *
     * \param eqIdx The index of the equation which is used to set the
     *              constraint.
     */
    unsigned eqToPvIndex(unsigned eqIdx) const
    { return eq2pvIdx_[eqIdx]; }

private:
    // this is a bitfield structure!
    struct __attribute__((__packed__)) {
        unsigned char isConstraint : 1;
    } constraintInfo_[numEq];

    unsigned char eq2pvIdx_[numEq];
    unsigned char pv2eqIdx_[numEq];
};

} // namespace Ewoms

#endif
