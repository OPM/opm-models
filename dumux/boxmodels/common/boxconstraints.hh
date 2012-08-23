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
 * \brief Class to specify constraints in a box model.
 */
#ifndef DUMUX_BOX_CONSTRAINTS_HH
#define DUMUX_BOX_CONSTRAINTS_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NumEq);
}

/*!
 * \brief Class to specify constraints in a box model.
 */
template <class TypeTag>
class BoxConstraints : public GET_PROP_TYPE(TypeTag, PrimaryVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) ParentType;   
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    
public:
    BoxConstraints()
    { reset(); }

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

            Valgrind::SetDefined(constraintInfo_[eqIdx]);
        }
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

        Valgrind::SetDefined(constraintInfo_[eqIdx]);
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

        Valgrind::SetDefined(constraintInfo_[eqIdx]);
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

}

#endif
