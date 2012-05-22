// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief Implements a boundary vector for the fully implicit three-phase, three-component model.
 */
#ifndef DUMUX_BOX_3P3CNI_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_BOX_3P3CNI_BOUNDARY_RATE_VECTOR_HH

#include "3p3cniproperties.hh"

#include <dumux/boxmodels/3p3c/3p3cboundaryratevector.hh>
#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup 3P3CNIModel
 *
 * \brief Implements a boundary vector for the fully implicit three-phase, three-component model.
 */
template <class TypeTag>
class ThreePThreeCNIBoundaryRateVector
    : public ThreePThreeCBoundaryRateVector<TypeTag>
{
    typedef ThreePThreeCBoundaryRateVector<TypeTag> ParentType;
    friend class ThreePThreeCBoundaryRateVector<TypeTag>;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { energyEqIdx = Indices::energyEqIdx };

public:
    /*!
     * \brief Default constructor
     */
    ThreePThreeCNIBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \brief Constructor with assignment from scalar
     */
    ThreePThreeCNIBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Copy constructor
     */
    ThreePThreeCNIBoundaryRateVector(const ThreePThreeCNIBoundaryRateVector &value)
        : ParentType(value)
    { }


protected:
    template <class FluidState>
    void enthalpyFlux_(const FluxVariables &fluxVars,
                       const VolumeVariables &insideVolVars,
                       const FluidState &fs,
                       const typename FluidSystem::ParameterCache &paramCache,
                       int phaseIdx,
                       Scalar density)
    {
        Scalar enthalpy =
            (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
            ? FluidSystem::enthalpy(fs, paramCache, phaseIdx)
            : insideVolVars.fluidState().enthalpy(phaseIdx);
        
        (*this)[energyEqIdx] += 
            fluxVars.filterVelocityNormal(phaseIdx)
            * enthalpy
            * density;
    }
};

} // end namepace

#endif
