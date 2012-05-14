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
 * \brief Implements a boundary vector for the fully implicit non-isothermal Stokes model.
 */
#ifndef DUMUX_STOKES_NI_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_STOKES_NI_BOUNDARY_RATE_VECTOR_HH

#include "stokesniproperties.hh"

#include <dumux/freeflow/stokes/stokesboundaryratevector.hh>
#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup StokesNIModel
 *
 * \brief Implements a boundary vector for the fully implicit non-isothermal Stokes model.
 */
template <class TypeTag>
class StokesNIBoundaryRateVector
    : public StokesBoundaryRateVector<TypeTag>
{
    typedef StokesBoundaryRateVector<TypeTag> ParentType;
    friend class StokesBoundaryRateVector<TypeTag>;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { energyEqIdx = Indices::energyEqIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    /*!
     * \brief Default constructor
     */
    StokesNIBoundaryRateVector()
        : ParentType()
    { };

    /*!
     * \brief Constructor with assignment from scalar
     */
    StokesNIBoundaryRateVector(Scalar value)
        : ParentType(value)
    { };

    /*!
     * \brief Copy constructor
     */
    StokesNIBoundaryRateVector(const StokesNIBoundaryRateVector &value)
        : ParentType(value)
    { };


protected:
    template <class FluidState>
    void enthalpyFlux_(const DimVector &velocity,
                       Scalar density,
                       const DimVector &normal,
                       const VolumeVariables &insideVolVars,
                       const FluidState &fs,
                       const typename FluidSystem::ParameterCache &paramCache)
    {
        Scalar vTimesN = normal * velocity;
        Scalar enthalpy =
            (vTimesN > 0)
            ? FluidSystem::enthalpy(fs, paramCache, phaseIdx)
            : insideVolVars.fluidState().enthalpy(phaseIdx);
       
        (*this)[energyEqIdx] = 
            vTimesN
            * enthalpy
            * density;
        Valgrind::CheckDefined((*this)[energyEqIdx]);
    }
};

} // end namepace

#endif
