/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Gives the phase state.
 */
#ifndef DUMUX_2P_PHASE_STATE_HH
#define DUMUX_2P_PHASE_STATE_HH

#include <dumux/new_material/fluidstate.hh>
#include <dumux/boxmodels/2p/2pproperties.hh>

namespace Dune
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        2p model.
 */
template <class TypeTag>
class TwoPFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                           TwoPFluidState<TypeTag> >
{
    typedef TwoPFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        wPhaseIdx = Indices::wPhase,
        nPhaseIdx = Indices::nPhase,
    };

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    
public:
    void update (Scalar pressW, Scalar pressN, Scalar temperature)
    {
        phasePressure_[wPhaseIdx] = pressW;
        phasePressure_[nPhaseIdx] = pressN;
        temperature_=temperature;
    }

    void update (Scalar temperature)
    {
        phasePressure_[wPhaseIdx] = 1e5;
        phasePressure_[nPhaseIdx] = 1e5;
        temperature_ = temperature;
    }

    /*!
     * \brief Returns the pressure of a fluid phase [Pa].
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids [K].
     */
    Scalar temperature() const
    { return temperature_; };

private:
    Scalar phasePressure_[numPhases];
    Scalar temperature_;
};

} // end namepace

#endif
