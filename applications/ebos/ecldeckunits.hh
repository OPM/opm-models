// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
 * \copydoc Ewoms::EclDeckUnits
 */
#ifndef EWOMS_ECL_DECK_UNITS_HH
#define EWOMS_ECL_DECK_UNITS_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/simulator.hh>

#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/Units/Dimension.hpp>

#include <vector>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Scalar);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class converts quantities between the SI units used by eWoms and the unit
 *         system used by the ECL deck file.
 *
 * This is a simpler and "typo safe" alternative to the unit system class of opm-parser.
 */
template <class TypeTag>
class EclDeckUnits
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    enum Dimension {
        temperature,
        pressure,
        saturation,
        length,
        time,
        permeability,
        transmissibility,
        gasDissolutionFactor,
        oilDissolutionFactor,
        liquidSurfaceVolume,
        gasSurfaceVolume,
        reservoirVolume,
        density,
        viscosity,

        liquidRate,
        gasRate,

        gasOilRatio,

        numDimensions
    };

    EclDeckUnits(const Simulator& simulator)
    {
        const auto eclDeck = simulator.gridManager().deck();
        const auto deckUnitSystem = eclDeck->getActiveUnitSystem();

        deckToSiFactor_.resize(numDimensions);
        deckToSiOffset_.resize(numDimensions, 0.0);

        deckToSiFactor_[temperature] = deckUnitSystem->getDimension("Temperature")->getSIScaling();
        deckToSiOffset_[temperature] = deckUnitSystem->getDimension("Temperature")->getSIOffset();

        deckToSiFactor_[pressure] = deckUnitSystem->getDimension("Pressure")->getSIScaling();;
        deckToSiFactor_[saturation] = 1.0;
        deckToSiFactor_[length] = deckUnitSystem->getDimension("Length")->getSIScaling();;
        deckToSiFactor_[time] = deckUnitSystem->getDimension("Time")->getSIScaling();;
        deckToSiFactor_[permeability] = deckUnitSystem->getDimension("Permeability")->getSIScaling();;
        deckToSiFactor_[transmissibility] = deckUnitSystem->getDimension("Transmissibility")->getSIScaling();;
        deckToSiFactor_[gasDissolutionFactor] = deckUnitSystem->getDimension("GasDissolutionFactor")->getSIScaling();;
        deckToSiFactor_[oilDissolutionFactor] = deckUnitSystem->getDimension("OilDissolutionFactor")->getSIScaling();;
        deckToSiFactor_[liquidSurfaceVolume] = deckUnitSystem->getDimension("LiquidSurfaceVolume")->getSIScaling();;
        deckToSiFactor_[gasSurfaceVolume] = deckUnitSystem->getDimension("GasSurfaceVolume")->getSIScaling();;
        deckToSiFactor_[reservoirVolume] = deckUnitSystem->getDimension("ReservoirVolume")->getSIScaling();;
        deckToSiFactor_[density] = deckUnitSystem->getDimension("Density")->getSIScaling();;
        deckToSiFactor_[viscosity] = deckUnitSystem->getDimension("Viscosity")->getSIScaling();;

        deckToSiFactor_[liquidRate] = deckToSiFactor_[liquidSurfaceVolume] / deckToSiFactor_[time];
        deckToSiFactor_[gasRate] = deckToSiFactor_[gasSurfaceVolume] / deckToSiFactor_[time];
        deckToSiFactor_[gasOilRatio] = deckToSiFactor_[gasSurfaceVolume] / deckToSiFactor_[liquidSurfaceVolume];
    }

    /*!
     * \brief Convert a single scalar value from deck to SI units.
     */
    Scalar deckToSi(Scalar value, Dimension dimens) const
    { return value*deckToSiFactor_[dimens] + deckToSiOffset_[dimens]; }

    /*!
     * \brief Convert an array of scalar quanities from deck to SI units.
     */
    template <class DeckScalar>
    void deckToSi(std::vector<DeckScalar>& values, Dimension dimens) const
    {
        for (unsigned i = 0; i < values.size(); ++i)
            values[i] = values[i]*deckToSiFactor_[dimens] + deckToSiOffset_[dimens];
    }

    /*!
     * \brief Convert a single scalar value from deck to SI units.
     */
    Scalar siToDeck(Scalar value, Dimension dimens) const
    { return (value - deckToSiOffset_[dimens])/deckToSiFactor_[dimens]; }

    /*!
     * \brief Convert an array of scalar quanities from SI to deck units.
     */
    template <class DeckScalar>
    void siToDeck(std::vector<DeckScalar>& values, Dimension dimens) const
    {
        for (unsigned i = 0; i < values.size(); ++i)
            values[i] = (values[i] - deckToSiOffset_[dimens])/deckToSiFactor_[dimens];
    }

private:
    std::vector<Scalar> deckToSiFactor_;
    std::vector<Scalar> deckToSiOffset_;
};

} // namespace Ewoms

#endif
