/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Represents the primary variables used in the 2-phase box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_2P_PRIMARY_VARIABLES_HH
#define DUMUX_2P_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/immiscibleflash.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

namespace Dumux
{
/*!
 * \ingroup 2PModel
 *
 * \brief Represents the primary variables used in the M-phase,
 *        N-component box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class TwoPPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

    // primary variable indices
    enum { pressureIdx = Indices::pressureIdx };
    enum { saturationIdx = Indices::saturationIdx };

    enum {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        
        pwSn = Indices::pwSn,
        pnSw = Indices::pnSw,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) EnergyModule;

    typedef Dumux::ImmiscibleFlash<Scalar, FluidSystem> ImmiscibleFlash;

public:
    /*!
     * \brief Default constructor
     */
    TwoPPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
    };

    /*!
     * \brief Constructor with assignment from scalar
     */
    TwoPPrimaryVariables(Scalar value)
        : ParentType(value)
    { };

    /*!
     * \brief Copy constructor
     */
    TwoPPrimaryVariables(const TwoPPrimaryVariables &value)
        : ParentType(value)
    { };

    /*!
     * \brief Set the primary variables with the wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pw The pressure of the wetting phase [Pa]
     * \param Sw The saturation of the wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    template <class MaterialLaw>
    void assignImmiscibleFromWetting(Scalar T,
                                     Scalar pw,
                                     Scalar Sw,
                                     const typename MaterialLaw::Params &matParams)
    {
        ImmiscibleFluidState<Scalar, FluidSystem> fs;
        
        fs.setTemperature(T);
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1 - Sw);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(wPhaseIdx, pw);
        fs.setPressure(nPhaseIdx, pw + (pC[nPhaseIdx] - pC[wPhaseIdx]));
        
        assignNaive(fs);
    };

    /*!
     * \brief Set the primary variables with the non-wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pn The pressure of the non-wetting phase [Pa]
     * \param Sn The saturation of the non-wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    template <class MaterialLaw>
    void assignImmiscibleFromNonWetting(Scalar T,
                                        Scalar pn,
                                        Scalar Sn,
                                        const typename MaterialLaw::Params &matParams)
    {
        ImmiscibleFluidState<Scalar, FluidSystem> fs;
        
        fs.setTemperature(T);
        fs.setSaturation(wPhaseIdx, 1 - Sn);
        fs.setSaturation(nPhaseIdx, Sn);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(nPhaseIdx, pn);
        fs.setPressure(nPhaseIdx, pn + (pC[wPhaseIdx] - pC[nPhaseIdx]));
        
        assignNaive(fs);
    };

    /*!
     * \brief Set the primary variables from an arbitrary fluid state
     *        in a mass conservative way.
     *
     * If an energy equation is included, the fluid temperatures are
     * the same as the one given in the fluid state, *not* the
     * enthalpy.
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     * \param matParams The capillary pressure law parameters
     * \param isInEquilibrium If true, the fluid state expresses
     *                        thermodynamic equilibrium assuming the
     *                        relations expressed by the fluid
     *                        system. This implies that in addition to
     *                        the quantities mentioned above, the
     *                        fugacities are also defined.
     */
    template <class MaterialLaw, class FluidState>
    void assignMassConservative(const FluidState &fluidState,
                                const typename MaterialLaw::Params &matParams,
                                bool isInEquilibrium = false)
    {
        ComponentVector globalMolarities(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                globalMolarities[compIdx] += 
                    fluidState.molarity(phaseIdx, compIdx)
                    * fluidState.saturation(phaseIdx);
            }
        }
        
        ImmiscibleFluidState<Scalar, FluidSystem> fsFlash;
        fsFlash.assign(fluidState);
        typename FluidSystem::ParameterCache paramCache;
        ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams, globalMolarities);
        
        assignNaive(fsFlash);
    };

    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);
        
        if (formulation == pwSn) {
            (*this)[pressureIdx] = fluidState.pressure(wPhaseIdx);
            (*this)[saturationIdx] = fluidState.saturation(nPhaseIdx);
        }
        else if (formulation == pnSw) {
            (*this)[pressureIdx] = fluidState.pressure(nPhaseIdx);
            (*this)[saturationIdx] = fluidState.saturation(wPhaseIdx);
        }
        else {
            // invalid formulation
            assert(false);
        }
    }
};

} // end namepace

#endif
