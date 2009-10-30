// $Id$
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
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
 * \ingroup TwoPTwoCBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_VERTEX_DATA_HH
#define DUMUX_2P2C_VERTEX_DATA_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/material/multicomponentrelations.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class TwoPTwoCVertexData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;

    // this is a bit hacky: the Vertex data might not be identical to
    // the implementation.
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   Implementation;

    enum {
        numEq         = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases     = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        formulation   = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    typedef typename SolutionTypes::PrimaryVarVector               PrimaryVarVector;
    typedef Dune::FieldVector<Scalar, numPhases>                   PhasesVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MultiComp))       MultiComp;

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const PrimaryVarVector  &sol,
                const Element           &element,
                const FVElementGeometry &elemGeom,
                int                      vertIdx,
                const Problem           &problem,
                bool                     isOldSol) 
    {
        typedef Indices I;

        asImp().updateTemperature_(sol,
                                   element, 
                                   elemGeom,
                                   vertIdx,
                                   problem);
        
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                dim);
        int globalVertIdx = problem.model().dofEntityMapper().map(element, vertIdx, dim);
        int phaseState = problem.model().localJacobian().phaseState(globalVertIdx, isOldSol);

        // update saturations, pressures and mass fractions
        updateSaturations(sol, phaseState);
        Scalar pC = problem.materialLaw().pC(saturation[I::wPhase],
                                             global,
                                             element,
                                             local,
                                             temperature);
        updatePressures(sol, pC);
        updateMassFracs(sol, problem.multicomp(), phaseState, temperature);

        // Densities
        density[I::wPhase] = problem.wettingPhase().density(temperature,
                                                            pressure[I::wPhase],
                                                            massfrac[I::nComp][I::wPhase]);
        density[I::nPhase] = problem.nonwettingPhase().density(temperature,
                                                               pressure[I::nPhase],
                                                               massfrac[I::wComp][I::nPhase]);

        // Mobilities
        mobility[I::wPhase] = problem.materialLaw().mobW(saturation[I::wPhase],
                                                         global,
                                                         element,
                                                         local,
                                                         temperature,
                                                         pressure[I::wPhase]);
        mobility[I::nPhase] = problem.materialLaw().mobN(saturation[I::nPhase],
                                                         global,
                                                         element,
                                                         local,
                                                         temperature,
                                                         pressure[I::nPhase]);

        // binary diffusion coefficents
        diffCoeff[I::wPhase] =
            problem.wettingPhase().diffCoeff(temperature, pressure[I::wPhase]);
        diffCoeff[I::nPhase] =
            problem.nonwettingPhase().diffCoeff(temperature, pressure[I::nPhase]);
        
        // porosity
        porosity = problem.soil().porosity(global,
                                           element,
                                                 local);
    }

    void updateTemperature_(const PrimaryVarVector  &sol,
                            const Element           &element,
                            const FVElementGeometry &elemGeom,
                            int                      vertIdx,
                            Problem                 &problem) 
    {
        temperature = problem.temperature(element, elemGeom, vertIdx);
    }

    /*!
     * \brief Update saturations.
     */
    void updateSaturations(const PrimaryVarVector &sol,
                           int phaseState)
    {
        typedef Indices I;

        // update saturations
        if (formulation == I::pWsN)
        {
            if (phaseState == I::bothPhases)
                saturation[I::nPhase] = sol[I::switchIdx];
            else if (phaseState == I::wPhaseOnly)
                saturation[I::nPhase] = 0.0;
            else if (phaseState == I::nPhaseOnly)
                saturation[I::nPhase] = 1.0;
            else
                DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            saturation[I::wPhase] = 1.0 - saturation[I::nPhase];
        }
        else if (formulation == I::pNsW)
        {
            if (phaseState == I::bothPhases)
                saturation[I::wPhase] = sol[I::switchIdx];
            else if (phaseState == I::wPhaseOnly)
                saturation[I::wPhase] = 1.0;
            else if (phaseState == I::nPhaseOnly)
                saturation[I::wPhase] = 0.0;
            else
                DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            saturation[I::nPhase] = 1.0 - saturation[I::wPhase];
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
    }

    /*!
     * \brief Update phase pressures.
     *
     * Requires up to date saturations.
     */
    void updatePressures(const PrimaryVarVector &sol,
                         Scalar capillaryPressure)
    {
        typedef Indices I;

        // update pressures
        pC = capillaryPressure;

        if (formulation == I::pWsN) {
            pressure[I::wPhase] = sol[I::pressureIdx];
            pressure[I::nPhase] = pressure[I::wPhase] + pC;
        }
        else if (formulation == I::pNsW) {
            pressure[I::nPhase] = sol[I::pressureIdx];
            pressure[I::wPhase] = pressure[I::nPhase] - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
    }

    /*!
     * \brief Update mass fraction matrix.
     *
     * Requires up to date pressures and saturations.
     */
    void updateMassFracs(const PrimaryVarVector &sol,
                         const MultiComp &multicomp,
                         int phaseState,
                         Scalar temperature)
    {
        typedef Indices I;

        // Solubilities of components in phases
        if (phaseState == I::bothPhases) {
            massfrac[I::wComp][I::nPhase] = multicomp.xWN(pressure[I::nPhase], temperature);
            massfrac[I::nComp][I::wPhase] = multicomp.xAW(pressure[I::nPhase], temperature);
        }
        else if (phaseState == I::wPhaseOnly) {
            massfrac[I::wComp][I::nPhase] = 0.0;
            massfrac[I::nComp][I::wPhase] = sol[I::switchIdx];
        }
        else if (phaseState == I::nPhaseOnly) {
            massfrac[I::wComp][I::nPhase] = sol[I::switchIdx];
            massfrac[I::nComp][I::wPhase] = 0.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

        massfrac[I::wComp][I::wPhase] = 1.0 - massfrac[I::nComp][I::wPhase];
        massfrac[I::nComp][I::nPhase] = 1.0 - massfrac[I::wComp][I::nPhase];
    };

public:
    PhasesVector saturation;//!< Effective phase saturations within the control volume
    PhasesVector pressure;  //!< Effective phase pressures within the control volume
    Scalar temperature;     //!< Temperature within the control volume
    Scalar pC;              //!< Effective capillary pressure within the control volume
    Scalar porosity;        //!< Effective porosity within the control volume
    PhasesVector mobility;  //!< Effective mobility within the control volume
    PhasesVector density;   //!< Effective density within the control volume
    PhasesVector diffCoeff; //!< Binary diffusion coefficients for the phases

    //! Mass fractions of each component within each phase
    Dune::FieldMatrix<Scalar, numComponents, numPhases> massfrac;

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }
    
    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
