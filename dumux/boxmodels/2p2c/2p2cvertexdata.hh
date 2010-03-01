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

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include "2p2cphasestate.hh"

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

    // this is a bit hacky: the vertex data might not be identical to
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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))           Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))   Indices;

    enum {
        lCompIdx  = Indices::lCompIdx,
        gCompIdx  = Indices::gCompIdx,

        lPhaseIdx  = Indices::lPhaseIdx,
        gPhaseIdx  = Indices::gPhaseIdx
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))        SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector               PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem))     FluidSystem;
    typedef TwoPTwoCPhaseState<TypeTag>                            PhaseState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw))       MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    void update(const PrimaryVarVector  &sol,
                const Element           &element,
                const FVElementGeometry &elemGeom,
                int                      scvIdx,
                const Problem           &problem,
                bool                     isOldSol) 
    {
        asImp().updateTemperature_(sol,
                                   element, 
                                   elemGeom,
                                   scvIdx,
                                   problem);

        // capillary pressure parameters
        const MaterialLawParams &materialParams = 
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);
        
        int globalVertIdx = problem.model().dofEntityMapper().map(element, scvIdx, dim);
        int phasePresence = problem.model().localJacobian().phasePresence(globalVertIdx, isOldSol);

        // calculate phase state
        phaseState_.update(sol, materialParams, temperature(), phasePresence);
        Valgrind::CheckDefined(phaseState_);

        // Mobilities
        Scalar muL = FluidSystem::phaseViscosity(lPhaseIdx, 
                                                 phaseState().temperature(),
                                                 phaseState().phasePressure(lPhaseIdx),
                                                 phaseState());
        Scalar muG = FluidSystem::phaseViscosity(gPhaseIdx,
                                                 phaseState().temperature(),
                                                 phaseState().phasePressure(gPhaseIdx),
                                                 phaseState());
        mobility_[lPhaseIdx] = Scalar(1.0) / muL * MaterialLaw::krw(materialParams, saturation(lPhaseIdx));
        mobility_[gPhaseIdx] = Scalar(1.0) / muG * MaterialLaw::krn(materialParams, saturation(lPhaseIdx));
        Valgrind::CheckDefined(mobility_[lPhaseIdx]);
        Valgrind::CheckDefined(mobility_[gPhaseIdx]);

        // binary diffusion coefficents
        diffCoeff_[lPhaseIdx] = 
            FluidSystem::diffCoeff(lPhaseIdx, 
                                   lCompIdx,
                                   gCompIdx,
                                   phaseState_.temperature(),
                                   phaseState_.phasePressure(lPhaseIdx),
                                   phaseState_);
        diffCoeff_[gPhaseIdx] = 
            FluidSystem::diffCoeff(gPhaseIdx, 
                                   lCompIdx, 
                                   gCompIdx,
                                   phaseState_.temperature(),
                                   phaseState_.phasePressure(gPhaseIdx),
                                   phaseState_);
        Valgrind::CheckDefined(diffCoeff_);
        
        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);
   }

    void updateTemperature_(const PrimaryVarVector  &sol,
                            const Element           &element,
                            const FVElementGeometry &elemGeom,
                            int                      scvIdx,
                            const Problem           &problem) 
    {
        temperature_ = problem.temperature(element, elemGeom, scvIdx);
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const PhaseState &phaseState() const
    { return phaseState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     */
    Scalar saturation(int phaseIdx) const
    { return phaseState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     */
    Scalar density(int phaseIdx) const
    { return phaseState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     */
    Scalar pressure(int phaseIdx) const
    { return phaseState_.phasePressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     */
    Scalar mobility(int phaseIdx) const
    { 
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return phaseState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase
     */
    Scalar diffCoeff(int phaseIdx) const
    { return diffCoeff_[phaseIdx]; }


protected:
    Scalar temperature_;     //!< Temperature within the control volume
    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    PhaseState phaseState_;

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }
    
    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
