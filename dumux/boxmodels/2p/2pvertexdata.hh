// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Quantities required by the twophase box model defined on a vertex.
 */
#ifndef DUMUX_2P_VERTEX_DATA_HH
#define DUMUX_2P_VERTEX_DATA_HH

#include "2pproperties.hh"

namespace Dune
{

/*!
 * \ingroup TwoPBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPVertexData
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

        formulation   = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))     SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem))        FluidSystem;
    typedef TwoPPhaseState<TypeTag>                           PhaseState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw))        MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams))  MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename SolutionTypes::PrimaryVarVector  PrimaryVarVector;
    typedef Dune::FieldVector<Scalar, numPhases>      PhasesVector;

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
        
        // material law parameters
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, vertIdx);

        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);

        if (formulation == I::pWsN) {
            satN = sol[I::saturationIdx];
            satW = 1.0 - satN;
            pC = MaterialLaw::pC(materialParams, satW);
            pressure[I::wPhase] = sol[I::pressureIdx];
            pressure[I::nPhase] = pressure[I::wPhase] + pC;
        }
        else if (formulation == I::pNsW) {
            satW = sol[I::saturationIdx];
            satN = 1.0 - satW;
            pC = MaterialLaw::pC(materialParams, satW);
            pressure[I::nPhase] = sol[I::pressureIdx];
            pressure[I::wPhase] = pressure[I::nPhase] - pC;
        }

        phaseState.update(pressure[I::wPhase], pressure[I::nPhase], temperature);

        density[I::wPhase] = FluidSystem::phaseDensity(I::wPhase, phaseState);
        density[I::nPhase] = FluidSystem::phaseDensity(I::nPhase, phaseState);

        mobility[I::wPhase] =MaterialLaw::krw(materialParams, satW)/FluidSystem::phaseViscosity(I::wPhase, phaseState);
        mobility[I::nPhase] = MaterialLaw::krn(materialParams, satW)/FluidSystem::phaseViscosity(I::nPhase, phaseState);

        // porosity
        porosity = problem.spatialParameters().porosity(global,
                                           element);
    }

    void updateTemperature_(const PrimaryVarVector  &sol,
                            const Element           &element,
                            const FVElementGeometry &elemGeom,
                            int                      vertIdx,
                            const Problem           &problem) 
    {
        temperature = problem.temperature(element, elemGeom, vertIdx);
    }
        

    PhaseState phaseState;
    Scalar satW;
    Scalar satN;
    Scalar pC;
    Scalar porosity;
    Scalar temperature;

    PhasesVector density;
    PhasesVector pressure;
    PhasesVector mobility;

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }
    
    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
