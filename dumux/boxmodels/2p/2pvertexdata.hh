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
        
        // coordinates of the vertex
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                dim);

        if (formulation == I::pWsN) {
            satN = sol[I::saturationIdx];
            satW = 1.0 - satN;
            pC = problem.materialLaw().pC(satW,
                                          global,
                                          element,
                                          local);
            pressure[I::wPhase] = sol[I::pressureIdx];
            pressure[I::nPhase] = pressure[I::wPhase] + pC;
        }
        else if (formulation == I::pNsW) {
            satW = sol[I::saturationIdx];
            satN = 1.0 - satW;
            pC =problem.materialLaw().pC(satW,
                                               global,
                                               element,
                                               local);
            pressure[I::nPhase] = sol[I::pressureIdx];
            pressure[I::wPhase] = pressure[I::nPhase] - pC;
        }

        density[I::wPhase] = problem.wettingPhase().density(temperature,
                                                           pressure[I::wPhase]);
        density[I::nPhase] = problem.nonwettingPhase().density(temperature,
                                                              pressure[I::nPhase]);

        mobility[I::wPhase] =problem.materialLaw().mobW(satW,
                                                        global,
                                                        element,
                                                        local,
                                                        temperature,
                                                        pressure[I::wPhase]);
        mobility[I::nPhase] =problem.materialLaw().mobN(satN,
                                                        global,
                                                        element,
                                                        local,
                                                        temperature,
                                                        pressure[I::nPhase]);

        // porosity
        porosity = problem.soil().porosity(global,
                                           element,
                                           local);
    }

    void updateTemperature_(const PrimaryVarVector  &sol,
                            const Element           &element,
                            const FVElementGeometry &elemGeom,
                            int                      vertIdx,
                            const Problem           &problem) 
    {
        temperature = problem.temperature(element, elemGeom, vertIdx);
    }
        
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
