// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Melanie Darcis 								 *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
#ifndef DUMUX_2P2CNI_VERTEX_DATA_HH
#define DUMUX_2P2CNI_VERTEX_DATA_HH

#include <dumux/boxmodels/2p2c/2p2cvertexdata.hh>

namespace Dune
{

/*!
 * \ingroup TwoPTwoCNIBoxModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
template <class TypeTag>
class TwoPTwoCNIVertexData : public TwoPTwoCVertexData<TypeTag>
{
    typedef TwoPTwoCVertexData<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum {
        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,

        numPhases     = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))                SolutionTypes;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements))::Container ReferenceElements;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
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
                Problem                 &problem,
                bool                     isOldSol) 
    {
        typedef Indices I;

        // vertex update data for the mass balance
        ParentType::update(sol,
                           element,
                           elemGeom,
                           vertIdx,
                           problem,
                           isOldSol);

        // data for the energy equation
        const LocalPosition &local =
            ReferenceElements::general(element.type()).position(vertIdx,
                                                                dim);
        const GlobalPosition &global =
            element.geometry().corner(vertIdx);

        heatCond = problem.soil().heatCond(global,
                                           element,
                                           local,
                                           this->saturation[I::wPhase]);
        
        enthalpy[I::wPhase] = problem.wettingPhase().enthalpy(this->temperature,
                                                              this->pressure[I::wPhase],
                                                              this->massfrac[I::nComp][I::wPhase]);
        enthalpy[I::nPhase] = problem.nonwettingPhase().enthalpy(this->temperature,
                                                                 this->pressure[I::nPhase],
                                                                 this->massfrac[I::wComp][I::nPhase]);
        intEnergy[I::wPhase] = problem.wettingPhase().intEnergy(this->temperature,
                                                                this->pressure[I::wPhase],
                                                                this->massfrac[I::nComp][I::wPhase]);
        intEnergy[I::nPhase] = problem.nonwettingPhase().intEnergy(this->temperature,
                                                                   this->pressure[I::nPhase],
                                                                   this->massfrac[I::wComp][I::nPhase]);
    }

    // this method gets called by the parent class
    void updateTemperature_(const PrimaryVarVector  &sol,
                            const Element           &element,
                            const FVElementGeometry &elemGeom,
                            int                      vertIdx,
                            const Problem           &problem)
    {
        typedef Indices I;
        this->temperature = sol[I::temperatureIdx];
    }


    PhasesVector intEnergy; //!< Internal energy.
    PhasesVector enthalpy;  //!< Enthalpy.
    Scalar       heatCond; //!< Total heat conductivity.
};

} // end namepace

#endif
