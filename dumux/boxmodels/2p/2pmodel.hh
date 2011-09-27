// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2007 by Bernd Flemisch                               *
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
* \brief Adaption of the box scheme to the two-phase flow model.
*/

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include "2pproperties.hh"
#include "2plocalresidual.hh"
#include "2pproblem.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPBoxModel
 * \brief A two-phase, isothermal flow model using the box scheme.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * These equations are discretized by a fully-coupled vertex centered finite volume
 * (box) scheme as spatial and the implicit Euler method as time
 * discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPCommonIndices::pWsN</tt> or <tt>TwoPCommonIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class TwoPModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::ctype CoordScalar;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global index of the vertex
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (pressureIdx == pvIdx) {
            Scalar absPv = 
                std::abs(this->solution(/*historyIdx=*/1)[globalVertexIdx][pvIdx]);
            return std::min(10.0/absPv, 1.0);
        }
        return 1;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The global solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        bool velocityOutput = GET_PROP_VALUE(TypeTag, EnableVelocityOutput);
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *pW = writer.allocateManagedBuffer(numVertices);
        ScalarField *pN = writer.allocateManagedBuffer(numVertices);
        ScalarField *pC = writer.allocateManagedBuffer(numVertices);
        ScalarField *Sw = writer.allocateManagedBuffer(numVertices);
        ScalarField *Sn = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoW = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoN = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobW = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobN = writer.allocateManagedBuffer(numVertices);
        ScalarField *poro = writer.allocateManagedBuffer(numVertices);
        ScalarField *Te = writer.allocateManagedBuffer(numVertices);
        ScalarField *cellNum =writer.allocateManagedBuffer (numVertices);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numVertices);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numVertices);

        if(velocityOutput) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (int i = 0; i < numVertices; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
                (*cellNum)[i] = Scalar(0.0);
            }
        }
        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();

            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*historyIdx=*/0);

            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, scvIdx, dim);

                const VolumeVariables &volVars = elemCtx.volVars(scvIdx, /*historyIdx=*/0);
                (*pW)[globalIdx] = volVars.pressure(wPhaseIdx);
                (*pN)[globalIdx] = volVars.pressure(nPhaseIdx);
                (*pC)[globalIdx] = volVars.capillaryPressure();
                (*Sw)[globalIdx] = volVars.saturation(wPhaseIdx);
                (*Sn)[globalIdx] = volVars.saturation(nPhaseIdx);
                (*rhoW)[globalIdx] = volVars.density(wPhaseIdx);
                (*rhoN)[globalIdx] = volVars.density(nPhaseIdx);
                (*mobW)[globalIdx] = volVars.mobility(wPhaseIdx);
                (*mobN)[globalIdx] = volVars.mobility(nPhaseIdx);
                (*poro)[globalIdx] = volVars.porosity();
                (*Te)[globalIdx] = volVars.temperature();
                if(velocityOutput)
                {
                    (*cellNum)[globalIdx] += 1;
                }
            };

        }

        writer.attachVertexData(*Sn, "Sn");
        writer.attachVertexData(*Sw, "Sw");
        writer.attachVertexData(*pN, "pn");
        writer.attachVertexData(*pW, "pw");
        writer.attachVertexData(*pC, "pc");
        writer.attachVertexData(*rhoW, "rhoW");
        writer.attachVertexData(*rhoN, "rhoN");
        writer.attachVertexData(*mobW, "mobW");
        writer.attachVertexData(*mobN, "mobN");
        writer.attachVertexData(*poro, "porosity");
        writer.attachVertexData(*Te, "temperature");

        writer.attachCellData(*rank, "process rank");
    }
};
}

#include "2ppropertydefaults.hh"

#endif
