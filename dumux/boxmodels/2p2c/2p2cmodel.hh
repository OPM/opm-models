// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 */
#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include "2p2cproperties.hh"
#include "2p2clocalresidual.hh"
#include "2p2cproblem.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{{\bf D}_{\alpha, pm}^\kappa \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^a_w<X^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^w_n<X^w_{n,max})\f$</li>
 * </ul>
 */
template<class TypeTag>
class TwoPTwoCModel: public BoxModel<TypeTag>
{
    typedef TwoPTwoCModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        staticVertexDat_.resize(this->gridView_().size(dim));

        setSwitched_(false);

        ElementContext elemCtx(this->problem_());

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator it = this->gridView_().template begin<0>();
        const ElementIterator &eendit = this->gridView_().template end<0>();
        for (; it != eendit; ++it) {
            // deal with the current element
            elemCtx.updateFVElemGeom(*it);

            // loop over all element vertices, i.e. sub control volumes
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
            {
                // map the local vertex index to the global one
                int globalIdx = this->vertexMapper().map(*it,
                                                         scvIdx,
                                                         dim);


                // initialize phase presence
                staticVertexDat_[globalIdx].phasePresence[/*historyIdx=*/0]
                    = this->problem_().initialPhasePresence(elemCtx, scvIdx);
                staticVertexDat_[globalIdx].wasSwitched = false;

                staticVertexDat_[globalIdx].phasePresence[/*historyIdx=*/1]
                    = staticVertexDat_[globalIdx].phasePresence[/*historyIdx=*/0];
            }
        }
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \param dest Contains the storage of each component for one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &dest, int phaseIdx)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*historyIdx=*/0);

            this->localResidual().addPhaseStorage(dest, elemCtx, /*historyIdx=*/0, phaseIdx);
        };

        this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        setSwitched_(false);
        resetPhasePresence_();
    };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return std::min(1.0/this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx], 1.0);
        return 1;
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();

        // update the phase state
        updateOldPhasePresence_();
        setSwitched_(false);
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Returns the phase presence of the current or the old solution of a vertex.
     *
     * \param globalVertexIdx The global vertex index
     * \param oldSol Evaluate function with solution of current or previous time step
     */
    int phasePresence(int globalVertexIdx, int historyIdx) const
    { return staticVertexDat_[globalVertexIdx].phasePresence[historyIdx]; }


    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param vert The vertex
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outStream << staticVertexDat_[vertIdx].phasePresence[/*historyIdx=*/0] << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param vert The vertex
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        // read phase presence
        int vertIdx = this->dofMapper().map(vert);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        inStream >> staticVertexDat_[vertIdx].phasePresence[/*historyIdx=*/0];
        staticVertexDat_[vertIdx].phasePresence[/*historyIdx=*/1]
            = staticVertexDat_[vertIdx].phasePresence[/*historyIdx=*/0];
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     *
     * \param curGlobalSol The current global solution
     * \param oldGlobalSol The previous global solution
     */
    void updateStaticData_(SolutionVector &curGlobalSol,
                           const SolutionVector &oldGlobalSol)
    {
        bool wasSwitched = false;

        for (unsigned i = 0; i < staticVertexDat_.size(); ++i)
            staticVertexDat_[i].visited = false;

        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            bool fvElemGeomUpdated = false;
            int numScv = elemIt->template count<dim>();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, scvIdx, dim);

                if (staticVertexDat_[globalIdx].visited)
                    continue;
                staticVertexDat_[globalIdx].visited = true;

                if (!fvElemGeomUpdated) {
                    fvElemGeomUpdated = true;
                    elemCtx.updateFVElemGeom(*elemIt);
                }

                // compute the volume variables of the current
                // sub-control volume
                elemCtx.updateScvVars(curGlobalSol[globalIdx],
                                       scvIdx,
                                       /*historyIdx=*/0);
                const VolumeVariables &volVars = elemCtx.volVars(scvIdx, /*historyIdx=*/0);

                const GlobalPosition &globalPos = elemCtx.pos(scvIdx);
                if (primaryVarSwitch_(curGlobalSol,
                                      volVars,
                                      globalIdx,
                                      globalPos))
                {
                    this->jacobianAssembler().markVertexRed(globalIdx);
                    wasSwitched = true;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        if (this->gridView_().comm().size() > 1)
            wasSwitched = this->gridView_().comm().max(wasSwitched);

        setSwitched_(wasSwitched);
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    };

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVars
    {
        int phasePresence[historySize];
        bool wasSwitched;
        bool visited;
    };

    /*!
     * \brief Reset the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence_()
    {
        int numVertices = this->gridView_().size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].phasePresence[/*historyIdx=*/0]
                = staticVertexDat_[i].phasePresence[/*historyIdx=*/1];
            staticVertexDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhasePresence_()
    {
        int numVertices = this->gridView_().size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].phasePresence[/*historyIdx=*/1]
                = staticVertexDat_[i].phasePresence[/*historyIdx=*/0];
            staticVertexDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    void setSwitched_(bool yesno)
    { switchFlag_ = yesno; }

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SolutionVector &globalSol,
                           const VolumeVariables &volVars,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = staticVertexDat_[globalIdx].phasePresence[/*historyIdx=*/0];
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == gPhaseOnly)
        {
            // calculate mole fraction in the hypothetic liquid phase
            Scalar xll = volVars.fluidState().moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
            Scalar xlg = volVars.fluidState().moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);

            Scalar xlMax = 1.0;
            if (xll + xlg > xlMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xlMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, liquid phase appears
            if (xll + xlg > xlMax)
            {
                // liquid phase appears
                std::cout << "liquid phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xll + xlg: "
                          << xll + xlg << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    globalSol[globalIdx][switchIdx] = 0.0;
                else if (formulation == plSg)
                    globalSol[globalIdx][switchIdx] = 1.0;
            };
        }
        else if (phasePresence == lPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic gas phase
            Scalar xgl = volVars.fluidState().moleFraction(/*phaseIdx=*/1, /*compIdx=*/0);
            Scalar xgg = volVars.fluidState().moleFraction(/*phaseIdx=*/1, /*compIdx=*/1);

            Scalar xgMax = 1.0;
            if (xgl + xgg > xgMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xgl + xgg > xgMax)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", x_gl + x_gg: "
                          << xgl + xgg << std::endl;
                std::cout.flush();
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    globalSol[globalIdx][switchIdx] = 0.999;
                else if (formulation == plSg)
                    globalSol[globalIdx][switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0;
            if (staticVertexDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.fluidState().saturation(/*phaseIdx=*/1) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sg: "
                          << volVars.fluidState().saturation(/*phaseIdx=*/1) << std::endl;
                newPhasePresence = lPhaseOnly;

                globalSol[globalIdx][switchIdx] =
                    volVars.fluidState().massFraction(/*phaseIdx=*/0, /*compIdx=*/1);
            }
            else if (volVars.fluidState().saturation(/*phaseIdx=*/0) <= Smin)
            {
                wouldSwitch = true;
                // liquid phase disappears
                std::cout << "Liquid phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sl: "
                          << volVars.fluidState().saturation(/*phaseIdx=*/0) << std::endl;
                newPhasePresence = gPhaseOnly;

                globalSol[globalIdx][switchIdx] =
                    volVars.fluidState().massFraction(/*phaseIdx=*/1, /*compIdx=*/0);
            }
        }

        staticVertexDat_[globalIdx].phasePresence[/*historyIdx=*/0] = newPhasePresence;
        staticVertexDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

protected:
    // parameters given in constructor
    std::vector<StaticVars> staticVertexDat_;
    bool switchFlag_;
};

}

#include "2p2cpropertydefaults.hh"

#endif
