// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief Some simple test problem for the black-oil box model
 *        inspired by a simple reservoir.
 */
#ifndef DUMUX_RESERVOIR_PROBLEM_HH
#define DUMUX_RESERVOIR_PROBLEM_HH

#include <dumux/boxmodels/blackoil/blackoilmodel.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

#include <iostream>
#include <string>

namespace Dumux
{
template <class TypeTag>
class ReservoirProblem;

namespace Properties
{
NEW_TYPE_TAG(ReservoirProblem, INHERITS_FROM(BoxBlackOil));

NEW_PROP_TAG(MaxDepth);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);

// Set the grid type
SET_TYPE_PROP(ReservoirProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ReservoirProblem, Problem, Dumux::ReservoirProblem<TypeTag>);

// Set the material Law
SET_PROP(ReservoirProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef MpLinearMaterial<FluidSystem::numPhases, Scalar> type;
};

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(ReservoirProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(ReservoirProblem, EnableGravity, true);

// Reuse Jacobian matrices if possible?
SET_BOOL_PROP(ReservoirProblem, EnableJacobianRecycling, true);

// Smoothen the upwinding method?
SET_BOOL_PROP(ReservoirProblem, EnableSmoothUpwinding, false);

// Enable constraint DOFs?
SET_BOOL_PROP(ReservoirProblem, EnableConstraints, true);

// set the defaults for some problem specific properties
SET_SCALAR_PROP(ReservoirProblem, MaxDepth, 2500);
SET_SCALAR_PROP(ReservoirProblem, Temperature, 293.15);
SET_STRING_PROP(ReservoirProblem, SimulationName, "reservoir");
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxTestProblems
 *
 * \brief Some simple test problem for the black-oil box model
 *        inspired by a simple reservoir.
 */
template <class TypeTag>
class ReservoirProblem 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,
        numComponents = FluidSystem::numComponents,

        gPhaseIdx = FluidSystem::gPhaseIdx,
        oPhaseIdx = FluidSystem::oPhaseIdx,
        wPhaseIdx = FluidSystem::wPhaseIdx,

        gCompIdx = FluidSystem::gCompIdx,
        oCompIdx = FluidSystem::oCompIdx,
        wCompIdx = FluidSystem::wCompIdx,

        conti0EqIdx = Indices::conti0EqIdx,
        pressure0Idx = Indices::pressure0Idx,
        saturation0Idx = Indices::saturation0Idx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) BlackOilFluidState;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ReservoirProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {      
        eps_ = 1e-6;

        temperature_ = GET_PARAM(TypeTag, Scalar, Temperature);
        maxDepth_ = GET_PARAM(TypeTag, Scalar, MaxDepth);
        name_ = GET_PARAM(TypeTag, std::string, SimulationName);

        FluidSystem::initBegin();
        std::vector<std::pair<Scalar, Scalar> > Bg = {
            { 1.013529e+05,  9.998450e-01 },
            { 2.757903e+06,  3.075500e-02 },
            { 5.515806e+06,  1.537947e-02 },
            { 8.273709e+06,  1.021742e-02 },
            { 1.103161e+07,  7.662783e-03 },
            { 1.378951e+07,  6.151899e-03 },
            { 1.654742e+07,  5.108709e-03 },
            { 1.930532e+07,  4.378814e-03 },
            { 2.206322e+07,  3.857780e-03 },
            { 2.482113e+07,  3.388401e-03 },
            { 2.757903e+07,  3.049842e-03 }
        };
        std::vector<std::pair<Scalar, Scalar> > Bo = {
            { 1.013529e+05, 1.000000e+00 },
            { 2.757903e+06, 1.012000e+00 },
            { 5.515806e+06, 1.025500e+00 },
            { 8.273709e+06, 1.038000e+00 },
            { 1.103161e+07, 1.051000e+00 },
            { 1.378951e+07, 1.063000e+00 },
            { 1.654742e+07, 1.075000e+00 },
            { 1.930532e+07, 1.087000e+00 },
            { 2.206322e+07, 1.098500e+00 },
            { 2.482113e+07, 1.110000e+00 },
            { 2.757903e+07, 1.120000e+00 }
        };
        std::vector<std::pair<Scalar, Scalar> > Rs = {
            { 1.013529e+05, 0.000000e+00 },
            { 2.757903e+06, 2.938776e+01 },
            { 5.515806e+06, 5.966605e+01 },
            { 8.273709e+06, 8.905380e+01 },
            { 1.103161e+07, 1.184416e+02 },
            { 1.378951e+07, 1.474731e+02 },
            { 1.654742e+07, 1.754360e+02 },
            { 1.930532e+07, 2.012616e+02 },
            { 2.206322e+07, 2.261967e+02 },
            { 2.482113e+07, 2.475696e+02 },
            { 2.757903e+07, 2.671614e+02 }
        };
        std::vector<std::pair<Scalar, Scalar> > muo = {
            { 1.013529e+05, 1.200000e-03 },
            { 2.757903e+06, 1.170000e-03 },
            { 5.515806e+06, 1.140000e-03 },
            { 8.273709e+06, 1.110000e-03 },
            { 1.103161e+07, 1.080000e-03 },
            { 1.378951e+07, 1.060000e-03 },
            { 1.654742e+07, 1.030000e-03 },
            { 1.930532e+07, 1.000000e-03 },
            { 2.206322e+07, 9.800000e-04 },
            { 2.482113e+07, 9.500000e-04 },
            { 2.757903e+07, 9.400000e-04 }
        };
        std::vector<std::pair<Scalar, Scalar> > mug = {
            { 1.013529e+05, 1.250000e-05 },
            { 2.757903e+06, 1.300000e-05 },
            { 5.515806e+06, 1.350000e-05 },
            { 8.273709e+06, 1.400000e-05 },
            { 1.103161e+07, 1.450000e-05 },
            { 1.378951e+07, 1.500000e-05 },
            { 1.654742e+07, 1.550000e-05 },
            { 1.930532e+07, 1.600000e-05 },
            { 2.206322e+07, 1.650000e-05 },
            { 2.482113e+07, 1.700000e-05 },
            { 2.757903e+07, 1.750000e-05 },
        };
        FluidSystem::setGasFormationVolumeFactor(Bg);
        FluidSystem::setOilFormationVolumeFactor(Bo);
        FluidSystem::setGasFormationFactor(Rs);
        FluidSystem::setOilViscosity(muo);
        FluidSystem::setGasViscosity(mug);
        FluidSystem::setWaterViscosity(9.6e-4);
        FluidSystem::setWaterCompressibility(1.450377e-10);
        FluidSystem::setSurfaceDensities(/*oil=*/720.51,
                                         /*water=*/1009.32,
                                         /*gas=*/1.1245);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            FluidSystem::setReferenceVolumeFactor(phaseIdx, 1.0);
        FluidSystem::initEnd();

        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-12);
        coarseK_ = this->toDimMatrix_(1e-11);

        // porosities
        finePorosity_ = 0.2;
        coarsePorosity_ = 0.3;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            fineMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            fineMaterialParams_.setPcMaxSat(phaseIdx, 0.0);
            fineMaterialParams_.setResidSat(phaseIdx, 0.0);

            coarseMaterialParams_.setPcMinSat(phaseIdx, 0.0);
            coarseMaterialParams_.setPcMaxSat(phaseIdx, 0.0);
            coarseMaterialParams_.setResidSat(phaseIdx, 0.0);
        }

        initFluidState_();
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \brief Return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Evaluate the boundary conditions for a boundary segment.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        // no flow on top and bottom
        values.setNoFlow();
    }
    
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Evaluate the boundary conditions for a boundary segment.
     */
    template <class Context>
    void constraints(Constraints &values,
                     const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        Scalar x = pos[0] - this->bboxMin()[0];
        Scalar y = pos[dim-1] - this->bboxMin()[dim-1];
        Scalar height = this->bboxMax()[dim-1] - this->bboxMin()[dim-1];
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        if ((onLeftBoundary_(pos)
             || onRightBoundary_(pos))
            && y < height/2)
        {
            auto fs = initialFluidState_;

            Scalar pwInj = fs.pressure(wPhaseIdx)*1.5;
            fs.setPressure(wPhaseIdx, pwInj);
            fs.setPressure(oPhaseIdx, pwInj);
            fs.setPressure(gPhaseIdx, pwInj);
            fs.setSaturation(wPhaseIdx, 1.0);
            fs.setSaturation(oPhaseIdx, 0.0);
            fs.setSaturation(gPhaseIdx, 0.0);

            values.setConstraint(pressure0Idx, conti0EqIdx, fs.pressure(/*phaseIdx=*/0));
            values.setConstraint(saturation0Idx, conti0EqIdx + 1, fs.saturation(/*phaseIdx=*/0));
            values.setConstraint(saturation0Idx + 1, conti0EqIdx + 2, fs.saturation(/*phaseIdx=*/1));
        }
        else if (width/2 - 1 < x && x < width/2 + 1 && y > height/2)
        {
            auto fs = initialFluidState_;

            Scalar pwInj = fs.pressure(wPhaseIdx)/1.5;
            fs.setPressure(wPhaseIdx, pwInj);
            fs.setPressure(oPhaseIdx, pwInj);
            fs.setPressure(gPhaseIdx, pwInj);
            fs.setSaturation(wPhaseIdx, 0.0);
            fs.setSaturation(oPhaseIdx, 1.0);
            fs.setSaturation(gPhaseIdx, 0.0);

            values.setConstraint(pressure0Idx, conti0EqIdx, fs.pressure(/*phaseIdx=*/0));
            values.setConstraint(saturation0Idx, conti0EqIdx + 1, fs.saturation(/*phaseIdx=*/0));
            values.setConstraint(saturation0Idx + 1, conti0EqIdx + 2, fs.saturation(/*phaseIdx=*/1));
        }
    }
    
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        //////
        // set the primary variables
        //////
        values.assignNaive(initialFluidState_);
    }

    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { values = 0; }

    // \}
private:
    void initFluidState_()
    {
        auto &fs = initialFluidState_;
        const auto &pos = this->bboxMin();

        //////
        // set temperatures
        //////
        fs.setTemperature(temperature_);
        
        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::oPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::wPhaseIdx, 0.0);
        fs.setSaturation(FluidSystem::gPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar densityL = 1e3;
        Scalar depth = maxDepth_ - pos[dim -1];
        Scalar pl = 1e5 - densityL*this->gravity()[dim - 1]*depth;

        PhaseVector pC;
        const auto &matParams = fineMaterialParams_;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(oPhaseIdx, pl + (pC[oPhaseIdx] - pC[wPhaseIdx]));
        fs.setPressure(wPhaseIdx, pl + (pC[wPhaseIdx] - pC[wPhaseIdx]));
        fs.setPressure(gPhaseIdx, pl + (pC[gPhaseIdx] - pC[wPhaseIdx]));
        
        // reset all mole fractions to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fs.setMoleFraction(phaseIdx, compIdx, 0.0);
        
        //////
        // set composition of the gas and water phases
        //////
        fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        fs.setMoleFraction(gPhaseIdx, gCompIdx, 1.0);

        //////
        // set composition of the oil phase
        //////
        
        // retrieve the relevant black-oil parameters from the fluid
        // system.
        Scalar Bo = FluidSystem::oilFormationVolumeFactor(fs.pressure(oPhaseIdx));
        Scalar Rs = FluidSystem::gasFormationFactor(fs.pressure(oPhaseIdx));
        Scalar rhoo = FluidSystem::surfaceDensity(oPhaseIdx)/Bo;
        Scalar rhogref = FluidSystem::surfaceDensity(gPhaseIdx);
        Scalar MG = FluidSystem::molarMass(gPhaseIdx);
        Scalar MO = FluidSystem::molarMass(oPhaseIdx);

        // calculate composition of oil phase in terms of mass
        // fractions.
        Scalar XoG = Rs*rhogref / rhoo;
        Scalar XoO = 1 - XoG;

        // convert to mole fractions
        Scalar avgMolarMass = MO*MG/(MG + XoO*(MO - MG));
        Scalar xoG = XoG*avgMolarMass/MG;
        Scalar xoO = 1 - XoG;
        
        // finally set the oil-phase composition. yeah!
        fs.setMoleFraction(oPhaseIdx, gCompIdx, xoG);
        fs.setMoleFraction(oPhaseIdx, oCompIdx, xoO);
    }

    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onRightBoundary_(pos) && (5 < pos[1]) && (pos[1] < 15); }

    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim-1] > layerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    BlackOilFluidState initialFluidState_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    std::string name_ ;
};
} //end namespace

#endif
