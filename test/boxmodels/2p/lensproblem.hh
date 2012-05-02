// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */
#ifndef DUMUX_LENS_PROBLEM_HH
#define DUMUX_LENS_PROBLEM_HH

#include "lensgridcreator.hh"

#include <dumux/boxmodels/2p/2pmodel.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>

#include <dune/common/fvector.hh>

#include <iostream>

namespace Dumux
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(BoxTwoP));

// declare the properties specific for the lens problem
NEW_PROP_TAG(LensLowerLeftX);
NEW_PROP_TAG(LensLowerLeftY);
NEW_PROP_TAG(LensLowerLeftZ);
NEW_PROP_TAG(LensUpperRightX);
NEW_PROP_TAG(LensUpperRightY);
NEW_PROP_TAG(LensUpperRightZ);

NEW_PROP_TAG(GridSizeX);
NEW_PROP_TAG(GridSizeY);
NEW_PROP_TAG(GridSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

// set the GridCreator property
SET_TYPE_PROP(LensProblem, GridCreator, LensGridCreator<TypeTag>);

// Retrieve the grid type from the grid creator
SET_TYPE_PROP(LensProblem, Grid, typename GET_PROP_TYPE(TypeTag, GridCreator)::Grid);

// Set the problem property
SET_TYPE_PROP(LensProblem, Problem, Dumux::LensProblem<TypeTag>);

#if 1
// Set the wetting phase
SET_PROP(LensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(LensProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleDNAPL<Scalar> > type;
};
#else
// OR: set the fluid system
SET_TYPE_PROP(LensProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar)>);
#endif

// Set the material Law
SET_PROP(LensProblem, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> TwoPMaterialLaw;
    
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

public:
    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};

// Enable partial reassembly of the jacobian matrix?
SET_BOOL_PROP(LensProblem, EnablePartialReassemble, true);

// Enable reuse of jacobian matrices?
SET_BOOL_PROP(LensProblem, EnableJacobianRecycling, true);

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(LensProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(LensProblem, NumericDifferenceMethod, +1);

// Linear solver settings
SET_TYPE_PROP(LensProblem, LinearSolver, Dumux::BoxBiCGStabILU0Solver<TypeTag> );
SET_INT_PROP(LensProblem, LinearSolverVerbosity, 0);
SET_INT_PROP(LensProblem, PreconditionerIterations, 1);
SET_SCALAR_PROP(LensProblem, PreconditionerRelaxation, 1.0);

// Enable gravity
SET_BOOL_PROP(LensProblem, EnableGravity, true);

// define the properties specific for the lens problem
SET_SCALAR_PROP(LensProblem, LensLowerLeftX, 1.0);
SET_SCALAR_PROP(LensProblem, LensLowerLeftY, 2.0);
SET_SCALAR_PROP(LensProblem, LensLowerLeftZ, 0.0);
SET_SCALAR_PROP(LensProblem, LensUpperRightX, 4.0);
SET_SCALAR_PROP(LensProblem, LensUpperRightY, 3.0);
SET_SCALAR_PROP(LensProblem, LensUpperRightZ, 1.0);

SET_SCALAR_PROP(LensProblem, GridSizeX, 6.0);
SET_SCALAR_PROP(LensProblem, GridSizeY, 4.0);
SET_SCALAR_PROP(LensProblem, GridSizeZ, 1.0);

SET_INT_PROP(LensProblem, CellsX, 48);
SET_INT_PROP(LensProblem, CellsY, 32);
SET_INT_PROP(LensProblem, CellsZ, 16);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p -parameterFile test_2p.input</tt>
 */
template <class TypeTag>
class LensProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // equation indices
        contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx,
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    LensProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 3e-6;
        FluidSystem::init();


        temperature_ = 273.15 + 20; // -> 20°C
        lensLowerLeft_[0] = GET_PARAM(TypeTag, Scalar, LensLowerLeftX);
        lensLowerLeft_[1] = GET_PARAM(TypeTag, Scalar, LensLowerLeftY);
        lensUpperRight_[0] = GET_PARAM(TypeTag, Scalar, LensUpperRightX);
        lensUpperRight_[1] = GET_PARAM(TypeTag, Scalar, LensUpperRightY);
        
        if (dimWorld == 3) {
            lensLowerLeft_[2] = GET_PARAM(TypeTag, Scalar, LensLowerLeftZ);
            lensUpperRight_[2] = GET_PARAM(TypeTag, Scalar, LensUpperRightZ);
        }

        // parameters for the Van Genuchten law
        // alpha and n
        lensMaterialParams_.setVgAlpha(0.00045);
        lensMaterialParams_.setVgN(7.3);
        outerMaterialParams_.setVgAlpha(0.0037);
        outerMaterialParams_.setVgN(4.7);

        lensK_ = this->toTensor_(9.05e-12);
        outerK_ = this->toTensor_(4.6e-10);
        
        this->gravity_ = 0;
        this->gravity_[1] = -9.81; // [m/s^2]
    }

    /*!
     * \name Soil parameters
     */
    // \{

    /*!
     * \brief Intrinsic permability
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Intrinsic permeability
     */
    template <class Context>
    const Tensor &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \brief Porosity
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return Porosity
     */
     template <class Context>
     Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
     { return 0.4; }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param element The current element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \brief Returns the temperature within the domain.
     */
    template <class Context>
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return temperature_; }

    // \}

    /*!
     * \name Auxiliary methods
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "lens"; }

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: " << storage << std::endl;
        }
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
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            // free flow boundary
            Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/1e5);

            Scalar T = temperature(context, spaceIdx, timeIdx);
            Scalar pw, Sw;
        
            // set wetting phase pressure and saturation
            if (onLeftBoundary_(globalPos))
            {
                Scalar height = this->bboxMax()[1] - this->bboxMin()[1];
                Scalar depth = this->bboxMax()[1] - globalPos[1];
                Scalar alpha = (1 + 1.5/height);

                // hydrostatic pressure scaled by alpha
                pw = 1e5 - alpha*densityW*this->gravity()[1]*depth;
                Sw = 1.0;
            }
            else
            {
                // if this triggers, something went wrong in boundaryTypes()!
                assert(onRightBoundary_(globalPos));

                Scalar depth = this->bboxMax()[1] - globalPos[1];

                // hydrostatic pressure
                pw = 1e5 - densityW*this->gravity()[1]*depth;
                Sw = 1.0;
            }
        
            // specify a full fluid state using pw and Sw
            const MaterialLawParams &matParams =
                this->materialLawParams(context, spaceIdx, timeIdx);
            
            Dumux::ImmiscibleFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
            fs.setSaturation(wPhaseIdx, Sw);
            fs.setSaturation(nPhaseIdx, 1 - Sw);
            fs.setTemperature(T);

            Scalar pC[numPhases];
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            fs.setPressure(wPhaseIdx, pw);
            fs.setPressure(nPhaseIdx, pw +  pC[nPhaseIdx] - pC[wPhaseIdx]);

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(globalPos)) {
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiNEqIdx] = -0.04; // kg / (m^2 * s)

            // impose a forced flow boundary
            values.setMassRate(massRate);
        }
        else {
            // no flow boundary
            values.setNoFlow();
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
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        Scalar depth = this->bboxMax()[1] - globalPos[1];
        
        ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setPressure(wPhaseIdx, /*pressure=*/1e5);
        fs.setTemperature(temperature_);
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fs, wPhaseIdx);
        Scalar densityW = FluidSystem::density(fs, paramCache, wPhaseIdx);

        // set the fluid state's temperature
        Scalar T = temperature(context, spaceIdx, timeIdx);
        
        // hydrostatic pressure
        Scalar pw = 1e5 - densityW*this->gravity()[1]*depth;
        Scalar Sw = 1.0;

        // assign primary variables
        const MaterialLawParams &matParams =
            this->materialLawParams(context, spaceIdx, timeIdx);
        values.assignImmiscibleFromWetting(T, pw, Sw, matParams);
    }


    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { values = 0; }
    // \}

private:   
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bboxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bboxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bboxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bboxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Tensor lensK_;
    Tensor outerK_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;

    Scalar temperature_;
    Scalar eps_;
};
} //end namespace

#endif
