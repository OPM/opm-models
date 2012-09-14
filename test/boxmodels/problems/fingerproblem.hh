// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
 *   Copyright (C) 2012 by Markus Wolff                                      *
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
#ifndef DUMUX_FINGER_PROBLEM_HH
#define DUMUX_FINGER_PROBLEM_HH

#include "fingergridcreator.hh"

#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/parkerlenhard.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>

#include <dumux/boxmodels/immiscible/immiscibleproperties.hh>

#include <dune/common/fvector.hh>

#include <iostream>

#define FINGER_USE_PARKER_LENHARD 1

namespace Dumux
{

template <class TypeTag>
class FingerProblem;

//////////
// Specify the properties for the finger problem
//////////
namespace Properties
{
NEW_TYPE_TAG(FingerBaseProblem);

// declare the properties specific for the finger problem
NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(InitialWaterSaturation);

// set the GridCreator property
SET_TYPE_PROP(FingerBaseProblem, GridCreator, FingerGridCreator<TypeTag>);

// Retrieve the grid type from the grid creator
SET_TYPE_PROP(FingerBaseProblem, Grid, typename GET_PROP_TYPE(TypeTag, GridCreator)::Grid);

// Set the problem property
SET_TYPE_PROP(FingerBaseProblem, Problem, Dumux::FingerProblem<TypeTag>);

// Set the wetting phase
SET_PROP(FingerBaseProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(FingerBaseProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};

// Set the material Law
SET_PROP(FingerBaseProblem, MaterialLaw)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
#if FINGER_USE_PARKER_LENHARD
    // use the parker-lenhard hysteresis law
    typedef Dumux::ParkerLenhard<Scalar> TwoPMaterialLaw;
    typedef Dumux::ParkerLenhard<Scalar> ParkerLenhard;
#else
    // define the material law which is parameterized by effective
    // saturations
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> TwoPMaterialLaw;
#endif
    
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};

// Enable partial reassembly of the jacobian matrix?
//SET_BOOL_PROP(FingerBaseProblem, EnablePartialReassemble, true);

// Enable reuse of jacobian matrices?
//SET_BOOL_PROP(FingerBaseProblem, EnableJacobianRecycling, true);

// Write the solutions of individual newton iterations?
SET_BOOL_PROP(FingerBaseProblem, NewtonWriteConvergence, false);

// Use forward differences instead of central differences
SET_INT_PROP(FingerBaseProblem, NumericDifferenceMethod, +1);

// Enable smooth upwinding
SET_INT_PROP(FingerBaseProblem, EnableSmoothUpwinding, true);

// Enable constraints
SET_INT_PROP(FingerBaseProblem, EnableConstraints, true);

// Enable gravity
SET_BOOL_PROP(FingerBaseProblem, EnableGravity, true);

// define the properties specific for the finger problem
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeX, 0.1);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeY, 1.0);
SET_SCALAR_PROP(FingerBaseProblem, DomainSizeZ, 0.1);

SET_SCALAR_PROP(FingerBaseProblem, InitialWaterSaturation, 0.01);

SET_INT_PROP(FingerBaseProblem, CellsX, 20);
SET_INT_PROP(FingerBaseProblem, CellsY, 200);
SET_INT_PROP(FingerBaseProblem, CellsZ, 1);
}

/*!
 * \ingroup BoxTestProblems
 * \brief Two-phase problem might sometimes featuring a saturation overshoot.
 *
 * This problem was proposed by Rainer Helmig.
 */
template <class TypeTag>
class FingerProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
//!\cond 0
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;

    enum {
        // number of phases
        numPhases = FluidSystem::numPhases,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // equation indices
        contiWEqIdx = Indices::conti0EqIdx + wPhaseIdx,
        contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;

#if FINGER_USE_PARKER_LENHARD
    typedef typename GET_PROP(TypeTag, MaterialLaw)::ParkerLenhard ParkerLenhard;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
//!\endcond

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FingerProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 3e-6;
        FluidSystem::init();

        temperature_ = 273.15 + 20; // -> 20°C
    }

    /*!
     * \brief Initialize the content of the problem object.
     *
     * At this place there is the guarantee that all objects that do
     * not depend on the model have been allocated. The model itself
     * and the grid is allocated and initialized, although the initial
     * solution is not yet set.
    */
    void init()
    {
#if FINGER_USE_PARKER_LENHARD
        // parameters for the Van Genuchten law of the main imbibition
        // and the main drainage curves.
        micParams_.setVgAlpha(0.0037);
        micParams_.setVgN(4.7);

        mdcParams_.setVgAlpha(0.0037);
        mdcParams_.setVgN(4.7);

        // initialize the material parameter objects of the individual
        // finite volumes
        int n = GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView().size(dimWorld);
        materialParams_.resize(n);
        for (int i = 0; i < n; ++i) {
            materialParams_[i].setMicParams(&micParams_);
            materialParams_[i].setMdcParams(&mdcParams_);
            materialParams_[i].setSwr(0.02);
            materialParams_[i].setSnr(0.0);
            materialParams_[i].setSnre(0.0);
        }
#else
        // parameters for the Van Genuchten law
        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgN(4.7);

        materialParams_.setSwr(0.02);
        materialParams_.setSnr(0.00);
#endif

        K_ = this->toDimMatrix_(4.6e-10);

        setupInitialFluidState_();

        ParentType::init();
    }

    /*!
     * \brief Called after each successful time integration
     */
    void postTimeStep()
    {
#if FINGER_USE_PARKER_LENHARD
        // update the history of the hysteresis law
        ElementContext elemCtx(*this);

        auto elemIt = this->gridView().template begin<0>();
        const auto &elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            elemCtx.updateAll(*elemIt);
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                int globalIdx = elemCtx.globalSpaceIndex(scvIdx, /*timeIdx=*/0);
                const auto &fs = elemCtx.volVars(scvIdx, /*timeIdx=*/0).fluidState();
                ParkerLenhard::update(materialParams_[globalIdx],
                                      fs.saturation(wPhaseIdx));
            }
        }
#endif
    }

    /*!
     * \name Soil parameters
     */
    // \{

    /*!
     * \brief Intrinsic permability
     *
     * \return Intrinsic permeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    { return K_; }

    /*!
     * \brief Porosity
     *
     * \return Porosity
     */
     template <class Context>
     Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
     { return 0.4; }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \return the material parameters object
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
#if FINGER_USE_PARKER_LENHARD
        int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return materialParams_[globalSpaceIdx];
#else
        return materialParams_;
#endif
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
    { return "finger"; }

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
        const GlobalPosition &pos = context.cvCenter(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos) || onRightBoundary_(pos) || onLowerBoundary_(pos)) {
            values.setNoFlow();
        }
        else {
            assert(onUpperBoundary_(pos));

            values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidState_);
        };

        // override the value for the liquid phase by forced
        // imbibition of water on inlet boundary segments
        if (onInlet_(pos)) { 
            values[contiWEqIdx] = - 0.001; // [kg/(m^2 s)]
        }
    }

    template <class Context>
    void constraints(Constraints &values,
                     const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos) && !onInlet_(pos)) {
            values.setAllConstraint();
            values.assignNaive(initialFluidState_);
        }
        else if (onLowerBoundary_(pos)) {
            values.setAllConstraint();
            values.assignNaive(initialFluidState_);
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
     * \param pos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        // assign the primary variables
        values.assignNaive(initialFluidState_);
    }


    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { values = 0; }
    // \}

private:   
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - pos[0])/width;
        
        if (!onUpperBoundary_(pos))
            return false;
        
        Scalar xInject[] = { 0.25, 0.75 };
        Scalar injectLen[] = { 0.1, 0.1 };
        for (unsigned i = 0; i < sizeof(xInject)/sizeof(Scalar); ++ i) {
            if (xInject[i] - injectLen[i]/2 < lambda &&  lambda < xInject[i] + injectLen[i]/2)
                return true;
        }
        return false;
    }

    void setupInitialFluidState_()
    {
        auto &fs = initialFluidState_;
        fs.setPressure(wPhaseIdx, /*pressure=*/1e5);

        Scalar Sw = GET_PARAM(TypeTag, Scalar, InitialWaterSaturation);
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1 - Sw);

        fs.setTemperature(temperature_);

        // set the absolute pressures
        Scalar pn = 1e5;
        fs.setPressure(nPhaseIdx, pn);
        fs.setPressure(wPhaseIdx, pn);
    }

    DimMatrix K_;

#if FINGER_USE_PARKER_LENHARD
    typename MaterialLawParams::VanGenuchtenParams micParams_;
    typename MaterialLawParams::VanGenuchtenParams mdcParams_;

    std::vector<MaterialLawParams> materialParams_;
#else
    MaterialLawParams materialParams_;
#endif
    ImmiscibleFluidState<Scalar, FluidSystem> initialFluidState_;

    Scalar temperature_;
    Scalar eps_;
};

} //end namespace

#endif
