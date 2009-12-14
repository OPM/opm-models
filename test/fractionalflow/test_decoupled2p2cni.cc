// $Id$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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

#define DUNE_DEVEL_MODE
#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties_waterair.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "dumux/transport/problems/testproblem2p2cni.hh"
#include "dumux/transport/fv/decoupled2p2cni.hh"
#include <dumux/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // create a grid object
        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(20); N[0] = 20;
        FieldVector L(0);
        FieldVector H(300); H[0] = 300;
        Grid grid(N,L,H);

        grid.globalRefine(0);

        double tStart = 0;
        double tEnd = 2.5e5;
        int modulo = 1;
        double cFLFactor = 0.9;

        Dune::Liq_WaterAir wetmat;
        Dune::Gas_WaterAir nonwetmat;
        Dune::HomogeneousSoil<Grid, Scalar> soil;
        //    Dune::HeterogeneousSoil<Grid, Scalar> soil(grid, "permeab.dat", true);
        Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, wetmat, nonwetmat);

        Dune::VariableClass2p2cni<Grid,Scalar> var(grid);

        typedef Dune::TestProblem2p2cni<Grid, Scalar> TransProb;
        TransProb problem(var, wetmat, nonwetmat, soil, materialLaw, false);

        Dune::DiffusivePart<Grid, Scalar> diffPart;
        const Dune::Upwind<Scalar> numFl;

        typedef Dune::Decoupled2p2cni<Grid, Scalar> ModelType;
        ModelType model(grid, problem, grid.maxLevel(), diffPart, false, 0.8, numFl, "CG");

        Dune::ExplicitEulerStep<Grid, ModelType> timestep;
        Dune::TimeLoop<Grid, ModelType > timeloop(tStart, tEnd, "2p2cni", modulo, cFLFactor, 1e100, 1e100, timestep);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(model);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
