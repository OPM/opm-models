// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/fluids/uniform.hh"
#include "dumux/transport/fv/implicitfvtransport.hh"
#include "dumux/transport/fv/computeupwind.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "simplenonlinearproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 4e9;
        double dt = 5e7;
	double firstDt = dt;
	double maxDt = dt;
        int modulo = 1;

        // create a grid object
        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Grid::LevelGridView GridView;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(1); N[0] = 16;
        FieldVector L(0);
        FieldVector H(300); H[0] = 600;
        Grid grid(N,L,H);

        grid.globalRefine(0);

        GridView gridView(grid.levelView(0));

        Dune::Uniform mat(0.2);
        //Dune::HomogeneousLinearSoil<Grid, Scalar> soil;
        Dune::HomogeneousNonlinearSoil<Grid, Scalar> soil;
        Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridView, Scalar> VariableClass;

        double initsat=0;
        Dune::FieldVector<double,dim> vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VariableClass variables(gridView,initsat,vel);

        //Dune::SimpleProblem<GridView, Scalar, VariableClass> ProblemType;
        typedef Dune::SimpleNonlinearProblem<GridView, Scalar, VariableClass> ProblemType;
	ProblemType problem(variables, mat, mat , soil, materialLaw,L,H);

        //Dune::DiffusivePart<GridView, Scalar> diffusivePart;
        Dune::CapillaryDiffusion<GridView, Scalar, VariableClass, ProblemType> capillaryDiffusion(problem, soil, false);
        Dune::ComputeUpwind<GridView, Scalar, VariableClass> computeNumFlux(problem);

        typedef Dune::ImplicitFVTransport<GridView, Scalar, VariableClass> Transport;

        //Transport transport(gridView, problem, dt, diffusivePart, computeNumFlux);
        Transport transport(gridView, problem, capillaryDiffusion, computeNumFlux);

	Dune::RungeKuttaStep<Grid, Transport> timeStep(1);
	//Dune::ImplicitEulerStep<Grid, Transport> timeStep;
        Dune::TimeLoop<Grid, Transport > timeloop(tStart, tEnd, dt, "implicitfvtransport", modulo, maxDt, firstDt, timeStep);

        timeloop.execute(transport);

        printvector(std::cout, variables.saturation(), "saturation", "row", 200, 1);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
