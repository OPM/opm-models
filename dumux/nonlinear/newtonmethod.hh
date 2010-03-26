// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
#ifndef DUMUX_NEWTONMETHOD_HH
#define DUMUX_NEWTONMETHOD_HH

#ifdef DUMUX_NEWTONMETHOD_DEPRECATED_HH
# error "Never use the old and the new newton method in one program!"
#endif

#include <limits>
#include <dumux/common/exceptions.hh>

#include <dumux/io/vtkmultiwriter.hh>


namespace Dumux
{
/*!
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template <class TypeTag>
class NewtonMethod
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::JacobianAssembler       JacobianAssembler;
public:
    NewtonMethod(Model &model)
    {
        model_ = NULL;
        uOld = NULL;
        f = NULL;
    }

    ~NewtonMethod()
    { }

    /*!
     * \brief Returns a reference to the current numeric model.
     */
    Model &model()
    { return *model_; }

    /*!
     * \brief Returns a reference to the current numeric model.
     */
    const Model &model() const
    { return *model_; }

    /*!
     * \brief Returns true iff the newton method should be verbose
     *        about what it is going on;
     */
    bool verbose() const
    { return (model().gridView().comm().rank() == 0); }

    /*!
     * \brief Run the newton method. The controller is responsible
     *        for all the strategic decisions.
     */
    template <class NewtonController>
    bool execute(Model &model, NewtonController &ctl)
    {
        try {
            return execute_(model, ctl);
        }
        catch (const Dune::ISTLError &e) {
            if (verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            model_ = NULL;
            return false;
        }
        catch (const Dumux::NumericalProblem &e) {
            if (verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            model_ = NULL;
            return false;
        };
    }

protected:
    template <class NewtonController>
    bool execute_(Model &model, NewtonController &ctl)
    {
    	if (!uOld)
    	{
#ifdef HAVE_DUNE_PDELAB
            uOld = new SolutionFunction(model.jacobianAssembler().gridFunctionSpace(), 0.0);
            f = new SolutionFunction(model.jacobianAssembler().gridFunctionSpace(), 0.0);
#else
            uOld = new SolutionFunction(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0);
            f = new SolutionFunction(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0);
#endif
    	}
        model_ = &model;

        // TODO (?): u shouldn't be hard coded to the model
        SolutionFunction  &u             = model.curSolFunction();
        LocalJacobian     &localJacobian = model.localJacobian();
        JacobianAssembler &jacobianAsm   = model.jacobianAssembler();
        
        // tell the controller that we begin solving
        ctl.newtonBegin(this, u);
        
        // execute the method as long as the controller thinks
        // that we should do another iteration
        while (ctl.newtonProceed(u))
        {
            // notify the controller that we're about to start
            // a new timestep
            ctl.newtonBeginStep();

            // make the current solution to the old one
            *(*uOld) = *u;
            *(*f) = 0;

            if (verbose()) {
                std::cout << "Assembling global jacobian";
                std::cout.flush();
            }
            // linearize the problem at the current solution
            jacobianAsm.assemble(localJacobian, u, *f);

            // solve the resultuing linear equation system
            if (verbose()) {
                std::cout << "\rSolve Mx = r              ";
                std::cout.flush();
            }

/*
            int i = 38;
            std::cout << (boost::format("*jacobianAsm[%d][%d]: \n")%i%i).str()
                      << (*jacobianAsm)[i][i];
            std::cout << "*u["<<i<<"]: "
                      << (*u)[i] << "\n";
            exit(1);
*/

/*
#if HAVE_DUNE_PDELAB
            printmatrix(std::cout, (*jacobianAsm).base(), "J PDELab", "row");
#else
            printmatrix(std::cout, *jacobianAsm, "J", "row");
#endif
            std::cout << "rhs:" << *(*f) << "\n";
            exit(1);
*/

            // ask the controller to solve the linearized system
            *u = 0; // set the delta vector to zero before solving the
                    // linear system!
            ctl.newtonSolveLinear(*jacobianAsm, u, *(*f));

            // update the current solution (i.e. uOld) with the delta
            // (i.e. u). The result is stored in u
            ctl.newtonUpdate(u, *uOld);

            // tell the controller that we're done with this iteration
            ctl.newtonEndStep(u, *uOld);
        }

        // tell the controller that we're done
        ctl.newtonEnd();
        
        if (!ctl.newtonConverged()) {
            ctl.newtonFail();
            model_ = NULL;
            return false;
        }
        
        ctl.newtonSucceed();
        model_ = NULL;
        return true;
    }


private:
    SolutionFunction *uOld;
    SolutionFunction *f;

    Model        *model_;
    bool          verbose_;
};
}

#endif
