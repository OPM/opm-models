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
#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

#ifdef DUNE_NEWTONMETHOD_DEPRECATED_HH
# error "Never use the old and the new newton method in one program!"
#endif

#include <limits>
#include <dumux/exceptions.hh>

#include <dumux/io/vtkmultiwriter.hh>


namespace Dune
{
/*!
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template <class ModelT>
class NewtonMethod
{
public:
    typedef ModelT Model;

    //    private:
    typedef typename Model::NewtonTraits::Function          Function;
    typedef typename Model::NewtonTraits::LocalJacobian     LocalJacobian;
    typedef typename Model::NewtonTraits::JacobianAssembler JacobianAssembler;
    typedef typename Model::NewtonTraits::Scalar            Scalar;
    typedef NewtonMethod<Model>                             ThisType;

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
        catch (const Dune::NumericalProblem &e) {
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
            uOld = new Function(model.jacobianAssembler().gridFunctionSpace(), 0.0);
            f = new Function(model.jacobianAssembler().gridFunctionSpace(), 0.0);
#else
            uOld = new Function(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0);
            f = new Function(model.gridView(), model.gridView(), model.gridView().overlapSize(0) == 0);
#endif
    	}
        model_ = &model;

        // TODO (?): u shouldn't be hard coded to the model
        Function          &u             = model.curSolFunction();
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
    Function       *uOld;
    Function       *f;

    Model        *model_;
    bool          verbose_;
};
}

#endif
