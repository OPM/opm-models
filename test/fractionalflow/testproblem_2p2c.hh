// $Id: testproblem_2p2c.hh 2172 2009-06-19 11:57:34Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Jochen Fritz                                 *
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

#ifndef TESTPROBLEM_2P2C_HH
#define TESTPROBLEM_2P2C_HH

#include "dumux/transport/decoupled2p2cproblem.hh"

namespace Dumux
{

//! Example problem class for decoupled 2p2c simulations
template<class GridView, class Scalar>
class Testproblem_2p2c
    : public DecoupledProblem2p2c<GridView, Scalar>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum {dim=GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Entity;

public:

    virtual const Dune::FieldVector<Scalar,dim> gravity()
    {
        Dune::FieldVector<Scalar,dim> gravity_(0);
        gravity_[2] = -10;
        return gravity_;
    }

    Scalar temperature(const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                           const Dune::FieldVector<Scalar, dim>& localPos, double t)
    {
        return 283.15;
    }

    BoundaryConditions2p2c::Flags bc_type (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                           const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {

        return BoundaryConditions2p2c::concentration;
    }

    BoundaryConditions2p2c::Flags initcond_type (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                                 const Dune::FieldVector<Scalar, dim>& localPos) const
    {
        return BoundaryConditions2p2c::concentration;
    }

    BoundaryConditions::Flags press_bc_type (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                             const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
            return Dumux::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dumux::BoundaryConditions::neumann;
    }


    Scalar dirichlet (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element, const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        return (globalPos[0] < 1e-6) ? (2.5e5 - 10000 * globalPos[2]) : (2e5 - 10000 * globalPos[2]);
    }

    Scalar dirichletConcentration (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                   const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        return 1;
    }

    Scalar dirichletSat (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                         const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        if (globalPos[0] < 15)
            return 0;
        else
            return 0;
    }

    virtual Dune::FieldVector<Scalar,2> neumann (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                           const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        Dune::FieldVector<Scalar,2> J_(0);
        return J_;
    }

    virtual Dune::FieldVector<Scalar,2> source (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                                          const Dune::FieldVector<Scalar, dim>& localPos, double t) const
    {
        Dune::FieldVector<Scalar,2> q_(0);
        if (fabs(globalPos[0] - 4.5) < 1 && fabs(globalPos[1] - 4.5) < 1) q_[1] = 0.0001;
        return q_;
    }

    Scalar initSat (const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                    const Dune::FieldVector<Scalar, dim>& localPos) const
    {
        return 0.999;
    }

    Scalar initConcentration(const Dune::FieldVector<Scalar, dim>& globalPos, const Entity& element,
                             const Dune::FieldVector<Scalar, dim>& localPos) const
    {
        return 1;
    }

    Testproblem_2p2c(GridView& gv, Dumux::VariableClass2p2c<GridView, Scalar>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<typename GridView::Grid, Scalar>& s,
                     int level, TwoPhaseRelations<typename GridView::Grid, Scalar>& law = *(new TwoPhaseRelations<typename GridView::Grid, Scalar>),const bool cap = false)
        : DecoupledProblem2p2c<GridView, Scalar>(var, liq, gas, s, law, cap), gridview(gv)
    {
    }

private:
    GridView& gridview;
};

}
#endif
