// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#ifndef DUMUX_BOX_PROPERTIES_HH
#define DUMUX_BOX_PROPERTIES_HH

#include <dumux/auxiliary/properties.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dune
{
namespace Properties
{
/*!
 * \addtogroup BoxScheme
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the box-scheme
NEW_TYPE_TAG(BoxScheme);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//!< Property tag for scalar values
NEW_PROP_TAG(Scalar);

//! Property tag for types associated with the solution of the PDE.
//! This means vectors of primary variables, solution functions on the
//! grid, and elements, and shape functions.
NEW_PROP_TAG(SolutionTypes);

NEW_PROP_TAG(Grid);     //!< The type of the DUNE grid
NEW_PROP_TAG(GridView); //!< The type of the grid view

NEW_PROP_TAG(ReferenceElements); //!< DUNE reference elements to be used
NEW_PROP_TAG(FVElementGeometry); //! The type of the finite-volume geometry in the box scheme

NEW_PROP_TAG(Problem);       //!< The type of the problem
NEW_PROP_TAG(Model);         //!< The type of the discretization
NEW_PROP_TAG(NumEq);         //!< Number of equations in the system of PDEs
NEW_PROP_TAG(LocalJacobian); //!< The type of the local jacobian operator

#if HAVE_DUNE_PDELAB
NEW_PROP_TAG(LocalFEMSpace); //!< The local finite element space used for the finite element interpolation
NEW_PROP_TAG(PDELabTypes); //!< The types used from dune-pdelab
#endif

NEW_PROP_TAG(VertexData);  //!< Data merging from constitutive relations defined on the vertices of the grid
NEW_PROP_TAG(ElementData); //!< Data merging from constitutive relations defined on the elements of the grid
NEW_PROP_TAG(FluxData);    //!< Data required to calculate a flux over a face

NEW_PROP_TAG(NewtonMethod);     //!< The type of the newton method
NEW_PROP_TAG(NewtonController); //!< The type of the newton controller
}
}

#include <dumux/fvgeometry/fvelementgeometry.hh>

#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#if HAVE_DUNE_PDELAB
#include <dumux/boxmodels/pdelab/assemblerpdelab.hh>
#include <dumux/boxmodels/pdelab/functionpdelab.hh>
#include <dumux/boxmodels/pdelab/boxdirichletconstraints.hh>
#else
#include <dune/disc/operators/p1operator.hh>
#endif

#include <dumux/boundarytypes.hh>

namespace Dune {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for scalar values to double
SET_PROP_DEFAULT(Scalar)
{ typedef double   type; };

//! Use the leaf grid view if not defined otherwise
SET_PROP_DEFAULT(GridView)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef typename Grid::LeafGridView type;
};

/*!
 * \brief Specify the reference elements which we ought to use.
 *
 * We use Dune::ReferenceElements by default (-> old entity
 * numbering).
 *
 * TODO: Some specialization if the grid only supports one kind of
 *       cells would be nice. this would be better fixed inside DUNE,
 *       though. something like:
 *       GenericReferenceElements<GeometryType<cube, dim> >
 */
SET_PROP_DEFAULT(ReferenceElements)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

    typedef typename Grid::ctype CoordScalar;
    static const int dim = Grid::dimension;

public:
#if HAVE_DUNE_PDELAB
    typedef Dune::GenericReferenceElements<CoordScalar, dim> Container;
    typedef Dune::GenericReferenceElement<CoordScalar, dim>  ReferenceElement;
#else
    typedef Dune::ReferenceElements<CoordScalar, dim> Container;
    typedef Dune::ReferenceElement<CoordScalar, dim>  ReferenceElement;
#endif
};

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Set the default for the FVElementGeometry
SET_PROP(BoxScheme, FVElementGeometry)
{
#if HAVE_DUNE_PDELAB
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace))  LocalFEMSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;

public:
    typedef Dune::FVElementGeometry<Grid, LocalFEMSpace>  type;

#else
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
public:
    typedef Dune::FVElementGeometry<Grid>  type;
#endif
};

//! use the plain newton method for the box scheme by default
SET_PROP(BoxScheme, NewtonMethod)
{public:
    typedef Dune::NewtonMethod<TypeTag> type;
};

//! use the plain newton controller for the box scheme by default
SET_PROP(BoxScheme, NewtonController)
{public:
    typedef Dune::NewtonController<TypeTag> type;
};

/*!
 * \brief Specifies the types which are assoicated with a solution.
 *
 * This means shape functions, solution vectors, etc.
 */
SET_PROP(BoxScheme, SolutionTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))    Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))  GridView;
    typedef typename GridView::Grid                          Grid;
    typedef typename Grid::ctype                             CoordScalar;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))
    };

    template<int dim>
    struct VertexLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == 0; }
    };

    template<int dim>
    struct ElementLayout {
        bool contains (Dune::GeometryType gt) const
        { return gt.dim() == dim; }
    };

public:
    //! A solution function. This is a function with the same domain as the grid.
#ifdef HAVE_DUNE_PDELAB
	typedef FunctionPDELab<TypeTag>               SolutionFunction;
#else
    typedef Dune::P1Function<GridView, Scalar, GridView, numEq>       SolutionFunction;
#endif

    /*!
     * \brief The type which maps an entity with an attached degree of
     *        freedom to an index in the solution.
     */
    typedef typename SolutionFunction::VM                              DofEntityMapper;

    /*!
     * \brief Mapper for the grid view's vertices.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout> VertexMapper;

    /*!
     * \brief Mapper for the grid view's elements.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, ElementLayout> ElementMapper;

    /*!
     * \brief The type of a solution at a fixed time.
     *
     * This is the representation of the solution function and defines
     * a primary variable vector at each degree of freedom.
     */
    typedef typename SolutionFunction::RepresentationType             Solution;
    /*!
     * \brief A vector of primary variables.
     */
    typedef typename SolutionFunction::BlockType                      PrimaryVarVector;

    /*!
     * \brief The solution for a single finite element.
     */
    typedef Dune::BlockVector<PrimaryVarVector>                       SolutionOnElement;

    /*!
     * \brief Vector of boundary types at a degree of freedom.
     */
    typedef Dune::BoundaryTypes<numEq>  BoundaryTypeVector;

    /*!
     * \brief Assembler for the global jacobian matrix.
     */
#ifdef HAVE_DUNE_PDELAB
	typedef AssemblerPDELab<TypeTag>      JacobianAssembler;
#else
    typedef Dune::P1OperatorAssembler<Grid, Scalar, GridView, GridView, numEq>  JacobianAssembler;
#endif
};



#ifdef HAVE_DUNE_PDELAB
//! use the local FEM space associated with cubes by default
SET_PROP(BoxScheme, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

SET_PROP(BoxScheme, PDELabTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
public:
    //typedef typename Dune::PDELab::NonoverlappingConformingDirichletConstraints Constraints;
    typedef typename Dune::NonoverlappingBoxDirichletConstraints<TypeTag> Constraints;
    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints,
							Dune::PDELab::ISTLVectorBackend<numEq> > ScalarGridFunctionSpace;
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq,
							Dune::PDELab::GridFunctionSpaceBlockwiseMapper> GridFunctionSpace;
	typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type ConstraintsTrafo;
    typedef BoxJacobianPDELab<TypeTag> LocalOperator;
    typedef Dune::PDELab::GridOperatorSpace<GridFunctionSpace,
                                            GridFunctionSpace,
                                            LocalOperator,
                                            ConstraintsTrafo,
                                            ConstraintsTrafo,
                                            Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
                                            true
                                           > GridOperatorSpace;
};
#endif


// \}

}
}

#endif
