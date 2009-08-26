#ifdef HAVE_DUNE_PDELAB

#ifndef DUMUX_ASSEMBLERPDELAB_HH
#define DUMUX_ASSEMBLERPDELAB_HH

#include<dune/pdelab/finiteelementmap/p1fem.hh>
#include<dune/pdelab/finiteelementmap/q1fem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include"boundarytypespdelab.hh"
#include"boxjacobianpdelab.hh"

namespace Dune {
namespace PDELab {

template<typename E, typename BTypes, typename GridFunctionSpace>
class BoxConstraintsTransformation: public ConstraintsTransformation<unsigned int, E>
{
public:
	BoxConstraintsTransformation(BTypes& bTypes, GridFunctionSpace& gridFunctionSpace)
	{
		constraints(bTypes, gridFunctionSpace, *this);
	}

	BoxConstraintsTransformation()
	{}
};


template<typename GV, typename LFEM, typename CE=NoConstraints,
         typename B=StdVectorBackend, typename P=GridFunctionGeneralMapper>
class BoxGridFunctionSpace : public GridFunctionSpace<GV,LFEM,CE,B,P>
{
public:
	typedef GridFunctionSpace<GV,LFEM,CE,B,P> BaseType;

	BoxGridFunctionSpace (const GV& gridView, const LFEM& lfem, const CE& ce_=CE())
	:  BaseType(gridView, lfem, ce_)
	{}
};

template<typename GV, typename LFEM, typename CE, typename B, typename IIS>
class BoxGridFunctionSpace<GV,LFEM,CE,B,GridFunctionStaticSize<IIS> >
: public GridFunctionSpace<GV,LFEM,CE,B,GridFunctionStaticSize<IIS> >
{
public:
	typedef GridFunctionSpace<GV,LFEM,CE,B,GridFunctionStaticSize<IIS> > BaseType;
	typedef typename BaseType::template VectorContainer<double>::Type Vector;
	typedef GhostDataHandle<BaseType, Vector> GDH;

	BoxGridFunctionSpace (const GV& gridView, const LFEM& lfem, const CE& ce_, std::vector<int>& intghost)
	:  BaseType(gridView, lfem, ce_)
	{
		Vector ghost(*this, 0.0);
		GDH ghostDataHandle(*this, ghost);
		if (gridView.comm().size() > 1)
			gridView.communicate(ghostDataHandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
		ghost.std_copy_to(intghost); // copy to std::vector as Constraints object is allocated before
	}
};
}
}

template<class TypeTag>
class AssemblerPDELab
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    enum{dim = GridView::dimension};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> FEM;
//    typedef typename Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
    typedef typename Dune::PDELab::OverlappingConformingDirichletConstraints CON;
    typedef Dune::PDELab::BoxGridFunctionSpace<GridView, FEM, CON,
      Dune::PDELab::ISTLVectorBackend<numEq>, Dune::PDELab::SimpleGridFunctionStaticSize> ScalarGridFunctionSpace;
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace,numEq,Dune::PDELab::GridFunctionSpaceBlockwiseMapper> GridFunctionSpace;
	typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type Vector;
	typedef Dune::PDELab::GhostDataHandle<GridFunctionSpace, Vector> GhostDataHandle;
    typedef BoundaryTypesPDELab<TypeTag, Problem> BTypes;
    typedef Dune::PDELab::BoxConstraintsTransformation<Scalar, BTypes, GridFunctionSpace> ConstraintsTrafo;
    typedef BoxJacobianPDELab<TypeTag> LocalOperator;
    typedef Dune::PDELab::GridOperatorSpace<GridFunctionSpace,
                                            GridFunctionSpace,
                                            LocalOperator,
                                            ConstraintsTrafo,
                                            ConstraintsTrafo,
                                            Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
                                            true
                                           > GridOperatorSpace;
    typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef Matrix RepresentationType;

    AssemblerPDELab(Problem& problem)
    : problem_(problem), /*cn_(intghost_), */scalarGridFunctionSpace_(problem.gridView(), fem_, cn_, intghost_),
    gridFunctionSpace_(scalarGridFunctionSpace_/*, scalarGridFunctionSpace_*/), bTypes_(problem),
    constraintsTrafo_(bTypes_, gridFunctionSpace_), localOperator_(problem),
    gridOperatorSpace_(gridFunctionSpace_, constraintsTrafo_, gridFunctionSpace_, constraintsTrafo_, localOperator_),
    matrix_(gridOperatorSpace_)
    {
        matrix_ = 0;
    }

    //! return const reference to matrix
    const Matrix& operator* () const
    {
        return matrix_;
    }

    //! return reference to matrix
    Matrix& operator* ()
    {
        return matrix_;
    }

    void assemble(LocalJacobian& loc, SolutionFunction& u, SolutionFunction& f)
    {
    	matrix_ = 0;
    	gridOperatorSpace_.jacobian(*u, matrix_);
    	*f = 0;
    	gridOperatorSpace_.residual(*u, *f);
		set_constrained_dofs(constraintsTrafo_, 0.0, *f);
		set_constrained_dofs(constraintsTrafo_, 0.0, *u);
    }

    const GridFunctionSpace& gridFunctionSpace() const
    {
        return gridFunctionSpace_;
    }

    const ConstraintsTrafo& constraintsTrafo() const
    {
        return constraintsTrafo_;
    }

private:
	Problem& problem_;
	std::vector<int> intghost_;
	CON cn_;
    FEM fem_;
    ScalarGridFunctionSpace scalarGridFunctionSpace_;
    GridFunctionSpace gridFunctionSpace_;
    BTypes bTypes_;
    ConstraintsTrafo constraintsTrafo_;
    LocalOperator localOperator_;
    GridOperatorSpace gridOperatorSpace_;
    Matrix matrix_;
};

#endif

#endif
