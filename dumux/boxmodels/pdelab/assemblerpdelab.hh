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
    typedef typename Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
//    typedef typename Dune::PDELab::OverlappingConformingDirichletConstraints CON;
    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, CON,
      Dune::PDELab::ISTLVectorBackend<numEq>/*, Dune::PDELab::SimpleGridFunctionStaticSize*/> ScalarGridFunctionSpace;
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace,numEq,Dune::PDELab::GridFunctionSpaceBlockwiseMapper> GridFunctionSpace;
	typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type Vector;
    typedef BoundaryTypesPDELab<TypeTag> BTypes;
//    typedef Dune::PDELab::PowerGridFunction<ScalarBTypes, numEq> BTypes;
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
    typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef Matrix RepresentationType;

    AssemblerPDELab(Problem& problem)
    : problem_(problem)
    {
    	fem_ = new FEM();
    	cn_ = new CON(intghost_);
    	scalarGridFunctionSpace_ = new ScalarGridFunctionSpace(problem_.gridView(), *fem_, *cn_);
    	gridFunctionSpace_ = new GridFunctionSpace(*scalarGridFunctionSpace_);

    	Vector ghost(*gridFunctionSpace_, 0.0);
    	Dune::PDELab::GhostDataHandle<GridFunctionSpace, Vector> ghostDataHandle(*gridFunctionSpace_, ghost);
    	if (problem_.gridView().comm().size() > 1)
    		problem_.gridView().communicate(ghostDataHandle,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
    	ghost.std_copy_to(intghost_); // copy to std::vector as Constraints object is allocated before

//    	scalarBTypes0_ = new ScalarBTypes(problem_, 0);
//    	scalarBTypes1_ = new ScalarBTypes(problem_, 1);
    	bTypes_ = new BTypes(problem_);//*scalarBTypes0_, *scalarBTypes1_);
    	constraintsTrafo_ = new ConstraintsTrafo();
    	Dune::PDELab::constraints(*bTypes_, *gridFunctionSpace_, *constraintsTrafo_);

    	localOperator_ = new LocalOperator(problem_);
    	gridOperatorSpace_ = new GridOperatorSpace(*gridFunctionSpace_, *constraintsTrafo_,
												   *gridFunctionSpace_, *constraintsTrafo_, *localOperator_);

    	matrix_ = new Matrix(*gridOperatorSpace_);
    	*matrix_ = 0;
    }

    //! return const reference to matrix
    const Matrix& operator* () const
    {
        return *matrix_;
    }

    //! return reference to matrix
    Matrix& operator* ()
    {
        return *matrix_;
    }

    void assemble(LocalJacobian& loc, SolutionFunction& u, SolutionFunction& f)
    {
    	*matrix_ = 0;
    	gridOperatorSpace_->jacobian(*u, *matrix_);
    	*f = 0;
    	gridOperatorSpace_->residual(*u, *f);
		set_constrained_dofs(*constraintsTrafo_, 0.0, *f);
		set_constrained_dofs(*constraintsTrafo_, 0.0, *u);
    }

    const GridFunctionSpace& gridFunctionSpace() const
    {
        return *gridFunctionSpace_;
    }

    const ConstraintsTrafo& constraintsTrafo() const
    {
        return *constraintsTrafo_;
    }

private:
	Problem& problem_;
	std::vector<int> intghost_;
	CON *cn_;
    FEM *fem_;
    ScalarGridFunctionSpace *scalarGridFunctionSpace_;
    GridFunctionSpace *gridFunctionSpace_;
//    ScalarBTypes *scalarBTypes0_;
//    ScalarBTypes *scalarBTypes1_;
    BTypes *bTypes_;
    ConstraintsTrafo *constraintsTrafo_;
    LocalOperator *localOperator_;
    GridOperatorSpace *gridOperatorSpace_;
    Matrix *matrix_;
};

#endif

#endif
