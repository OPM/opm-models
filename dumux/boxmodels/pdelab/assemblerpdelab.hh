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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))    Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    enum{dim = GridView::dimension};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalFEMSpace)) FEM;
    typedef typename GET_PROP(TypeTag, PTAG(PDELabTypes)) PDELabTypes;
    typedef typename PDELabTypes::Constraints Constraints;
    typedef typename PDELabTypes::ScalarGridFunctionSpace ScalarGridFunctionSpace;
    typedef typename PDELabTypes::GridFunctionSpace GridFunctionSpace;
	typedef typename PDELabTypes::ConstraintsTrafo ConstraintsTrafo;
    typedef typename PDELabTypes::LocalOperator LocalOperator;
    typedef typename PDELabTypes::GridOperatorSpace GridOperatorSpace;
    typedef BoundaryIndexHelperPDELab<TypeTag> BoundaryFunction;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
	typedef typename GridFunctionSpace::template VectorContainer<Scalar>::Type Vector;
    typedef typename GridOperatorSpace::template MatrixContainer<Scalar>::Type Matrix;
    typedef Matrix RepresentationType;

    AssemblerPDELab(Model &model, Problem& problem)
    : problem_(problem)
    {
    	fem_ = 0;
    	cn_ = 0;
    	scalarGridFunctionSpace_ = 0;
    	gridFunctionSpace_ = 0;
    	bTypes_ = 0;
    	constraintsTrafo_ = 0;
    	localOperator_ = 0;
    	gridOperatorSpace_ = 0;
    	matrix_ = 0;

    	fem_ = new FEM();
    	//cn_ = new Constraints(intghost_);
        cn_ = new Constraints(problem_);
    	scalarGridFunctionSpace_ = new ScalarGridFunctionSpace(problem_.gridView(), *fem_, *cn_);
    	gridFunctionSpace_ = new GridFunctionSpace(*scalarGridFunctionSpace_);

    	Vector ghost(*gridFunctionSpace_, 0.0);
    	Dune::PDELab::GhostDataHandle<GridFunctionSpace, Vector> ghostDataHandle(*gridFunctionSpace_, ghost);
    	if (problem_.gridView().comm().size() > 1)
    		problem_.gridView().communicate(ghostDataHandle,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
    	ghost.std_copy_to(intghost_);

    	bTypes_ = new BoundaryFunction();
    	constraintsTrafo_ = new ConstraintsTrafo();
    	Dune::PDELab::constraints(*bTypes_, *gridFunctionSpace_, *constraintsTrafo_, true);

    	localOperator_ = new LocalOperator(model);
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

    	typedef typename Matrix::RowIterator RowIterator;
    	typedef typename Matrix::ColIterator ColIterator;
    	const typename Matrix::block_type::size_type rowsInBlock = Matrix::block_type::rows;
    	const typename Matrix::block_type::size_type colsInBlock = Matrix::block_type::cols;
    	RowIterator endIBlock = matrix_->end();
    	for (RowIterator iBlock = matrix_->begin(); iBlock != endIBlock; ++iBlock)
    		for (typename Matrix::block_type::size_type iLocal = 0; iLocal < rowsInBlock; iLocal++)
    		{
    			Scalar diagonalEntry = 0;

    			ColIterator endJBlock = iBlock->end();
    			for (ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
    				if (iBlock.index() == jBlock.index())
    				{
    					diagonalEntry = (*jBlock)[iLocal][iLocal];
    					break;
    				}

    			if (std::abs(diagonalEntry) > 1e-12)
    			{
    				for (ColIterator jBlock = iBlock->begin(); jBlock != endJBlock; ++jBlock)
    					for (typename Matrix::block_type::size_type jLocal = 0; jLocal < colsInBlock; jLocal++)
    						(*jBlock)[iLocal][jLocal] /= diagonalEntry;

    				(*f)[iBlock.index()][iLocal] /= diagonalEntry;
    			}
//    			else
//    				std::cout << "assemblerpdelab.hh:assemble(): iBlock = " << iBlock.index() << ", iLocal = " << iLocal << ": ZERO diagonal!" << std::endl;
    		}
    }

    const GridFunctionSpace& gridFunctionSpace() const
    {
        return *gridFunctionSpace_;
    }

    const ConstraintsTrafo& constraintsTrafo() const
    {
        return *constraintsTrafo_;
    }

    ~AssemblerPDELab()
    {
    	delete matrix_;
    	delete gridOperatorSpace_;
    	delete localOperator_;
    	delete constraintsTrafo_;
    	delete bTypes_;
    	delete gridFunctionSpace_;
    	delete scalarGridFunctionSpace_;
    	delete fem_;
    	delete cn_;
    }

private:
	Problem& problem_;
	std::vector<int> intghost_;
	Constraints *cn_;
    FEM *fem_;
    ScalarGridFunctionSpace *scalarGridFunctionSpace_;
    GridFunctionSpace *gridFunctionSpace_;
    BoundaryFunction *bTypes_;
    ConstraintsTrafo *constraintsTrafo_;
    LocalOperator *localOperator_;
    GridOperatorSpace *gridOperatorSpace_;
    Matrix *matrix_;
};

#endif

#endif
