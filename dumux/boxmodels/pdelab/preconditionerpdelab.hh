#ifndef DUMUX_PRECONDITIONERPDELAB_HH
#define DUMUX_PRECONDITIONERPDELAB_HH

#include<dune/pdelab/backend/istlsolverbackend.hh>

namespace Dune {
namespace PDELab {

// wrapped sequential preconditioner
template<class CC, class GFS, class P>
class NonoverlappingWrappedPreconditioner
: public Dune::Preconditioner<typename P::domain_type,typename P::range_type>
{
public:
    //! \brief The domain type of the preconditioner.
    typedef typename P::domain_type domain_type;
    //! \brief The range type of the preconditioner.
    typedef typename P::range_type range_type;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::nonoverlapping
    };

    /*! \brief Constructor.

    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    \param n The number of iterations to perform.
    \param w The relaxation factor.
    */
    NonoverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_,
            const ParallelISTLHelper<GFS>& helper_)
    : gfs(gfs_), prec(prec_), cc(cc_), helper(helper_)
    {}

    /*!
    \brief Prepare the preconditioner.

    \copydoc Preconditioner::pre(domain_type&,range_type&)
    */
    virtual void pre (domain_type& x, range_type& b)
    {
        prec.pre(x,b);
    }

    /*!
    \brief Apply the precondioner.

    \copydoc Preconditioner::apply(domain_type&,const range_type&)
    */
    virtual void apply (domain_type& v, const range_type& d)
    {
        range_type dd(d);
        set_constrained_dofs(cc,0.0,dd);
        prec.apply(v,dd);
        Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
        if (gfs.gridview().comm().size()>1)
            gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
    }

    /*!
    \brief Clean up.

    \copydoc Preconditioner::post(domain_type&)
    */
    virtual void post (domain_type& x)
    {
        prec.post(x);
    }

private:
    const GFS& gfs;
    P& prec;
    const CC& cc;
    const ParallelISTLHelper<GFS>& helper;
};

template<class TypeTag>
class ISTLBackend_NoOverlap_BCGS_ILU
{
    typedef typename GET_PROP(TypeTag, PTAG(PDELabTypes))  PDELabTypes;
    typedef typename PDELabTypes::GridFunctionSpace GridFunctionSpace;
    typedef typename PDELabTypes::ConstraintsTrafo ConstraintsTrafo;
	typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> PHELPER;

public:
	/*! \brief make a linear solver object

	\param[in] gfs a grid function space
	\param[in] maxiter maximum number of iterations to do
	\param[in] verbose print messages if true
	*/
	explicit ISTLBackend_NoOverlap_BCGS_ILU (const GridFunctionSpace& gfs_, const ConstraintsTrafo& constraintsTrafo,
			unsigned maxiter_=5000, int verbose_=1)
	: gfs(gfs_), phelper(gfs), constraintsTrafo_(constraintsTrafo), maxiter(maxiter_), verbose(verbose_)
	{}

	/*! \brief compute global norm of a vector

	\param[in] v the given vector
	*/
	template<class V>
	typename V::ElementType norm (const V& v) const
	{
		V x(v); // make a copy because it has to be made consistent
		typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,V> PSP;
		PSP psp(gfs,phelper);
		psp.make_consistent(x);
		return psp.norm(x);
	}

	/*! \brief solve the given linear system

	\param[in] A the given matrix
	\param[out] z the solution vector to be computed
	\param[in] r right hand side
	\param[in] reduction to be achieved
	*/
	template<class M, class V, class W>
	void apply(M& A, V& z, W& r, typename V::ElementType reduction)
	{
		typedef Dune::PDELab::NonoverlappingOperator<GridFunctionSpace,M,V,W> POP;
		POP pop(gfs,A,phelper);
		typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,V> PSP;
		PSP psp(gfs,phelper);
    	typedef Dune::SeqILU0<M,V,V> SeqPreCond;
    	SeqPreCond seqPreCond(A, 0.9);
    	typedef Dune::PDELab::NonoverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
    	ParPreCond parPreCond(gfs, seqPreCond, constraintsTrafo_, phelper);

//		typedef Dune::PDELab::NonoverlappingRichardson<GridFunctionSpace,V,W> PRICH;
//		PRICH prich(gfs,phelper);
		int verb=0;
		if (gfs.gridview().comm().rank()==0) verb=verbose;
		Dune::BiCGSTABSolver<V> solver(pop,psp,parPreCond,reduction,maxiter,verb);
		Dune::InverseOperatorResult stat;
		solver.apply(z,r,stat);
		res.converged  = stat.converged;
		res.iterations = stat.iterations;
		res.elapsed    = stat.elapsed;
		res.reduction  = stat.reduction;
	}

	/*! \brief Return access to result data */
	const Dune::PDELab::LinearSolverResult<double>& result() const
    		  {
		return res;
    		  }

private:
	const GridFunctionSpace& gfs;
	PHELPER phelper;
	const ConstraintsTrafo& constraintsTrafo_;
	Dune::PDELab::LinearSolverResult<double> res;
	unsigned maxiter;
	int verbose;
};

template<class TypeTag>
class ISTLBackend_NoOverlap_Loop_Pardiso
{
    typedef typename GET_PROP(TypeTag, PTAG(PDELabTypes))  PDELabTypes;
    typedef typename PDELabTypes::GridFunctionSpace GridFunctionSpace;
    typedef typename PDELabTypes::ConstraintsTrafo ConstraintsTrafo;
	typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> PHELPER;

public:
	/*! \brief make a linear solver object

	\param[in] gfs a grid function space
	\param[in] maxiter maximum number of iterations to do
	\param[in] verbose print messages if true
	*/
	explicit ISTLBackend_NoOverlap_Loop_Pardiso (const GridFunctionSpace& gfs_, const ConstraintsTrafo& constraintsTrafo,
			unsigned maxiter_=5000, int verbose_=1)
	: gfs(gfs_), phelper(gfs), constraintsTrafo_(constraintsTrafo), maxiter(maxiter_), verbose(verbose_)
	{}

	/*! \brief compute global norm of a vector

	\param[in] v the given vector
	*/
	template<class V>
	typename V::ElementType norm (const V& v) const
	{
		V x(v); // make a copy because it has to be made consistent
		typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,V> PSP;
		PSP psp(gfs,phelper);
		psp.make_consistent(x);
		return psp.norm(x);
	}

	/*! \brief solve the given linear system

	\param[in] A the given matrix
	\param[out] z the solution vector to be computed
	\param[in] r right hand side
	\param[in] reduction to be achieved
	*/
	template<class M, class V, class W>
	void apply(M& A, V& z, W& r, typename V::ElementType reduction)
	{
		typedef Dune::PDELab::NonoverlappingOperator<GridFunctionSpace,M,V,W> POP;
		POP pop(gfs,A,phelper);
		typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,V> PSP;
		PSP psp(gfs,phelper);
    	typedef Dune::SeqPardiso<M,V,V> SeqPreCond;
    	SeqPreCond seqPreCond(A);
    	typedef Dune::PDELab::NonoverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
    	ParPreCond parPreCond(gfs, seqPreCond, constraintsTrafo_, phelper);

//		typedef Dune::PDELab::NonoverlappingRichardson<GridFunctionSpace,V,W> PRICH;
//		PRICH prich(gfs,phelper);
		int verb=0;
		if (gfs.gridview().comm().rank()==0) verb=verbose;
		Dune::BiCGSTABSolver<V> solver(pop,psp,parPreCond,reduction,maxiter,verb);
		Dune::InverseOperatorResult stat;
		solver.apply(z,r,stat);
		res.converged  = stat.converged;
		res.iterations = stat.iterations;
		res.elapsed    = stat.elapsed;
		res.reduction  = stat.reduction;
	}

	/*! \brief Return access to result data */
	const Dune::PDELab::LinearSolverResult<double>& result() const
    		  {
		return res;
    		  }

private:
	const GridFunctionSpace& gfs;
	PHELPER phelper;
	const ConstraintsTrafo& constraintsTrafo_;
	Dune::PDELab::LinearSolverResult<double> res;
	unsigned maxiter;
	int verbose;
};



} // namespace PDELab
} // namespace Dune

#endif
