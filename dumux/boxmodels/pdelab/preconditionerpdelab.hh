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



} // namespace PDELab
} // namespace Dune

#endif
