// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUMUX_BOXDIRICHLETCONSTRAINTS_HH
#define DUMUX_BOXDIRICHLETCONSTRAINTS_HH

#include <cstddef>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/grid.hh>
#include<dune/common/geometrytype.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>

#include<dune/disc/operators/boundaryconditions.hh>

namespace Dune {
//! Constraints construction
// works in any dimension and on all element types
template <class TypeTag>
class BoxDirichletConstraints // : public Dune::PDELab::ConformingDirichletConstraints
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

	enum {numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    enum {dim = GridView::dimension};

    typedef Dune::FieldVector<Dune::BoundaryConditions::Flags, numEq> BoundaryTypeVector;
    Problem &problem_;

public:
    enum { doBoundary = true };
    enum { doProcessor = false };
    enum { doSkeleton = false };
    enum { doVolume = false };

    BoxDirichletConstraints(Problem& problem)
        : problem_(problem)
	{}

    //! boundary constraints
    /**
     * \tparam F   grid function returning boundary condition type (this is ignored. we ask the problem directly...)
     * \tparam IG  intersection geometry
     * \tparam LFS local function space
     * \tparam T   TransformationType
     */
    template<typename F, typename I, typename LFS, typename T>
    void boundary (const F& f, const Dune::PDELab::IntersectionGeometry<I>& ig, 
                   const LFS& lfs, T& trafo) const
    {
        typename F::Traits::RangeType bctype;
        
		const Element& element = *(ig.inside());
        FVElementGeometry fvElemGeom;
        fvElemGeom.update(element);
		BoundaryTypeVector values;
       
        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename Dune::PDELab::IntersectionGeometry<I>::ctype DT;
        const Dune::GenericReferenceElement<DT,dim>& refElem = Dune::GenericReferenceElements<DT,dim>::general(gt);

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        int faceIdx = ig.indexInInside();
        for (int faceVertIdx = 0; faceVertIdx < ig.geometry().corners(); faceVertIdx++)
        {
        	int scvIdx = refElem.subEntity(faceIdx, 1, faceVertIdx, dim);
        	int boundaryFaceIdx = fvElemGeom.boundaryFaceIndex(faceIdx, faceVertIdx);

        	problem_.boundaryTypes(values, element, fvElemGeom, ig.intersection(), scvIdx, boundaryFaceIdx);
        	for (unsigned int comp = 0; comp < numEq; comp++)
        		if (values[comp] == Dune::BoundaryConditions::dirichlet) {
                    trafo[scvIdx] = empty; // TODO mixed boundaries
                    break;
                }
        }
    }
};
}

#endif
