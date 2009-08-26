#ifdef HAVE_DUNE_PDELAB

#ifndef DUMUX_BOUNDARYTYPESPDELAB_HH
#define DUMUX_BOUNDARYTYPESPDELAB_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>

template<typename TypeTag, typename Problem>
class BoundaryTypesPDELab
: public Dune::PDELab::BoundaryGridFunctionBase<
	Dune::PDELab::BoundaryGridFunctionTraits<
		typename GET_PROP_TYPE(TypeTag, PTAG(GridView)),
		int,
		GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
		Dune::FieldVector<int,GET_PROP_VALUE(TypeTag, PTAG(NumEq))> >,
	BoundaryTypesPDELab<TypeTag, Problem> >
{

public:
	enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    enum{dim = GridView::dimension};
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef typename GridView::IntersectionIterator                    IntersectionIterator;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) RefElemProp;
    typedef typename RefElemProp::Container                     ReferenceElements;
    typedef typename RefElemProp::ReferenceElement              ReferenceElement;

    typedef Dune::PDELab::BoundaryGridFunctionTraits<GridView,int,numEq,Dune::FieldVector<int,numEq> > Traits;
	typedef Dune::PDELab::BoundaryGridFunctionBase<Traits,BoundaryTypesPDELab<TypeTag, Problem> > BaseT;

	BoundaryTypesPDELab (Problem& problem)
	: problem_(problem)
	{}

	template<typename I>
	inline void evaluate (const I& ig, const typename Traits::DomainType& x,
			typename Traits::RangeType& y) const
	{
		const Element& element = *(ig.inside());
		problem_.model().localJacobian().setCurrentElement(element);
		const FVElementGeometry& fvGeom =  problem_.model().localJacobian().curFvElementGeometry();

		BoundaryTypeVector values;

        IntersectionIterator isIt = element.ileafbegin();
        const IntersectionIterator &endIt = element.ileafend();
        for (; isIt != endIt; ++isIt)
        {
        	if (isIt->indexInInside() == ig.indexInInside())
        		break;
        }

        y = 1; // Dirichlet
        Dune::GeometryType geoType = element.geometry().type();
        const ReferenceElement& refElem = ReferenceElements::general(geoType);
        int faceIdx = isIt->indexInInside();
        for (int faceVertIdx = 0; faceVertIdx < isIt->geometry().corners(); faceVertIdx++)
        {
        	int scvIdx = refElem.subEntity(faceIdx, 1, faceVertIdx, dim);
        	int boundaryFaceIdx = fvGeom.boundaryFaceIndex(faceIdx, faceVertIdx);

        	problem_.boundaryTypes(values, element, fvGeom, isIt, scvIdx, boundaryFaceIdx);

        	for (unsigned int comp = 0; comp < numEq; comp++)
        		if (values[comp] == Dune::BoundaryConditions::neumann)
        			y[comp] = 0; // Neumann
        }
	}

	inline const GridView& getGridView () {
		return problem_.gridView();
	}
private:
	Problem& problem_;
};

#endif

#endif
