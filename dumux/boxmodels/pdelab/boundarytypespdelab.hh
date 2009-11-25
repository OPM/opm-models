#if HAVE_DUNE_PDELAB && !defined DUMUX_BOUNDARYTYPESPDELAB_HH
#define DUMUX_BOUNDARYTYPESPDELAB_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>

template <class TypeTag>
class BoundaryIndexHelperPDELab
{
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
public:
    template <int i>
    struct Child {
        struct Type {
            enum { isLeaf = true };
            enum { eqIdx = i };
        };
    };
    
    template <int i>
    const typename Child<i>::Type &getChild() const
    {
        static typename Child<i>::Type dummy;
        return dummy;
    };
    
    enum { isLeaf = false };
    enum { CHILDREN = numEq };

    BoundaryIndexHelperPDELab()
	{}
};

#endif
