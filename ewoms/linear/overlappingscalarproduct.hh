/*
  Copyright (C) 2011-2013 by Andreas Lauser
  Copyright (C) 2012 by Bernd Flemisch

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::OverlappingScalarProduct
 */
#ifndef EWOMS_OVERLAPPING_SCALAR_PRODUCT_HH
#define EWOMS_OVERLAPPING_SCALAR_PRODUCT_HH

#if HAVE_MPI
#include <mpi.h>
#endif

#include <dune/istl/scalarproducts.hh>

namespace Ewoms {
namespace Linear {

/*!
 * \brief An overlap aware ISTL scalar product.
 */
template <class OverlappingBlockVector, class Overlap>
class OverlappingScalarProduct
    : public Dune::ScalarProduct<OverlappingBlockVector>
{
public:
    typedef typename OverlappingBlockVector::field_type field_type;

    enum { category = Dune::SolverCategory::overlapping };

    OverlappingScalarProduct(const Overlap &overlap) : overlap_(overlap)
    {}

    field_type dot(const OverlappingBlockVector &x,
                   const OverlappingBlockVector &y)
    {
        double sum = 0;
        int numDom = overlap_.numDomestic();
        for (int domIdx = 0; domIdx < numDom; ++domIdx) {
            if (overlap_.iAmMasterOf(domIdx))
                sum += x[domIdx] * y[domIdx];
        }

        // compute the global sum
        double sumGlobal = 0.0;
#if HAVE_MPI
        MPI_Allreduce(&sum,            // source buffer
                      &sumGlobal,      // destination buffer
                      1,               // number of objects in buffers
                      MPI_DOUBLE,      // data type
                      MPI_SUM,         // operation
                      MPI_COMM_WORLD); // communicator
#else
        sumGlobal = sum;
#endif // HAVE_MPI

        return sumGlobal;
    }

    double norm(const OverlappingBlockVector &x)
    { return std::sqrt(dot(x, x)); }

private:
    const Overlap &overlap_;
};

} // namespace Linear
} // namespace Ewoms

#endif
