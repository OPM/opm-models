// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::FemSparseMatrixAdapter
 */
#ifndef OPM_FEM_SPARSE_MATRIX_ADAPTER_HH
#define OPM_FEM_SPARSE_MATRIX_ADAPTER_HH

// this code only works with dune-fem available
#if HAVE_DUNE_FEM

// the following implementation of FemSparseMatrixAdapter only works for
// dune-fem version 2.7 or higher
#if DUNE_VERSION_NEWER(DUNE_FEM, 2, 7)
#include <dune/fem/function/blockvectorfunction.hh>

#if HAVE_PETSC
#include <dune/fem/operator/linear/petscoperator.hh>
#endif

#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/operator/linear/spoperator.hh>


namespace Opm {
namespace Linear {

/*!
 * \ingroup Linear
 * \brief A sparse matrix interface backend for linear operators from dune-fem.
 *
 * \note LinearOperators from dune-fem implement most methods needed for SparseMatrixAdapter
 *       and here we simply add a few forwarding methods.
 */
template <class LinearOperator>
struct FemSparseMatrixAdapter : public LinearOperator
{
    typedef LinearOperator  ParentType;
    typedef typename LinearOperator :: MatrixType    Matrix;
    typedef typename ParentType :: MatrixBlockType   MatrixBlock;

    typedef typename LinearOperator :: RangeFunctionType :: RangeFieldType  Scalar;

    template <class Simulator>
    FemSparseMatrixAdapter( const Simulator& simulator )
        : LinearOperator("Opm::Jacobian", simulator.model().space(), simulator.model().space() )
    {}

    void commit()
    {
      this->flushAssembly();
    }

    template< class LocalBlock >
    void addToBlock ( const size_t row, const size_t col, const LocalBlock& block )
    {
      this->addBlock( row, col, block );
    }

    void clearRow( const size_t row, const Scalar diag = 1.0 ) { this->unitRow( row ); }
};

template <class DiscreteFunction>
using FemSparseRowMatrixAdapter = FemSparseMatrixAdapter< Dune::Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction > >;

#if HAVE_PETSC
template <class DiscreteFunction>
using FemPetscMatrixAdapter = FemSparseMatrixAdapter< Dune::Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction > >;
#endif

#if HAVE_DUNE_ISTL
template <class DiscreteFunction>
using FemISTLMatrixAdapter = FemSparseMatrixAdapter< Dune::Fem::ISTLLinearOperator< DiscreteFunction, DiscreteFunction > >;
#endif

}} // namespace Linear, Opm

#endif // DUNE_VERSION_NEWER(DUNE_FEM, 2, 7)
#endif // HAVE_DUNE_FEM
#endif
