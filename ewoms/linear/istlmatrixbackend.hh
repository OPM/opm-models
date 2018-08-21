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
 * \copydoc Ewoms::Linear::ISTLMatrixBackend
 */
#ifndef EWOMS_ISTL_MATRIX_BACKEND_HH
#define EWOMS_ISTL_MATRIX_BACKEND_HH

#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Ewoms {
namespace Linear {

/*!
 * \ingroup Linear
 * \brief A backend for BCRSMatrix from dune-istl.
 */
template<class Block, class A=std::allocator< Block > >
class ISTLMatrixBackend
{
    typedef A AllocatorType;
public:
    //! \brief Implementation of matrix
    typedef Dune::BCRSMatrix< Block, AllocatorType >   Matrix;

    //! \brief block type forming the matrix entries (same as Block)
    typedef typename Matrix :: block_type          block_type;

    //! \brief type of scalar
    typedef typename block_type :: field_type      Scalar;

    /*!
     * \brief Constructor creating an empty matrix.
     */
    ISTLMatrixBackend( const size_t rows, const size_t columns )
        : rows_( rows )
        , columns_( columns )
        , matrix_()
    {}

    /*!
     * \brief Constructor taking simulator and creating an empty matrix .
     */
    template < class Simulator >
    ISTLMatrixBackend( const Simulator& simulator )
        : ISTLMatrixBackend( simulator.model().numTotalDof(), simulator.model().numTotalDof() )
    {}

    /*!
     * \brief Allocate matrix structure give a sparsity pattern.
     */
    template <class Set>
    void reserve( const std::vector< Set >& sparsityPattern )
    {
        // allocate raw matrix
        matrix_.reset( new Matrix(rows_, columns_, Matrix::random) );

        // make sure sparsityPattern is consistent with number of rows
        assert( rows_ == sparsityPattern.size() );

        // allocate space for the rows of the matrix
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx)
        {
            matrix_->setrowsize(dofIdx, sparsityPattern[dofIdx].size());
        }

        matrix_->endrowsizes();

        // fill the rows with indices. each degree of freedom talks to
        // all of its neighbors. (it also talks to itself since
        // degrees of freedom are sometimes quite egocentric.)
        for (size_t dofIdx = 0; dofIdx < rows_; ++ dofIdx)
        {
            auto nIt    = sparsityPattern[dofIdx].begin();
            auto nEndIt = sparsityPattern[dofIdx].end();
            for (; nIt != nEndIt; ++nIt)
            {
                matrix_->addindex(dofIdx, *nIt);
            }
        }
        matrix_->endindices();
    }

    /*!
     * \brief Return constant reference to matrix implementation.
     */
    Matrix& matrix() { return *matrix_; }
    const Matrix& matrix() const { return *matrix_; }

    /*!
     * \brief Return number of rows of the matrix.
     */
    size_t rows () const { return rows_; }

    /*!
     * \brief Return number of columns of the matrix.
     */
    size_t cols () const { return columns_; }

    /*!
     * \brief Set all matrix entries to zero.
     */
    void clear()
    {
        (*matrix_) = typename block_type :: field_type(0);
    }

    /*!
     * \brief Set given row to zero except for the diagonal entry which is set to one.
     */
    void clearRow( const size_t row, const Scalar diag = 1.0 )
    {
        block_type idBlock( Scalar(0) );
        for (int i = 0; i < idBlock.rows; ++i)
            idBlock[i][i] = diag;

        auto& matRow = (*matrix_)[ row ];
        auto colIt = matRow.begin();
        const auto& colEndIt = matRow.end();
        for (; colIt != colEndIt; ++colIt)
        {
            if( colIt.index() == row )
                *colIt = idBlock;
            else
                *colIt = Scalar(0);
        }
    }

    /*!
     * \brief Fill given block with entries stored in the matrix.
     */
    void block( const size_t row, const size_t col, block_type& entry ) const
    {
        entry = (*matrix_)[ row ][ col ];
    }

    /*!
     * \brief Set matrix block to given block.
     */
    void setBlock( const size_t row, const size_t col, const block_type& entry )
    {
        (*matrix_)[ row ][ col ] = entry;
    }

    /*!
     * \brief Add block to matrix block.
     */
    void addBlock( const size_t row, const size_t col, const block_type& entry )
    {
        (*matrix_)[ row ][ col ] += entry;
    }

    /*!
     * \brief Flush matrix from local caches into matrix structure Synchronize matrix (empty here)
     */
    void flush( )
    {

    }

    /*!
     * \brief Synchronize matrix and finalize building stage.
     */
    void finalize( )
    {
        // nothing to do here
        // may call compress when implicit build mode is used
    }

protected:
    size_t rows_;
    size_t columns_;

    std::unique_ptr< Matrix > matrix_;
};

} // namespace Linear
} // namespace Ewoms

#endif
