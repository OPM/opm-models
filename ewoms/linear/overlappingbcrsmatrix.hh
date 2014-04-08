/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::Linear::OverlappingBCRSMatrix
 */
#ifndef EWOMS_OVERLAPPING_BCRS_MATRIX_HH
#define EWOMS_OVERLAPPING_BCRS_MATRIX_HH

#include <ewoms/linear/domesticoverlapfrombcrsmatrix.hh>
#include <ewoms/linear/globalindices.hh>
#include <ewoms/parallel/mpibuffer.hh>

#include <opm/material/Valgrind.hpp>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/io.hh>

#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <vector>
#include <memory>

namespace Ewoms {
namespace Linear {

/*!
 * \brief An overlap aware block-compressed row storage (BCRS) matrix.
 */
template <class BCRSMatrix>
class OverlappingBCRSMatrix : public BCRSMatrix
{
    typedef BCRSMatrix ParentType;

public:
    typedef Ewoms::Linear::DomesticOverlapFromBCRSMatrix<BCRSMatrix> Overlap;

private:
    typedef std::vector<std::set<Index> > Entries;

public:
    typedef typename ParentType::ColIterator ColIterator;
    typedef typename ParentType::ConstColIterator ConstColIterator;
    typedef typename ParentType::block_type block_type;

    // no real copying done at the moment
    OverlappingBCRSMatrix(const OverlappingBCRSMatrix &M) : ParentType(M)
    {}

    OverlappingBCRSMatrix(const BCRSMatrix &M, const BorderList &borderList,
                          const std::set<Index> &blackList, int overlapSize)
    {
        overlap_ = std::shared_ptr<Overlap>(
            new Overlap(M, borderList, blackList, overlapSize));
        myRank_ = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
#endif // HAVE_MPI

        // build the overlapping matrix from the non-overlapping
        // matrix and the overlap
        build_(M);
    }

    ~OverlappingBCRSMatrix()
    {
        if (overlap_.use_count() == 0)
            return;

        // delete all MPI buffers
        const PeerSet &peerSet = overlap_->peerSet();
        typename PeerSet::const_iterator peerIt = peerSet.begin();
        typename PeerSet::const_iterator peerEndIt = peerSet.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            delete rowSizesRecvBuff_[peerRank];
            delete rowIndicesRecvBuff_[peerRank];
            delete entryColIndicesRecvBuff_[peerRank];
            delete entryValuesRecvBuff_[peerRank];

            delete numRowsSendBuff_[peerRank];
            delete rowSizesSendBuff_[peerRank];
            delete rowIndicesSendBuff_[peerRank];
            delete entryColIndicesSendBuff_[peerRank];
            delete entryValuesSendBuff_[peerRank];
        }
    }

    /*!
     * \brief Returns the domestic overlap for the process.
     */
    const Overlap &overlap() const
    { return *overlap_; }

    /*!
     * \brief Assign and syncronize the overlapping matrix from a
     *       non-overlapping one.
     */
    void assignAdd(const BCRSMatrix &M)
    {

        // copy the native entries
        assignFromNative_(M);

        // communicate and add the contents of overlapping rows
        syncAdd_();
    }

    /*!
     * \brief Assign and syncronize the overlapping matrix from a
     *       non-overlapping one.
     *
     * The non-master entries are copied from the master
     */
    void assignCopy(const BCRSMatrix &M)
    {
        // copy the native entries
        assignFromNative_(M);

        // communicate and add the contents of overlapping rows
        syncCopy_();
    }

    /*!
     * \brief Set the identity matrix on the main diagonal of front indices.
     */
    void resetFront()
    {
        // create an identity matrix
        block_type idMatrix(0.0);
        for (unsigned i = 0; i < idMatrix.size(); ++i)
            idMatrix[i][i] = 1.0;

        int numLocal = overlap_->numLocal();
        int numDomestic = overlap_->numDomestic();
        for (int domRowIdx = numLocal; domRowIdx < numDomestic; ++domRowIdx) {
            if (overlap_->isFront(domRowIdx)) {
                // set the front rows to a diagonal 1
                (*this)[domRowIdx] = 0.0;
                (*this)[domRowIdx][domRowIdx] = idMatrix;
            }
        }
    }

    void print() const
    {
        overlap_->print();

        for (int i = 0; i < this->N(); ++i) {
            if (overlap_->isLocal(i))
                std::cout << " ";
            else
                std::cout << "*";
            std::cout << "row " << i << " ";

            typedef typename BCRSMatrix::ConstColIterator ColIt;
            ColIt colIt = (*this)[i].begin();
            ColIt colEndIt = (*this)[i].end();
            for (int j = 0; j < this->M(); ++j) {
                if (colIt != colEndIt && j == colIt.index()) {
                    ++colIt;
                    if (overlap_->isBorder(j))
                        std::cout << "|";
                    else if (overlap_->isLocal(j))
                        std::cout << "X";
                    else
                        std::cout << "*";
                }
                else
                    std::cout << " ";
            }
            std::cout << "\n" << std::flush;
        };
        Dune::printSparseMatrix(std::cout,
                                *static_cast<const BCRSMatrix *>(this), "M",
                                "row");
    }

private:
    void assignFromNative_(const BCRSMatrix &M)
    {
        // first, set everything to 0,
        BCRSMatrix::operator=(0.0);

        // then copy the local entries of M to the overlapping matrix
        for (unsigned nativeRowIdx = 0; nativeRowIdx < M.N(); ++nativeRowIdx) {
            int localRowIdx = overlap_->nativeToDomestic(nativeRowIdx);
            if (localRowIdx < 0) {
                continue; // row corresponds to a black-listed entry
            }

            ConstColIterator nativeColIt = M[nativeRowIdx].begin();
            const ConstColIterator &nativeColEndIt = M[nativeRowIdx].end();
            ColIterator localColIt = (*this)[localRowIdx].begin();
            for (; nativeColIt != nativeColEndIt; ++nativeColIt) {
                int localColIdx = overlap_->nativeToDomestic(nativeColIt.index());
                if (localColIdx < 0) {
                    continue; // column corresponds to a black-listed entry
                }

                while (true) {
                    // go to the local column which corresponds the
                    // the same entry of the native matrix
                    if (static_cast<int>(localColIt.index()) == localColIdx) {
                        (*localColIt) = *nativeColIt;
                        break;
                    }

                    ++localColIt;
                }
            }
        }
    }

    void build_(const BCRSMatrix &M)
    {
        int numDomestic = overlap_->numDomestic();

        // allocate the rows
        this->setSize(numDomestic, numDomestic);
        this->setBuildMode(ParentType::random);

        // communicate the entries
        buildIndices_(M);
    }

    int numDomesticEntriesInRowFor_(const BCRSMatrix &M, int peerRank, int rowIdx)
    {
        int numEntries = 0;

        typedef typename BCRSMatrix::ConstColIterator ColIt;
        ColIt colIt = M[rowIdx].begin();
        ColIt colEndIt = M[rowIdx].end();
        for (; colIt != colEndIt; ++colIt) {
            if (overlap_->isDomesticIndexFor(peerRank, colIt.index()))
                ++numEntries;
        }

        return numEntries;
    }

    void buildIndices_(const BCRSMatrix &M)
    {
        /////////
        // first, add all local matrix entries
        /////////
        entries_.resize(overlap_->numDomestic());
        for (unsigned nativeRowIdx = 0; nativeRowIdx < M.N(); ++nativeRowIdx) {
            int localRowIdx = overlap_->nativeToDomestic(nativeRowIdx);
            if (localRowIdx < 0)
                continue;

            ConstColIterator colIt = M[nativeRowIdx].begin();
            ConstColIterator colEndIt = M[nativeRowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int localColIdx = overlap_->nativeToDomestic(colIt.index());
                if (localColIdx < 0)
                    continue;
                entries_[localRowIdx].insert(localColIdx);
            }
        }

        /////////
        // add the indices for all additional entries
        /////////

        // first, send all our indices to all peers
        const PeerSet &peerSet = overlap_->peerSet();
        typename PeerSet::const_iterator peerIt = peerSet.begin();
        typename PeerSet::const_iterator peerEndIt = peerSet.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendIndices_(M, peerRank);
        }

        // then recieve all indices from the peers
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveIndices_(peerRank);
        }

        // wait until all send operations are completed
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            numRowsSendBuff_[peerRank]->wait();
            rowSizesSendBuff_[peerRank]->wait();
            rowIndicesSendBuff_[peerRank]->wait();
            entryColIndicesSendBuff_[peerRank]->wait();

            // convert the global indices in the send buffers to domestic
            // ones
            globalToDomesticBuff_(*rowIndicesSendBuff_[peerRank]);
            globalToDomesticBuff_(*entryColIndicesSendBuff_[peerRank]);
        }

        /////////
        // actually initialize the BCRS matrix structure
        /////////

        // set the row sizes
        int numDomestic = overlap_->numDomestic();
        for (int rowIdx = 0; rowIdx < numDomestic; ++rowIdx) {
            this->setrowsize(rowIdx, entries_[rowIdx].size());
        }
        this->endrowsizes();

        // set the indices
        for (int rowIdx = 0; rowIdx < numDomestic; ++rowIdx) {
            const auto &colIndices = entries_[rowIdx];

            auto colIdxIt = colIndices.begin();
            const auto &colIdxEndIt = colIndices.end();
            for (; colIdxIt != colIdxEndIt; ++colIdxIt)
                this->addindex(rowIdx, *colIdxIt);
        }
        this->endindices();

        // free the memory occupied by the array of the matrix entries
        entries_.resize(0);
    }

    // send the overlap indices to a peer
    void sendIndices_(const BCRSMatrix &M, int peerRank)
    {
#if HAVE_MPI
        // send size of foreign overlap to peer
        int numOverlapRows = overlap_->foreignOverlapSize(peerRank);
        numRowsSendBuff_[peerRank] = new MpiBuffer<Index>(1);
        (*numRowsSendBuff_[peerRank])[0] = numOverlapRows;
        numRowsSendBuff_[peerRank]->send(peerRank);

        rowSizesSendBuff_[peerRank] = new MpiBuffer<Index>(numOverlapRows);
        rowIndicesSendBuff_[peerRank] = new MpiBuffer<Index>(numOverlapRows);

        for (int rowIdx = 0; rowIdx < numOverlapRows; ++ rowIdx)
            (*rowSizesSendBuff_[peerRank])[rowIdx] = 0;

        // create the row size and index MPI buffers for the peer
        int numEntries = 0;
        for (int overlapOffset = 0; overlapOffset < numOverlapRows; ++overlapOffset) {
            int localRowIdx = overlap_->foreignOverlapOffsetToDomesticIdx(peerRank, overlapOffset);
            int nativeRowIdx = overlap_->domesticToNative(localRowIdx);

            typedef typename BCRSMatrix::ConstColIterator ColIt;
            ColIt colIt = M[nativeRowIdx].begin();
            ColIt colEndIt = M[nativeRowIdx].end();
            int numCols = 0;
            for (; colIt != colEndIt; ++colIt) {
                int localColIdx = overlap_->nativeToDomestic(colIt.index());
                if (localColIdx < 0)
                    continue;
                else if (!overlap_->peerHasIndex(peerRank, localColIdx)) {
                    continue;
                }
                ++numCols;
            }

            (*rowSizesSendBuff_[peerRank])[overlapOffset] = numCols;
            (*rowIndicesSendBuff_[peerRank])[overlapOffset]
                = overlap_->domesticToGlobal(localRowIdx);
            numEntries += numCols;

            Valgrind::CheckDefined((*rowSizesSendBuff_[peerRank])[overlapOffset]);
            Valgrind::CheckDefined(numCols);
        }

        // actually communicate with the peer
        rowSizesSendBuff_[peerRank]->send(peerRank);
        rowIndicesSendBuff_[peerRank]->send(peerRank);

        // create and fill the MPI buffer for the indices of the
        // matrix entries
        entryColIndicesSendBuff_[peerRank] = new MpiBuffer<Index>(numEntries);
        int buffIdx = 0;
        for (int overlapOffset = 0; overlapOffset < numOverlapRows; ++overlapOffset) {
            int localRowIdx = overlap_->foreignOverlapOffsetToDomesticIdx(peerRank, overlapOffset);
            int nativeRowIdx = overlap_->domesticToNative(localRowIdx);

            typedef typename BCRSMatrix::ConstColIterator ColIt;
            ColIt colIt = M[nativeRowIdx].begin();
            ColIt colEndIt = M[nativeRowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int localColIdx = overlap_->nativeToDomestic(colIt.index());
                if (localColIdx < 0)
                    continue;
                else if (!overlap_->peerHasIndex(peerRank, localColIdx))
                    continue;

                int globalColIdx = overlap_->domesticToGlobal(localColIdx);
                (*entryColIndicesSendBuff_[peerRank])[buffIdx] = globalColIdx;
                ++ buffIdx;
            }
        }

        // actually send the column indices for the overlap rows to
        // the peer
        entryColIndicesSendBuff_[peerRank]->send(peerRank);

        // create the send buffers for the values of the matrix
        // entries
        entryValuesSendBuff_[peerRank] = new MpiBuffer<block_type>(numEntries);
#endif // HAVE_MPI
    }

    // receive the overlap indices to a peer
    void receiveIndices_(int peerRank)
    {
#if HAVE_MPI
        // receive size of foreign overlap to peer
        Index numOverlapRows;
        auto &numRowsRecvBuff = numRowsRecvBuff_[peerRank];
        numRowsRecvBuff.resize(1);
        numRowsRecvBuff.receive(peerRank);
        numOverlapRows = numRowsRecvBuff[0];

        // create receive buffer for the row sizes and receive them
        // from the peer
        rowSizesRecvBuff_[peerRank] = new MpiBuffer<Index>(numOverlapRows);
        rowIndicesRecvBuff_[peerRank] = new MpiBuffer<Index>(numOverlapRows);
        rowSizesRecvBuff_[peerRank]->receive(peerRank);
        rowIndicesRecvBuff_[peerRank]->receive(peerRank);

        // calculate the total number of indices which are send by the
        // peer
        int totalIndices = 0;
        for (Index i = 0; i < numOverlapRows; ++i) {
            totalIndices += (*rowSizesRecvBuff_[peerRank])[i];
        }

        // create the buffer to store the column indices of the matrix entries
        entryColIndicesRecvBuff_[peerRank] = new MpiBuffer<Index>(totalIndices);
        entryValuesRecvBuff_[peerRank] = new MpiBuffer<block_type>(totalIndices);

        // communicate with the peer
        entryColIndicesRecvBuff_[peerRank]->receive(peerRank);

        // convert the global indices in the receive buffers to
        // domestic ones
        globalToDomesticBuff_(*rowIndicesRecvBuff_[peerRank]);
        globalToDomesticBuff_(*entryColIndicesRecvBuff_[peerRank]);

        // add the entries to the global entry map
        int k = 0;
        for (Index i = 0; i < numOverlapRows; ++i) {
            int domRowIdx = (*rowIndicesRecvBuff_[peerRank])[i];
            for (Index j = 0; j < (*rowSizesRecvBuff_[peerRank])[i]; ++j) {
                int domColIdx = (*entryColIndicesRecvBuff_[peerRank])[k];
                entries_[domRowIdx].insert(domColIdx);
                ++k;
            }
        }
#endif // HAVE_MPI
    }

    // communicates and adds up the contents of overlapping rows
    void syncAdd_()
    {
        // first, send all entries to the peers
        const PeerSet &peerSet = overlap_->peerSet();
        typename PeerSet::const_iterator peerIt = peerSet.begin();
        typename PeerSet::const_iterator peerEndIt = peerSet.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            sendEntries_(peerRank);
        }

        // then, receive entries from the peers
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            receiveAddEntries_(peerRank);
        }

        // finally, make sure that everything which we send was
        // received by the peers
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            entryValuesSendBuff_[peerRank]->wait();
        }
    }

    // communicates and copies the contents of overlapping rows from
    // the master
    void syncCopy_()
    {
        // first, send all entries to the peers
        const PeerSet &peerSet = overlap_->foreignOverlap().peerSet();
        typename PeerSet::const_iterator peerIt = peerSet.begin();
        typename PeerSet::const_iterator peerEndIt = peerSet.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            sendEntries_(peerRank);
        }

        // then, receive entries from the peers
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            receiveCopyEntries_(peerRank);
        }

        // finally, make sure that everything which we send was
        // received by the peers
        peerIt = peerSet.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            entryValuesSendBuff_[peerRank]->wait();
        }
    }

    void sendEntries_(int peerRank)
    {
#if HAVE_MPI
        auto &mpiSendBuff = *entryValuesSendBuff_[peerRank];

        auto &mpiRowIndicesSendBuff = *rowIndicesSendBuff_[peerRank];
        auto &mpiRowSizesSendBuff = *rowSizesSendBuff_[peerRank];
        auto &mpiColIndicesSendBuff = *entryColIndicesSendBuff_[peerRank];

        // fill the send buffer
        unsigned k = 0;
        for (unsigned i = 0; i < mpiRowIndicesSendBuff.size(); ++i) {
            Index domRowIdx = mpiRowIndicesSendBuff[i];

            typedef typename ParentType::ConstColIterator ColIt;
            ColIt colIt = (*this)[domRowIdx].begin();
            for (Index j = 0;
                 j < static_cast<Index>(mpiRowSizesSendBuff[i]);
                 ++j)
            {
                // move to the next column which is in the overlap
                Index domColIdx = mpiColIndicesSendBuff[k];
                for (; static_cast<Index>(colIt.index()) < domColIdx; ++colIt)
                { }
                assert(static_cast<Index>(colIt.index()) == domColIdx);

                // add the values of this column to the send buffer
                mpiSendBuff[k] = (*colIt);
                ++k;
            }
        }

        mpiSendBuff.send(peerRank);
#endif // HAVE_MPI
    }

    void receiveAddEntries_(int peerRank)
    {
#if HAVE_MPI
        auto &mpiRecvBuff = *entryValuesRecvBuff_[peerRank];

        auto &mpiRowIndicesRecvBuff = *rowIndicesRecvBuff_[peerRank];
        auto &mpiRowSizesRecvBuff = *rowSizesRecvBuff_[peerRank];
        auto &mpiColIndicesRecvBuff = *entryColIndicesRecvBuff_[peerRank];

        mpiRecvBuff.receive(peerRank);

        // retrieve the values from the receive buffer
        int k = 0;
        for (Index i = 0; i < static_cast<Index>(mpiRowIndicesRecvBuff.size()); ++i) {
            Index domRowIdx = mpiRowIndicesRecvBuff[i];
            for (Index j = 0; j < static_cast<Index>(mpiRowSizesRecvBuff[i]); ++j) {
                Index domColIdx = mpiColIndicesRecvBuff[k];

                (*this)[domRowIdx][domColIdx] += mpiRecvBuff[k];
                ++k;
            }
        }
#endif // HAVE_MPI
    }

    void receiveCopyEntries_(int peerRank)
    {
#if HAVE_MPI
        MpiBuffer<block_type> &mpiRecvBuff = *entryValuesRecvBuff_[peerRank];

        MpiBuffer<Index> &mpiRowIndicesRecvBuff = *rowIndicesRecvBuff_[peerRank];
        MpiBuffer<Index> &mpiRowSizesRecvBuff = *rowSizesRecvBuff_[peerRank];
        MpiBuffer<Index> &mpiColIndicesRecvBuff
            = *entryColIndicesRecvBuff_[peerRank];

        mpiRecvBuff.receive(peerRank);

        // retrieve the values from the receive buffer
        int k = 0;
        for (int i = 0; i < mpiRowIndicesRecvBuff.size(); ++i) {
            Index domRowIdx = mpiRowIndicesRecvBuff[i];
            for (int j = 0; j < mpiRowSizesRecvBuff[i]; ++j) {
                Index domColIdx = mpiColIndicesRecvBuff[k];

                if (!overlap_->iAmMasterOf(domRowIdx)
                    || !overlap_->iAmMasterOf(domColIdx)) {
                    (*this)[domRowIdx][domColIdx] = mpiRecvBuff[k];
                }

                ++k;
            }
        }
#endif // HAVE_MPI
    }

    void globalToDomesticBuff_(MpiBuffer<Index> &idxBuff)
    {
        for (unsigned i = 0; i < idxBuff.size(); ++i) {
            idxBuff[i] = overlap_->globalToDomestic(idxBuff[i]);
        }
    }

    int myRank_;
    Entries entries_;
    std::shared_ptr<Overlap> overlap_;

    std::map<ProcessRank, MpiBuffer<Index> *> numRowsSendBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> rowSizesSendBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> rowIndicesSendBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> entryColIndicesSendBuff_;
    std::map<ProcessRank, MpiBuffer<block_type> *> entryValuesSendBuff_;

    std::map<ProcessRank, MpiBuffer<Index> > numRowsRecvBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> rowSizesRecvBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> rowIndicesRecvBuff_;
    std::map<ProcessRank, MpiBuffer<Index> *> entryColIndicesRecvBuff_;
    std::map<ProcessRank, MpiBuffer<block_type> *> entryValuesRecvBuff_;
};

} // namespace Linear
} // namespace Ewoms

#endif
