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
 * \copydoc Ewoms::Linear::OverlappingBlockVector
 */
#ifndef EWOMS_OVERLAPPING_BLOCK_VECTOR_HH
#define EWOMS_OVERLAPPING_BLOCK_VECTOR_HH

#include "overlaptypes.hh"

#include <ewoms/parallel/mpibuffer.hh>
#include <opm/material/Valgrind.hpp>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <memory>
#include <map>
#include <iostream>

namespace Ewoms {
namespace Linear {

/*!
 * \brief An overlap aware block vector.
 */
template <class FieldVector, class Overlap>
class OverlappingBlockVector : public Dune::BlockVector<FieldVector>
{
    typedef Dune::BlockVector<FieldVector> ParentType;
    typedef Dune::BlockVector<FieldVector> BlockVector;

public:
    /*!
     * \brief Given a domestic overlap object, create an overlapping
     *        block vector coherent to it.
     */
    OverlappingBlockVector(const Overlap &overlap)
        : ParentType(overlap.numDomestic()), overlap_(&overlap)
    { createBuffers_(); }

    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector(const OverlappingBlockVector &obv)
        : ParentType(obv), numIndicesSendBuff_(obv.numIndicesSendBuff_),
          indicesSendBuff_(obv.indicesSendBuff_),
          indicesRecvBuff_(obv.indicesRecvBuff_),
          valuesSendBuff_(obv.valuesSendBuff_),
          valuesRecvBuff_(obv.valuesRecvBuff_),
          numFrontIndicesSendBuff_(obv.numFrontIndicesSendBuff_),
          frontIndicesSendBuff_(obv.frontIndicesSendBuff_),
          frontIndicesRecvBuff_(obv.frontIndicesRecvBuff_),
          frontValuesSendBuff_(obv.frontValuesSendBuff_),
          frontValuesRecvBuff_(obv.frontValuesRecvBuff_), overlap_(obv.overlap_)
    {}

    /*!
     * \brief Default constructor.
     */
    OverlappingBlockVector()
    {}

    //! \cond SKIP
    /*!
     * \brief Recycle the assignment operators of Dune::BlockVector
     */
    using ParentType::operator=;
    //! \endcond

    /*!
     * \brief Assignment operator.
     */
    OverlappingBlockVector &operator=(const OverlappingBlockVector &obv)
    {
        ParentType::operator=(obv);
        numIndicesSendBuff_ = obv.numIndicesSendBuff_;
        indicesSendBuff_ = obv.indicesSendBuff_;
        indicesRecvBuff_ = obv.indicesRecvBuff_;
        valuesSendBuff_ = obv.valuesSendBuff_;
        valuesRecvBuff_ = obv.valuesRecvBuff_;
        numFrontIndicesSendBuff_ = obv.numFrontIndicesSendBuff_;
        frontIndicesSendBuff_ = obv.frontIndicesSendBuff_;
        frontIndicesRecvBuff_ = obv.frontIndicesRecvBuff_;
        frontValuesSendBuff_ = obv.frontValuesSendBuff_;
        frontValuesRecvBuff_ = obv.frontValuesRecvBuff_;
        overlap_ = obv.overlap_;
        return *this;
    }

    /*!
     * \brief Assign an overlapping block vector from a
     *        non-overlapping one, border entries are added.
     */
    void assignAddBorder(const BlockVector &nativeBlockVector)
    {
        int numDomestic = overlap_->numDomestic();

        // assign the local rows from the non-overlapping block vector
        for (int domRowIdx = 0; domRowIdx < numDomestic; ++domRowIdx) {
            int nativeRowIdx = overlap_->domesticToNative(domRowIdx);
            if (nativeRowIdx < 0)
                (*this)[domRowIdx] = 0.0;
            else
                (*this)[domRowIdx] = nativeBlockVector[nativeRowIdx];
        }

        // add up the contents of border rows, for the remaining rows,
        // get the values from their respective master process.
        syncAddBorder();
    }

    /*!
     * \brief Assign the local values to a non-overlapping block
     *        vector.
     */
    void assignTo(BlockVector &nativeBlockVector) const
    {
        // assign the local rows
        int numNative = overlap_->numNative();
        nativeBlockVector.resize(numNative);
        for (int nativeRowIdx = 0; nativeRowIdx < numNative; ++nativeRowIdx) {
            int domRowIdx = overlap_->nativeToDomestic(nativeRowIdx);

            if (domRowIdx < 0)
                nativeBlockVector[nativeRowIdx] = 0.0;
            else
                nativeBlockVector[nativeRowIdx] = (*this)[domRowIdx];
        }
    }

    /*!
     * \brief Syncronize all values of the block vector from their
     *        master process.
     */
    void sync()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendEntries_(peerRank);
        }

        // recieve all entries to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveFromMaster_(peerRank);
        }

        // wait until we have send everything
        waitSendFinished_();
    }

    /*!
     * \brief Syncronize all values of the block vector by adding up
     *        the values of all peer ranks.
     */
    void syncAdd()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendEntries_(peerRank);
        }

        // recieve all entries to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveAdd_(peerRank);
        }

        // wait until we have send everything
        waitSendFinished_();
    }

    /*!
     * \brief Syncronize the values of the front indices by copying
     *        them from their respective master.
     */
    void syncFront()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendFrontEntries_(peerRank);
        }

        // recieve all entries to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveFront_(peerRank);
        }

        // wait until we have send everything
        waitSendFinished_();
    }

    /*!
     * \brief Syncronize all values of the block vector from the
     *        master rank, but add up the entries on the border.
     */
    void syncAddBorder()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendEntries_(peerRank);
        }

        // recieve all entries to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveAddBorder_(peerRank);
        }

        // wait until we have send everything
        waitSendFinished_();

    }

    /*!
     * \brief Set all front indices to 0
     */
    void resetFront()
    {
        // set all front rows to 0
        int nLocal = overlap_->numLocal();
        int nDomestic = overlap_->numDomestic();
        for (int j = nLocal; j < nDomestic; ++j) {
            if (overlap_->isFront(j))
                (*this)[j] = 0.0;
        }
    }

    void print() const
    {
        for (int i = 0; i < this->size(); ++i) {
            std::cout << "row " << i << (overlap_->isLocal(i) ? " " : "*")
                      << ": " << (*this)[i] << "\n" << std::flush;
        }
    }

private:
    void createBuffers_()
    {
#if HAVE_MPI
        // create array for the front indices
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all indices to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            int numEntries = overlap_->foreignOverlapSize(peerRank);
            numIndicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<int> >(1);
            indicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<Index> >(numEntries);
            valuesSendBuff_[peerRank] = std::make_shared<MpiBuffer<FieldVector> >(numEntries);

            // fill the indices buffer with global indices
            MpiBuffer<Index> &indicesSendBuff = *indicesSendBuff_[peerRank];
            for (int i = 0; i < numEntries; ++i) {
                int domRowIdx = overlap_->foreignOverlapOffsetToDomesticIdx(peerRank, i);
                indicesSendBuff[i] = overlap_->domesticToGlobal(domRowIdx);
            }

            // first, send the number of indices
            (*numIndicesSendBuff_[peerRank])[0] = numEntries;
            numIndicesSendBuff_[peerRank]->send(peerRank);

            // then, send the indices themselfs
            indicesSendBuff.send(peerRank);
        }

        // receive the indices from the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            // receive size of overlap to peer
            MpiBuffer<int> numRowsRecvBuff(1);
            numRowsRecvBuff.receive(peerRank);
            int numRows = numRowsRecvBuff[0];

            // then, create the MPI buffers
            indicesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(
                new MpiBuffer<Index>(numRows));
            valuesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(
                new MpiBuffer<FieldVector>(numRows));
            MpiBuffer<Index> &indicesRecvBuff = *indicesRecvBuff_[peerRank];

            // next, receive the actual indices
            indicesRecvBuff.receive(peerRank);

            // finally, translate the global indices to domestic ones
            for (int i = 0; i != numRows; ++i) {
                int globalRowIdx = indicesRecvBuff[i];
                int domRowIdx = overlap_->globalToDomestic(globalRowIdx);

                indicesRecvBuff[i] = domRowIdx;
            }
        }

        // wait for all send operations to complete
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            numIndicesSendBuff_[peerRank]->wait();
            indicesSendBuff_[peerRank]->wait();

            // convert the global indices of the send buffer to
            // domestic ones
            MpiBuffer<Index> &indicesSendBuff = *indicesSendBuff_[peerRank];
            for (unsigned i = 0; i < indicesSendBuff.size(); ++i) {
                indicesSendBuff[i] = overlap_->globalToDomestic(indicesSendBuff[i]);
            }
        }

        // send all front indices to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            int numFrontEntries = overlap_->numFront(peerRank);
            numFrontIndicesSendBuff_[peerRank] =
                std::shared_ptr<MpiBuffer<int> >(new MpiBuffer<int>(1));
            frontIndicesSendBuff_[peerRank] =
                std::shared_ptr<MpiBuffer<Index> >(new MpiBuffer<Index>(numFrontEntries));
            frontValuesSendBuff_[peerRank] =
                std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numFrontEntries));

            // fill the indices buffer with global indices
            MpiBuffer<Index> &frontIndicesSendBuff = *frontIndicesSendBuff_[peerRank];
            int numOverlapRows = overlap_->foreignOverlapSize(peerRank);
            int mpiBuffOffset = 0;
            for (int overlapOffset = 0; overlapOffset < numOverlapRows; ++overlapOffset) {
                int rowIdx = overlap_->foreignOverlapOffsetToDomesticIdx(peerRank, overlapOffset);
                if (!overlap_->isFrontFor(peerRank, rowIdx))
                    continue;

                frontIndicesSendBuff[mpiBuffOffset] = overlap_->domesticToGlobal(rowIdx);
                ++ mpiBuffOffset;
            }

            // first, send the number of indices
            (*numFrontIndicesSendBuff_[peerRank])[0] = numFrontEntries;
            numFrontIndicesSendBuff_[peerRank]->send(peerRank);

            // then, send the indices themselfs
            frontIndicesSendBuff.send(peerRank);
        }

        // receive the indices from the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            // receive size of overlap to peer
            MpiBuffer<int> numFrontRowsRecvBuff(1);
            numFrontRowsRecvBuff.receive(peerRank);
            int numFrontRows = numFrontRowsRecvBuff[0];

            // then, create the MPI buffers
            frontIndicesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(
                new MpiBuffer<Index>(numFrontRows));
            frontValuesRecvBuff_[peerRank] = std::make_shared<MpiBuffer<FieldVector> >(numFrontRows);
            MpiBuffer<Index> &frontIndicesRecvBuff = *frontIndicesRecvBuff_[peerRank];

            // next, receive the actual indices
            frontIndicesRecvBuff.receive(peerRank);

            // finally, translate the global indices to domestic ones
            for (int i = 0; i < numFrontRows; ++i) {
                int globalRowIdx = frontIndicesRecvBuff[i];
                int domRowIdx = overlap_->globalToDomestic(globalRowIdx);

                frontIndicesRecvBuff[i] = domRowIdx;
            }
        }

        // wait for all send operations to complete
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            numFrontIndicesSendBuff_[peerRank]->wait();
            frontIndicesSendBuff_[peerRank]->wait();

            // convert the global indices of the send buffer to
            // domestic ones
            MpiBuffer<Index> &frontIndicesSendBuff = *frontIndicesSendBuff_[peerRank];
            for (unsigned i = 0; i < frontIndicesSendBuff.size(); ++i) {
                frontIndicesSendBuff[i] = overlap_->globalToDomestic(frontIndicesSendBuff[i]);
            }
        }
#endif // HAVE_MPI
    }

    void sendEntries_(int peerRank)
    {
        // copy the values into the send buffer
        const MpiBuffer<Index> &indices = *indicesSendBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesSendBuff_[peerRank];
        for (unsigned i = 0; i < indices.size(); ++i)
            values[i] = (*this)[indices[i]];

        values.send(peerRank);
    }

    void waitSendFinished_()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            valuesSendBuff_[peerRank]->wait();
        }
    }

    void sendFrontEntries_(int peerRank)
    {
        // copy the values into the send buffer
        const MpiBuffer<Index> &indices = *frontIndicesSendBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *frontValuesSendBuff_[peerRank];
        for (Index i = 0; i < static_cast<Index>(indices.size()); ++i)
            values[i] = (*this)[indices[i]];

        values.send(peerRank);
    }

    void waitSendFrontFinished_()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            frontValuesSendBuff_[peerRank]->wait();
        }
    }

    void receiveFromMaster_(int peerRank)
    {
        const MpiBuffer<Index> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // copy them into the block vector
        for (unsigned j = 0; j < indices.size(); ++j) {
            Index domRowIdx = indices[j];
            if (overlap_->masterRank(domRowIdx) == peerRank) {
                (*this)[domRowIdx] = values[j];
            }
        }
    }

    void receiveFront_(int peerRank)
    {
        const MpiBuffer<Index> &frontIndices = *frontIndicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &frontValues = *frontValuesRecvBuff_[peerRank];

        // receive the values from the peer
        frontValues.receive(peerRank);

        // copy them into the block vector
        for (unsigned j = 0; j < frontIndices.size(); ++j) {
            Index domRowIdx = frontIndices[j];
            (*this)[domRowIdx] = frontValues[j];
        }
    }

    void receiveAddBorder_(int peerRank)
    {
        const MpiBuffer<Index> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // add up the values of rows on the shared boundary
        for (unsigned j = 0; j < indices.size(); ++j) {
            int domRowIdx = indices[j];
            if (overlap_->isBorderWith(domRowIdx, peerRank))
                (*this)[domRowIdx] += values[j];
            else
                (*this)[domRowIdx] = values[j];
        }
    }

    void receiveAdd_(int peerRank)
    {
        const MpiBuffer<Index> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // add up the values of rows on the shared boundary
        for (int j = 0; j < indices.size(); ++j) {
            int domRowIdx = indices[j];
            (*this)[domRowIdx] += values[j];
        }
    }

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<int> > > numIndicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<Index> > > indicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<Index> > > indicesRecvBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesRecvBuff_;

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<int> > > numFrontIndicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<Index> > > frontIndicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<Index> > > frontIndicesRecvBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > >
    frontValuesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > >
    frontValuesRecvBuff_;

    const Overlap *overlap_;
};

} // namespace Linear
} // namespace Ewoms

#endif
