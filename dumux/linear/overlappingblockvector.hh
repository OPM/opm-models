// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Dumux::Linear::OverlappingBlockVector
 */
#ifndef DUMUX_OVERLAPPING_BLOCK_VECTOR_HH
#define DUMUX_OVERLAPPING_BLOCK_VECTOR_HH

#include "overlaptypes.hh"

#include <dumux/parallel/mpibuffer.hh>
#include <dumux/common/valgrind.hh>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <memory>
#include <vector>
#include <map>
#include <iostream>

namespace Dumux {
namespace Linear {

/*!
 * \brief An overlap aware block vector.
 */
template <class FieldVector, class Overlap>
class OverlappingBlockVector
    : public Dune::BlockVector<FieldVector>
{
    typedef Dune::BlockVector<FieldVector> ParentType;
    typedef Dune::BlockVector<FieldVector> BlockVector;


public:
    /*!
     * \brief Given a domestic overlap object, create an overlapping
     *        block vector coherent to it.
     */
    OverlappingBlockVector(const Overlap &overlap)
        : ParentType(overlap.numDomestic())
        , overlap_(&overlap)
    {
        createBuffers_();
    }

    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector(const OverlappingBlockVector &obv)
        : ParentType(obv)
        , numIndicesSendBuff_(obv.numIndicesSendBuff_)
        , indicesSendBuff_(obv.indicesSendBuff_)
        , indicesRecvBuff_(obv.indicesRecvBuff_)
        , valuesSendBuff_(obv.valuesSendBuff_)
        , valuesRecvBuff_(obv.valuesRecvBuff_)
        , numFrontIndicesSendBuff_(obv.numFrontIndicesSendBuff_)
        , frontIndicesSendBuff_(obv.frontIndicesSendBuff_)
        , frontIndicesRecvBuff_(obv.frontIndicesRecvBuff_)
        , frontValuesSendBuff_(obv.frontValuesSendBuff_)
        , frontValuesRecvBuff_(obv.frontValuesRecvBuff_)
        , overlap_(obv.overlap_)
    {
    }

    /*!
     * \brief Default constructor.
     */
    OverlappingBlockVector()
    {}

    /*!
     * \brief Assignment operator.
     */
    using ParentType::operator=;
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
    void assignAddBorder(const BlockVector &nbv)
    {
        int numLocal = overlap_->numLocal();
        int numDomestic = overlap_->numDomestic();

        // assign the local rows from the non-overlapping block vector
        const auto &foreignOverlap = overlap_->foreignOverlap();
        for (int localRowIdx = 0; localRowIdx < numLocal; ++localRowIdx)  {
            int nativeRowIdx = foreignOverlap.localToNative(localRowIdx);
            (*this)[localRowIdx] = nbv[nativeRowIdx];
        };

        // set the remote indices to 0 (strictly speaking, that's not
        // necessary because they will be overwritten by the values
        // from their respective master process, but setting them to 0
        // does not cost us much and keeps valgrind from complaining).
        for (int rowIdx = numLocal; rowIdx < numDomestic; ++rowIdx)
            (*this)[rowIdx] = 0.0;

        // add up the contents of border rows, for the remaining rows,
        // get the values from their respective master process.
        syncAddBorder();
    }

    /*!
     * \brief Assign the local values to a non-overlapping block
     *        vector.
     */
    void assignTo(BlockVector &nbv) const
    {
        // assign the local rows
        const auto &foreignOverlap = overlap_->foreignOverlap();
        int numNative = foreignOverlap.numNative();
        for (int nativeRowIdx = 0; nativeRowIdx < numNative; ++nativeRowIdx) {
            int localRowIdx = foreignOverlap.nativeToLocal(nativeRowIdx);

            if (localRowIdx < 0)
                nbv[nativeRowIdx] = 0.0;
            else
                nbv[nativeRowIdx] = (*this)[localRowIdx];
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

        // copy all entries from the master process
        sync();
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
            std::cout << "row " << i << (overlap_->isLocal(i)?" ":"*") << ": " << (*this)[i] << "\n";
        };
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

            const auto &foreignOverlap = overlap_->foreignOverlapWithPeer(peerRank);
            int numEntries = foreignOverlap.size();
            numIndicesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<int> >(new MpiBuffer<int>(1));
            indicesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(new MpiBuffer<Index>(numEntries));
            valuesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numEntries));

            // fill the indices buffer with global indices
            MpiBuffer<Index> &indicesSendBuff = *indicesSendBuff_[peerRank];
            auto ovlpIt = foreignOverlap.begin();
            const auto &ovlpEndIt = foreignOverlap.end();
            for (int i = 0; ovlpIt != ovlpEndIt; ++ovlpIt, ++i) {
                int rowIdx = ovlpIt->index;
                indicesSendBuff[i] = overlap_->domesticToGlobal(rowIdx);
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
            indicesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(new MpiBuffer<Index>(numRows));
            valuesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numRows));
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

            int numFrontEntries = overlap_->foreignOverlap().numFront(peerRank);
            const auto &foreignOverlap = overlap_->foreignOverlapWithPeer(peerRank);
            numFrontIndicesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<int> >(new MpiBuffer<int>(1));
            frontIndicesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(new MpiBuffer<Index>(numFrontEntries));
            frontValuesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numFrontEntries));

            // fill the indices buffer with global indices
            MpiBuffer<Index> &frontIndicesSendBuff = *frontIndicesSendBuff_[peerRank];
            auto ovlpIt = foreignOverlap.begin();
            const auto &ovlpEndIt = foreignOverlap.end();
            int i = 0;
            for (; ovlpIt != ovlpEndIt; ++ovlpIt) {
                int rowIdx = ovlpIt->index;
                if (!overlap_->foreignOverlap().isFrontFor(peerRank, rowIdx))
                    continue;

                frontIndicesSendBuff[i] = overlap_->domesticToGlobal(rowIdx);
                ++i;
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
            frontIndicesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<Index> >(new MpiBuffer<Index>(numFrontRows));
            frontValuesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numFrontRows));
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
        for (unsigned i = 0; i < indices.size(); ++ i) {
            values[i] = (*this)[indices[i]];
        }

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
        for (Index i = 0; i < indices.size(); ++ i)
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
            if (overlap_->masterRank(domRowIdx) == peerRank)
                (*this)[domRowIdx] = values[j];
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
            if (overlap_->isBorderWith(domRowIdx, peerRank)) {
                (*this)[domRowIdx] += values[j];
            }
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
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > frontValuesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > frontValuesRecvBuff_;

    const Overlap *overlap_;
};

} // namespace Linear
} // namespace Dumux

#endif
