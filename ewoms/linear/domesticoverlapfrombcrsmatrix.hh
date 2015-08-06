// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \copydoc Ewoms::Linear::DomesticOverlapFromBCRSMatrix
 */
#ifndef EWOMS_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH
#define EWOMS_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH

#include "foreignoverlapfrombcrsmatrix.hh"
#include "blacklist.hh"
#include "globalindices.hh"

#include <ewoms/parallel/mpibuffer.hh>

#include <algorithm>
#include <limits>
#include <set>
#include <map>
#include <vector>

namespace Ewoms {
namespace Linear {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
template <class BCRSMatrix>
class DomesticOverlapFromBCRSMatrix
{
    DomesticOverlapFromBCRSMatrix(const DomesticOverlapFromBCRSMatrix &A)
    {}

    typedef Ewoms::Linear::ForeignOverlapFromBCRSMatrix<BCRSMatrix> ForeignOverlap;
    typedef Ewoms::Linear::GlobalIndices<ForeignOverlap> GlobalIndices;

public:
    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    DomesticOverlapFromBCRSMatrix(const BCRSMatrix &A,
                                  const BorderList &borderList,
                                  const BlackList &blackList,
                                  int overlapSize)
        : foreignOverlap_(A, borderList, blackList, overlapSize)
        , blackList_(blackList)
        , globalIndices_(foreignOverlap_)
    {
        myRank_ = 0;
        worldSize_ = 1;

#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize_);
#endif // HAVE_MPI

        buildDomesticOverlap_();
        updateMasterRanks_();
        blackList_.updateNativeToDomesticMap(*this);

        setupDebugMapping_();
    }

    void check() const
    {
#ifndef NDEBUG
        // check consistency of global indices
        for (int domIdx = 0; domIdx < numDomestic(); ++domIdx) {
            assert(int(globalToDomestic(domesticToGlobal(domIdx))) == domIdx);
        }
#endif // NDEBUG

        // send the foreign overlap for which we are master to the
        // peers
        std::map<int, MpiBuffer<int> *> sizeBufferMap;

        auto peerIt = peerSet_.begin();
        const auto &peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            auto &buffer = *(new MpiBuffer<int>(1));
            sizeBufferMap[*peerIt] = &buffer;
            buffer[0] = foreignOverlap_.foreignOverlapWithPeer(*peerIt).size();
            buffer.send(*peerIt);
        }

        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            MpiBuffer<int> rcvBuffer(1);
            rcvBuffer.receive(*peerIt);

            assert(rcvBuffer[0] == domesticOverlapWithPeer_.find(*peerIt)->second.size());
        }

        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            sizeBufferMap[*peerIt]->wait();
            delete sizeBufferMap[*peerIt];
        }
    }

    /*!
     * \brief Returns the rank of the current process.
     */
    int myRank() const
    { return myRank_; }

    /*!
     * \brief Returns the number of processes in the global MPI communicator.
     */
    int worldSize() const
    { return worldSize_; }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet &peerSet() const
    { return peerSet_; }

    /*!
     * \brief Returns the number of indices on the border with a given
     *        peer rank.
     */
    int numBorder(int peerRank) const
    { return foreignOverlap_.numBorder(peerRank); }

    /*!
     * \brief Returns true iff a domestic index is a border index.
     */
    bool isBorder(int domesticIdx) const
    {
        return isLocal(domesticIdx)
            && foreignOverlap_.isBorder(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns true iff a domestic index is on the border with
     *        a given peer process.
     */
    bool isBorderWith(int domesticIdx, int peerRank) const
    {
        return isLocal(domesticIdx)
            && foreignOverlap_.isBorderWith(mapExternalToInternal_(domesticIdx),
                                            peerRank);
    }

    /*!
     * \brief Returns the number of indices on the front within a given
     *        peer rank's grid partition.
     */
    int numFront(int peerRank) const
    { return foreignOverlap_.numFront(peerRank); }

    /*!
     * \brief Returns true iff a domestic index is on the front.
     */
    bool isFront(int domesticIdx) const
    {
        if (isLocal(domesticIdx))
            return false;
        Index internalDomesticIdx = mapExternalToInternal_(domesticIdx);

        // check wether the border distance of the domestic overlap is
        // maximal for the index
        const auto &domOverlap = domesticOverlapByIndex_[internalDomesticIdx];
        return domOverlap.size() > 0
            && int(domOverlap.begin()->second) == foreignOverlap_.overlapSize();
    }

    /*!
     * \brief Returns the object which represents the black-listed native indices.
     */
    const BlackList& blackList() const
    { return blackList_; }

    /*!
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    int numPeers(int domesticIdx) const
    { return domesticOverlapByIndex_[mapExternalToInternal_(domesticIdx)].size(); }

    /*!
     * \brief Returns the size of the overlap region
     */
    int overlapSize() const
    { return foreignOverlap_.overlapSize(); }

    /*!
     * \brief Returns the number native indices
     *
     * I.e. the number of indices of the "raw" grid partition of the
     * local process (including the indices in ghost and overlap
     * elements).
     */
    int numNative() const
    { return foreignOverlap_.numNative(); }

    /*!
     * \brief Returns the number local indices
     *
     * I.e. indices in the interior or on the border of the process'
     * domain.
     */
    int numLocal() const
    { return foreignOverlap_.numLocal(); }

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its domestic overlap (i.e. indices for which it is not
     * neither master nor are on the process border).
     */
    int numDomestic() const
    { return globalIndices_.numDomestic(); }

    /*!
     * \brief Return true if a domestic index is local for the process
     *
     * I.e. the entity for this index is in the interior or on the
     * border of the process' domain.
     */
    bool isLocal(int domesticIdx) const
    { return mapExternalToInternal_(domesticIdx) < numLocal(); }

    /*!
     * \brief Return true iff the current process is the master of a
     *        given domestic index.
     */
    bool iAmMasterOf(int domesticIdx) const
    {
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.iAmMasterOf(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Return the rank of a master process for a domestic index
     */
    int masterRank(int domesticIdx) const
    { return masterRank_[mapExternalToInternal_(domesticIdx)]; }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    { globalIndices_.print(); }

    /*!
     * \brief Returns a domestic index given a global one
     */
    Index globalToDomestic(Index globalIdx) const
    {
        Index internalIdx = globalIndices_.globalToDomestic(globalIdx);
        if (internalIdx < 0)
            return -1;
        return mapInternalToExternal_(internalIdx);
    }

    /*!
     * \brief Returns a global index given a domestic one
     */
    Index domesticToGlobal(Index domIdx) const
    { return globalIndices_.domesticToGlobal(mapExternalToInternal_(domIdx)); }

    /*!
     * \brief Returns a native index given a domestic one
     */
    Index domesticToNative(Index domIdx) const
    {
        Index internalIdx = mapExternalToInternal_(domIdx);
        if (internalIdx >= numLocal())
            return -1;
        return foreignOverlap_.localToNative(internalIdx);
    }

    /*!
     * \brief Returns a domestic index given a native one
     */
    Index nativeToDomestic(Index nativeIdx) const
    {
        Index localIdx = foreignOverlap_.nativeToLocal(nativeIdx);
        if (localIdx < 0)
            return localIdx;
        return mapInternalToExternal_(localIdx);
    }

    /*!
     * \brief Returns true if a given domestic index is either in the
     *        foreign or in the domestic overlap.
     */
    bool isInOverlap(Index domesticIdx) const
    {
        return !this->isLocal(domesticIdx)
               || this->foreignOverlap_.isInOverlap(mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns true if a given domestic index is a front index
     *        for a peer rank.
     */
    bool isFrontFor(ProcessRank peerRank, Index domesticIdx) const
    {
        int internalIdx = mapExternalToInternal_(domesticIdx);
        return this->foreignOverlap_.isFrontFor(peerRank, internalIdx);
    }

    /*!
     * \brief Returns true iff a domestic index is seen by a peer rank.
     */
    bool peerHasIndex(int peerRank, int domesticIdx) const
    {
        return foreignOverlap_.peerHasIndex(peerRank,
                                            mapExternalToInternal_(domesticIdx));
    }

    /*!
     * \brief Returns number of indices which are contained in the
     *        foreign overlap with a peer.
     */
    int foreignOverlapSize(ProcessRank peerRank) const
    { return foreignOverlap_.foreignOverlapWithPeer(peerRank).size(); }

    /*!
     * \brief Returns the domestic index given an offset in the
     *        foreign overlap of a peer process with the local
     *        process.
     */
    int foreignOverlapOffsetToDomesticIdx(ProcessRank peerRank, int overlapOffset) const
    {
        int internalIdx =
            foreignOverlap_.foreignOverlapWithPeer(peerRank)[overlapOffset].index;
        return mapInternalToExternal_(internalIdx);
    }

    /*!
     * \brief Returns number of indices which are contained in the
     *        domestic overlap with a peer.
     */
    int domesticOverlapSize(ProcessRank peerRank) const
    { return domesticOverlapWithPeer_.at(peerRank).size(); }

    /*!
     * \brief Returns the domestic index given an offset in the
     *        domestic overlap of a peer process with the local
     *        process.
     */
    int domesticOverlapOffsetToDomesticIdx(ProcessRank peerRank, int overlapOffset) const
    {
        int internalIdx = domesticOverlapWithPeer_.at(peerRank)[overlapOffset];
        return mapInternalToExternal_(internalIdx);
    }

protected:
    void buildDomesticOverlap_()
    {
        // copy the set of peers from the foreign overlap
        peerSet_ = foreignOverlap_.peerSet();

        // resize the array which stores the number of peers for
        // each entry.
        domesticOverlapByIndex_.resize(numLocal());
        borderDistance_.resize(numLocal(), 0);

        PeerSet::const_iterator peerIt;
        PeerSet::const_iterator peerEndIt = peerSet_.end();

        // send the overlap indices to all peer processes
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendIndicesToPeer_(peerRank);
        }

        // receive our overlap from the processes to all peer processes
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveIndicesFromPeer_(peerRank);
        }

        // wait until all send operations complete
        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            waitSendIndices_(peerRank);
        }
    }

    void updateMasterRanks_()
    {
        unsigned nLocal = numLocal();
        unsigned nDomestic = numDomestic();
        masterRank_.resize(nDomestic);

        // take the master ranks for the local indices from the
        // foreign overlap
        for (size_t i = 0; i < nLocal; ++i) {
            masterRank_[i] = foreignOverlap_.masterRank(i);
        }

        // for non-local indices, initially use INT_MAX as their master
        // rank
        for (size_t i = nLocal; i < nDomestic; ++i)
            masterRank_[i] = std::numeric_limits<ProcessRank>::max();

        // for the non-local indices, take the peer process for which
        // a given local index is in the interior
        auto peerIt = peerSet_.begin();
        const auto &peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            const auto &overlapWithPeer = domesticOverlapWithPeer_.find(*peerIt)->second;

            auto idxIt = overlapWithPeer.begin();
            const auto &idxEndIt = overlapWithPeer.end();
            for (; idxIt != idxEndIt; ++idxIt) {
                if (*idxIt >= 0 && foreignOverlap_.isLocal(*idxIt))
                    continue; // ignore border indices

                masterRank_[*idxIt] = std::min(masterRank_[*idxIt], *peerIt);
            }
        }
    }

    void sendIndicesToPeer_(int peerRank)
    {
#if HAVE_MPI
        const auto &foreignOverlap = foreignOverlap_.foreignOverlapWithPeer(peerRank);

        // first, send a message containing the number of additional
        // indices stemming from the overlap (i.e. without the border
        // indices)
        int numIndices = foreignOverlap.size();
        numIndicesSendBuffer_[peerRank] = new MpiBuffer<size_t>(1);
        (*numIndicesSendBuffer_[peerRank])[0] = numIndices;
        numIndicesSendBuffer_[peerRank]->send(peerRank);

        // create MPI buffers
        indicesSendBuffer_[peerRank] = new MpiBuffer<IndexDistanceNpeers>(numIndices);

        // then send the additional indices themselfs
        auto overlapIt = foreignOverlap.begin();
        const auto &overlapEndIt = foreignOverlap.end();
        for (int i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            int localIdx = overlapIt->index;
            int borderDistance = overlapIt->borderDistance;
            int numPeers = foreignOverlap_.foreignOverlapByLocalIndex(localIdx).size();

            IndexDistanceNpeers tmp;
            tmp.index = globalIndices_.domesticToGlobal(localIdx);
            tmp.borderDistance = borderDistance;
            tmp.numPeers = numPeers;

            (*indicesSendBuffer_[peerRank])[i] = tmp;
        }

        indicesSendBuffer_[peerRank]->send(peerRank);
#endif // HAVE_MPI
    }

    void waitSendIndices_(int peerRank)
    {
        numIndicesSendBuffer_[peerRank]->wait();
        delete numIndicesSendBuffer_[peerRank];

        indicesSendBuffer_[peerRank]->wait();
        delete indicesSendBuffer_[peerRank];
    }

    void receiveIndicesFromPeer_(int peerRank)
    {
#if HAVE_MPI
        // receive the number of additional indices
        int numIndices = -1;
        MpiBuffer<size_t> numIndicesRecvBuff(1);
        numIndicesRecvBuff.receive(peerRank);
        numIndices = numIndicesRecvBuff[0];

        // receive the additional indices themselfs
        MpiBuffer<IndexDistanceNpeers> recvBuff(numIndices);
        recvBuff.receive(peerRank);
        for (Index i = 0; i < numIndices; ++i) {
            Index globalIdx = recvBuff[i].index;
            BorderDistance borderDistance = recvBuff[i].borderDistance;

            // if the index is not already known, add it to the
            // domestic indices
            if (!globalIndices_.hasGlobalIndex(globalIdx)) {
                Index newDomesticIdx = globalIndices_.numDomestic();
                globalIndices_.addIndex(newDomesticIdx, globalIdx);

                size_t newSize = globalIndices_.numDomestic();
                borderDistance_.resize(newSize, std::numeric_limits<int>::max());
                domesticOverlapByIndex_.resize(newSize);
            }

            // convert the global index into a domestic one
            Index domesticIdx = globalIndices_.globalToDomestic(globalIdx);

            // extend the domestic overlap
            domesticOverlapByIndex_[domesticIdx][peerRank] = borderDistance;
            domesticOverlapWithPeer_[peerRank].push_back(domesticIdx);

            assert(borderDistance >= 0);
            assert(globalIdx >= 0);
            assert(domesticIdx >= 0);
            assert(!(borderDistance == 0 && !foreignOverlap_.isLocal(domesticIdx)));
            assert(!(borderDistance > 0 && foreignOverlap_.isLocal(domesticIdx)));

            borderDistance_[domesticIdx] = std::min(borderDistance, borderDistance_[domesticIdx]);
        }
#endif // HAVE_MPI
    }

    // this method is intended to set up the code mapping code for
    // mapping domestic indices to the same ones used by a sequential
    // grid. this requires detailed knowledge about how a grid
    // distributes the degrees of freedom over multiple processes, but
    // it can simplify debugging considerably because the indices can
    // be made identical for the parallel and the sequential
    // computations.
    //
    // by default, this method does nothing
    void setupDebugMapping_()
    {}

    // this method is intended to map domestic indices to the ones
    // used by a sequential grid.
    //
    // by default, this method does nothing
    Index mapInternalToExternal_(Index internalIdx) const
    { return internalIdx; }

    // this method is intended to map the indices used by a sequential
    // to grid domestic indices ones.
    //
    // by default, this method does nothing
    Index mapExternalToInternal_(Index externalIdx) const
    { return externalIdx; }

    int myRank_;
    int worldSize_;
    ForeignOverlap foreignOverlap_;

    BlackList blackList_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    OverlapByIndex domesticOverlapByIndex_;
    std::vector<BorderDistance> borderDistance_;
    std::vector<ProcessRank> masterRank_;

    std::map<ProcessRank, MpiBuffer<size_t> *> numIndicesSendBuffer_;
    std::map<ProcessRank, MpiBuffer<IndexDistanceNpeers> *> indicesSendBuffer_;
    GlobalIndices globalIndices_;
    PeerSet peerSet_;
};

} // namespace Linear
} // namespace Ewoms

#endif
