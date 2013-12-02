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
                                  const std::set<Index> &blackList,
                                  int overlapSize)
        : foreignOverlap_(A, borderList, blackList, overlapSize),
          globalIndices_(foreignOverlap_)
    {
        myRank_ = 0;

#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
#endif // HAVE_MPI

        buildDomesticOverlap_();
        updateMasterRanks_();
    }

    void check() const
    {
        // check consistency of global indices
        for (int domIdx = 0; domIdx < numDomestic(); ++domIdx) {
            assert(int(globalToDomestic(domesticToGlobal(domIdx))) == domIdx);
        }

        // send the foreign overlap for which we are master to the
        // peers
        std::map<int, MpiBuffer<int> *> sizeBufferMap;

        auto peerIt = peerSet_.begin();
        const auto &peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            auto &buff = *(new MpiBuffer<int>(1));
            sizeBufferMap[*peerIt] = &buff;
            buff[0] = foreignOverlapWithPeer(*peerIt).size();
            buff.send(*peerIt);
        }

        peerIt = peerSet_.begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            MpiBuffer<int> rcvBuff(1);
            rcvBuff.receive(*peerIt);

            assert(rcvBuff[0]
                   == domesticOverlapWithPeer_.find(*peerIt)->second.size());
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
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet &peerSet() const
    { return peerSet_; }

    /*!
     * \brief Return the list of border indices.
     */
    const BorderList &borderList() const
    { return foreignOverlap_.borderList(); }

    /*!
     * \brief Returns true iff a domestic index is a border index.
     */
    bool isBorder(int domesticIdx) const
    { return isLocal(domesticIdx) && foreignOverlap_.isBorder(domesticIdx); }

    /*!
     * \brief Returns true iff a domestic index is on the border with
     *        a given peer process.
     */
    bool isBorderWith(int domesticIdx, int peerRank) const
    {
        return isLocal(domesticIdx)
               && foreignOverlap_.isBorderWith(domesticIdx, peerRank);
    }

    /*!
     * \brief Returns true iff a domestic index is on the front.
     */
    bool isFront(int domesticIdx) const
    {
        if (isLocal(domesticIdx))
            return false;

        // check wether the border distance of the domestic overlap is
        // maximal for the index
        const auto &domOverlap = domesticOverlapByIndex_[domesticIdx];
        return domOverlap.size() > 0 && int(domOverlap.begin()->second)
                                        == foreignOverlap_.overlapSize();
    }

    /*!
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    int numPeers(int domesticIdx) const
    { return domesticOverlapByIndex_[domesticIdx].size(); }

    /*!
     * \brief Returns the size of the overlap region
     */
    int overlapSize() const
    { return foreignOverlap_.overlapSize(); }

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
    { return domesticIdx < numLocal(); }

    /*!
     * \brief Return true iff the current process is the master of a
     *        given domestic index.
     */
    bool iAmMasterOf(int domesticIdx) const
    {
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.iAmMasterOf(domesticIdx);
    }

    /*!
     * \brief Return the rank of a master process for a domestic index
     */
    int masterRank(int domesticIdx) const
    { return masterRank_[domesticIdx]; }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    { globalIndices_.print(); }

    /*!
     * \brief Returns a domestic index given a global one
     */
    Index globalToDomestic(Index globalIdx) const
    { return globalIndices_.globalToDomestic(globalIdx); }

    /*!
     * \brief Returns a global index given a domestic one
     */
    Index domesticToGlobal(Index domIdx) const
    { return globalIndices_.domesticToGlobal(domIdx); }

    /*!
     * \brief Returns the foreign overlap of the other processes with
     *        the local process.
     */
    const ForeignOverlap &foreignOverlap() const
    { return foreignOverlap_; }

    /*!
     * \brief Returns true if a given domestic index is either in the foreign or
     * in the domestic overlap.
     */
    bool isInOverlap(Index domesticIdx) const
    {
        return !this->isLocal(domesticIdx)
               || this->foreignOverlap_.isInOverlap(domesticIdx);
    }

    /*!
     * \brief Returns the foreign overlap of a peer process with the
     *        local process.
     */
    const OverlapWithPeer &foreignOverlapWithPeer(ProcessRank peerRank) const
    { return foreignOverlap_.foreignOverlapWithPeer(peerRank); }

    /*!
     * \brief Returns true iff a local index is seen by a peer rank.
     */
    bool peerHasIndex(int peerRank, int localIdx) const
    { return foreignOverlap_.peerHasIndex(peerRank, localIdx); }

    /*!
     * \brief Returns the domestic overlap of a local process with a
     *        remote process.
     */
    const DomesticOverlapWithPeer &domesticOverlapWithPeer(ProcessRank peerRank) const
    {
        assert(domesticOverlapWithPeer_.find(peerRank)
               != domesticOverlapWithPeer_.end());
        return domesticOverlapWithPeer_.find(peerRank)->second;
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
        PeerSet::const_iterator peerEndIt = foreignOverlap_.peerSet().end();

        // send the overlap indices to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendIndicesToPeer_(peerRank);
        };

        // receive our overlap from the processes to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveIndicesFromPeer_(peerRank);
        };

        // receive our overlap from the processes to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            waitSendIndices_(peerRank);
        };
    }

    void updateMasterRanks_()
    {
        unsigned nLocal = numLocal();
        unsigned nDomestic = numDomestic();
        masterRank_.resize(nDomestic);

        // take the master ranks for the local indices from the
        // foreign overlap
        for (size_t i = 0; i < nLocal; ++i)
            masterRank_[i] = foreignOverlap_.masterRank(i);

        // for non-local indices, initialy use INT_MAX as their master
        // rank
        for (size_t i = nLocal; i < nDomestic; ++i)
            masterRank_[i] = std::numeric_limits<ProcessRank>::max();

        // for the non-local indices, take the peer process for which
        // a given local index is in the interior
        auto peerIt = peerSet_.begin();
        const auto &peerEndIt = peerSet_.end();
        for (; peerIt != peerEndIt; ++peerIt) {
            const auto &overlapWithPeer
                = domesticOverlapWithPeer_.find(*peerIt)->second;

            auto idxIt = overlapWithPeer.begin();
            const auto &idxEndIt = overlapWithPeer.end();
            for (; idxIt != idxEndIt; ++idxIt) {
                if (isLocal(*idxIt))
                    continue; // ignore border indices

                masterRank_[*idxIt] = std::min(masterRank_[*idxIt], *peerIt);
            }
        }
    }

    void sendIndicesToPeer_(int peerRank)
    {
#if HAVE_MPI
        const auto &foreignOverlap
            = foreignOverlap_.foreignOverlapWithPeer(peerRank);

        // first, send a message containing the number of additional
        // indices stemming from the overlap (i.e. without the border
        // indices)
        int numIndices = foreignOverlap.size();
        numIndicesSendBuff_[peerRank] = new MpiBuffer<size_t>(1);
        (*numIndicesSendBuff_[peerRank])[0] = numIndices;
        numIndicesSendBuff_[peerRank]->send(peerRank);

        // create MPI buffers
        indicesSendBuff_[peerRank]
            = new MpiBuffer<IndexDistanceNpeers>(numIndices);

        // then send the additional indices themselfs
        auto overlapIt = foreignOverlap.begin();
        const auto &overlapEndIt = foreignOverlap.end();
        for (int i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            int localIdx = overlapIt->index;
            int borderDistance = overlapIt->borderDistance;
            const auto &foreignIndexOverlap
                = foreignOverlap_.foreignOverlapByIndex(localIdx);
            int numPeers = foreignIndexOverlap.size();

            IndexDistanceNpeers tmp;
            tmp.index = globalIndices_.domesticToGlobal(localIdx);
            tmp.borderDistance = borderDistance;
            tmp.numPeers = numPeers;

            (*indicesSendBuff_[peerRank])[i] = tmp;
        };

        indicesSendBuff_[peerRank]->send(peerRank);

#endif // HAVE_MPI
    }

    void waitSendIndices_(int peerRank)
    {
        numIndicesSendBuff_[peerRank]->wait();
        delete numIndicesSendBuff_[peerRank];

        indicesSendBuff_[peerRank]->wait();
        delete indicesSendBuff_[peerRank];
    }

    void receiveIndicesFromPeer_(int peerRank)
    {
#if HAVE_MPI
        // receive the number of additional indices
        size_t numIndices = -1;
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
            assert(!(borderDistance == 0 && !isLocal(domesticIdx)));
            assert(!(borderDistance > 0 && isLocal(domesticIdx)));

            borderDistance_[domesticIdx]
                = std::min(borderDistance, borderDistance_[domesticIdx]);
        }
#endif // HAVE_MPI
    }

    int myRank_;
    ForeignOverlap foreignOverlap_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    OverlapByIndex domesticOverlapByIndex_;
    std::vector<BorderDistance> borderDistance_;
    std::vector<Index> masterRank_;

    std::map<ProcessRank, MpiBuffer<size_t> *> numIndicesSendBuff_;
    std::map<ProcessRank, MpiBuffer<IndexDistanceNpeers> *> indicesSendBuff_;
    GlobalIndices globalIndices_;
    PeerSet peerSet_;
};

} // namespace Linear
} // namespace Ewoms

#endif
