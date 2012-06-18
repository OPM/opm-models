// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
#ifndef DUMUX_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH
#define DUMUX_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH

#include "overlaptypes.hh"

#include <dune/grid/common/datahandleif.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <vector>
#include <tuple>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

namespace Dumux {
namespace Linear {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
template<class BCRSMatrix>
class ForeignOverlapFromBCRSMatrix
{
    ForeignOverlapFromBCRSMatrix(const ForeignOverlapFromBCRSMatrix &A)
    {}

public:
    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    ForeignOverlapFromBCRSMatrix(const BCRSMatrix &A,
                                 const BorderList &borderList,
                                 int overlapSize)
        : borderList_(borderList)
    {
        overlapSize_ = overlapSize;

        myRank_ = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
#endif

        numLocal_ = A.N();

        // update the set of indices at the border
        auto it = borderList.begin();
        const auto &endIt = borderList.end();
        for (; it != endIt; ++it)
            borderIndices_.insert(it->localIdx);

        // find the set of processes which have an overlap with the
        // local processes. (i.e. the set of processes which we will
        // have to communicate to.)
        peerSet_.update(borderList);

        // Create an initial seed list of indices which are in the
        // overlap.
        SeedList initialSeedList;
        initialSeedList.update(borderList);

        // calculate the foreign overlap for the local partition,
        // i.e. find the distance of each row from the seed set.
        foreignOverlapByIndex_.resize(A.N());
        extendForeignOverlap_(A, initialSeedList, overlapSize);

        computeMasterRanks_();

        // group foreign overlap by peer process rank
        groupForeignOverlapByRank_();
    }

    /*!
     * \brief Returns the size of the overlap region.
     */
    int overlapSize() const
    { return overlapSize_; }

    /*!
     * \brief Returns true iff a local index is a border index.
     */
    bool isBorder(int localIdx) const
    { return borderIndices_.count(localIdx) > 0; }

    /*!
     * \brief Returns true iff a local index is a border index shared with a given peer process.
     */
    bool isBorderWith(int localIdx, int peerRank) const
    { 
        const auto &indexOverlap = foreignOverlapByIndex_[localIdx];
        const auto &borderDistIt = indexOverlap.find(peerRank);
        if (borderDistIt == indexOverlap.end())
            return false;
                
        // border distance of the index needs to be 0
        return borderDistIt->second == 0;
    }

    /*!
     * \brief Return the rank of the master process of an
     *        index.
     */
    ProcessRank masterRank(int localIdx) const
    { return masterRank_[localIdx]; }

    /*!
     * \brief Return true if the current rank is the "master" of an
     *        index.
     *
     * If the index is at the interior of some process, we define this
     * process as its master, if the index is on the boundary, then
     * the master is defined as the process with the lowest rank.
     */
    bool iAmMasterOf(int localIdx) const
    { return masterRank_[localIdx] == myRank_; }

    /*!
     * \brief Returns the list of indices which intersect the process
     *        border.
     */
    const BorderList &borderList() const
    { return borderList_; }

    /*!
     * \brief Return the list of (local indices, border distance,
     *        number of processes) triples which are in the overlap of
     *        a given peer rank.
     */
    const OverlapWithPeer &foreignOverlapWithPeer(int peerRank) const
    {
        assert(foreignOverlapByRank_.find(peerRank) != foreignOverlapByRank_.end());
        return foreignOverlapByRank_.find(peerRank)->second;
    }

    /*!
     * \brief Return the map of (peer rank, border distance) for a given local index.
     */
    const std::map<ProcessRank, BorderDistance> &foreignOverlapByIndex(int localIdx) const
    {
        assert(isLocal(localIdx));
        return foreignOverlapByIndex_[localIdx];
    }

    /*!
     * \brief Returns true iff a local index is seen by a peer rank.
     */
    bool peerHasIndex(int peerRank, int localIdx) const
    { 
        const auto &idxOverlap = foreignOverlapByIndex_[localIdx];
        return idxOverlap.find(peerRank) != idxOverlap.end();
    }

    /*!
     * \brief Returns the number of front indices of a peer process in
     *        the local partition.
     */
    int numFront(int peerRank) const
    { 
        const auto &peerOverlap = foreignOverlapByRank_.find(peerRank)->second;

        int n = 0;
        auto it = peerOverlap.begin();
        const auto &endIt = peerOverlap.end();
        for (; it != endIt; ++it) {
            if (it->borderDistance == overlapSize_)
                ++ n;
        }
        return n;
    }

    /*!
     * \brief Returns whether a given local index is on the front of a
     *        given peer rank.
     */
    bool isFrontFor(int peerRank, int localIdx) const
    { 
        const auto &idxOverlap = foreignOverlapByIndex_[localIdx];

        auto it = idxOverlap.find(peerRank);
        if (it == idxOverlap.end())
            return false; // index is not in overlap

        return it->second == overlapSize_;
    }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet &peerSet() const
    { return peerSet_; }

    /*!
     * \brief Returns the number local indices
     */
    int numLocal() const
    { return numLocal_; }

    /*!
     * \brief Returns true iff a domestic index is local
     */
    bool isLocal(int domesticIdx) const
    { return domesticIdx < numLocal(); }

    /*!
     * \brief Return the number of peer ranks for which a given local
     *        index is visible.
     */
    int numPeers(int localIdx) const
    { return foreignOverlapByIndex_[localIdx].size(); }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    {
        auto it = foreignOverlapByRank_.begin();
        const auto &endIt = foreignOverlapByRank_.end();
        for (; it != endIt; ++it) {
            std::cout << "Overlap rows(distance) for rank " << it->first << ": ";

            auto rowIt = it->second.begin();
            const auto &rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                std::cout << rowIt->index << "(" << rowIt->borderDistance << ") ";
            };
            std::cout << "\n";
        }
    }

protected:
    // extend the foreign overlaps by one level. this uses a greedy
    // algorithm.
    void extendForeignOverlap_(const BCRSMatrix &A,
                               SeedList &seedList,
                               int overlapSize)
    {
        // add all processes in the seed rows of the current overlap
        // level
        int minOverlapDistance = overlapSize*2;
        auto it = seedList.begin();
        const auto &endIt = seedList.end();
        for (; it != endIt; ++it) {
            int localIdx = it->index;
            int peerRank = it->peerRank;
            int distance = it->borderDistance;
            if (foreignOverlapByIndex_[localIdx].count(peerRank) == 0) {
                foreignOverlapByIndex_[localIdx][peerRank] = distance;
            }
            else {
                foreignOverlapByIndex_[localIdx][peerRank] =
                    std::min(distance,
                             foreignOverlapByIndex_[localIdx][peerRank]);
            }

            minOverlapDistance = std::min(minOverlapDistance, distance);
        }

        // if we have reached the maximum overlap distance, we're
        // finished and break the recursion
        if (minOverlapDistance >= overlapSize)
            return;

        // find the seed list for the next overlap level using the
        // seed set for the current level
        SeedList nextSeedList;
        it = seedList.begin();
        for (; it != endIt; ++it) {
            int rowIdx = it->index;
            int peerRank = it->peerRank;
            int borderDist = it->borderDistance;

            // find all column indies in the row. The indices of the
            // columns are the additional indices of the overlap which
            // we would like to add
            typedef typename BCRSMatrix::ConstColIterator ColIterator;
            ColIterator colIt = A[rowIdx].begin();
            ColIterator colEndIt = A[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int newIdx = colIt.index();

                // if the process is already is in the overlap of the
                // column index, ignore this column index!
                if (foreignOverlapByIndex_[newIdx].count(peerRank) > 0)
                    continue;


                // check whether the new index is already in the overlap
                bool hasIndex = false;
                typename SeedList::iterator sIt = nextSeedList.begin();
                typename SeedList::iterator sEndIt = nextSeedList.end();
                for (; sIt != sEndIt; ++sIt) {
                    if (sIt->index == newIdx &&
                        sIt->peerRank == peerRank)
                    {
                        hasIndex = true;
                        sIt->borderDistance = std::min(sIt->borderDistance, borderDist + 1);
                        break;
                    }
                }
                if (hasIndex)
                    continue; // we already have this index

                // add the current processes to the seed list for the
                // next overlap level
                IndexRankDist newTuple;
                newTuple.index = newIdx;
                newTuple.peerRank = peerRank;
                newTuple.borderDistance =  borderDist + 1;
                nextSeedList.push_back(newTuple);
            }
        }

        // clear the old seed list to save some memory
        seedList.clear();

        // communicate the non-neigbor overlap indices
        addNonNeighborOverlapIndices_(A, nextSeedList, borderDist);

        // Perform the same excercise for the next overlap distance
        extendForeignOverlap_(A,
                              nextSeedList,
                              overlapSize);
    }

    void addNonNeighborOverlapIndices_(const BCRSMatrix &A,
                                       SeedList &seedList,
                                       int borderDist)
    {
        // get all indices in the border which have borderDist as
        // their distance to the closest border of their local process

        // communicate the (peerIdx, nonNeighborPeer, borderDist)
        // triples

        // filter out all indices which are already in the overlap
    };

    // given a list of border indices and provided that
    // borderListToSeedList_() was already called, calculate the
    // master process of each local index.
    void computeMasterRanks_()
    {
        // determine the minimum rank for all indices
        masterRank_.resize(numLocal_);
        for (int localIdx = 0; localIdx < numLocal_; ++localIdx) {
            int masterRank = myRank_;
            if (isBorder(localIdx)) {
                // if the local index is a border index, loop over all ranks
                // for which this index is also a border index. the lowest
                // rank wins!
                auto it = foreignOverlapByIndex_[localIdx].begin();
                const auto &endIt = foreignOverlapByIndex_[localIdx].end();
                for (; it != endIt; ++it) {
                    if (it->second == 0) {
                        // if the border distance is zero, the rank with the minimum
                        masterRank = std::min(masterRank, it->first);
                    }
                }
            }
            masterRank_[localIdx] = masterRank;
        }
    }

    // assuming that the foreign overlap has been created for each
    // local index, this method groups the foreign overlap by peer
    // process rank
    void groupForeignOverlapByRank_()
    {
        // loop over all indices which are in the overlap of some
        // process
        int nIndices = foreignOverlapByIndex_.size();
        for (int i = 0; i < nIndices; ++i)
        {
            // loop over the list of processes for the current index
            auto it = foreignOverlapByIndex_[i].begin();
            const auto &endIt = foreignOverlapByIndex_[i].end();
            int nRanks = foreignOverlapByIndex_[i].size();
            for (; it != endIt; ++it)  {
                IndexDistanceNpeers tmp;
                tmp.index = i;
                tmp.borderDistance = it->second;
                tmp.numPeers = nRanks;
                foreignOverlapByRank_[it->first].push_back(tmp);
            };
        };
    }

    // set of processes with which we have to communicate
    PeerSet peerSet_;

    // the list of indices on the border
    const BorderList &borderList_;

    // an array which contains the rank of the master process for each
    // index
    std::vector<ProcessRank> masterRank_;

    // set of all local indices which are on the border of some remote
    // process
    std::set<Index> borderIndices_;

    // stores the set of process ranks which are in the overlap for a
    // given row index "owned" by the current rank. The second value
    // store the distance from the nearest process border.
    OverlapByIndex foreignOverlapByIndex_;

    // stores a list of foreign overlap indices for each rank
    OverlapByRank foreignOverlapByRank_;

    // size of the overlap region
    BorderDistance overlapSize_;

    // number of local indices
    int numLocal_;

    // the MPI rank of the local process
    int myRank_;
};

} // namespace Linear
} // namespace Dumux

#endif
