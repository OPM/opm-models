// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief This files provides several data structures for storing
 *        tuples of indices of remote and/or local processes.
 */
#ifndef DUMUX_OVERLAP_TYPES_HH
#define DUMUX_OVERLAP_TYPES_HH

#include <set>
#include <list>
#include <vector>
#include <map>
#include <cstddef>

namespace Dumux {

namespace Linear {

/*!
 * \brief The type of an index of a degree of freedom.
 */
typedef size_t Index;

/*!
 * \brief The type of the rank of a process.
 */
typedef unsigned ProcessRank;

/*!
 * \brief The type representing the distance of an index to the border.
 */
typedef unsigned BorderDistance;

/*!
 * \brief This structure stores an index and a process rank
 */
struct IndexRank {
    Index index;
    ProcessRank rank;
};

/*!
 * \brief This structure stores a local index on a peer process and a
 *        global index.
 */
struct PeerIndexGlobalIndex {
    Index peerIdx;
    Index globalIdx;
};

/*!
 * \brief This structure stores an index, a process rank, and the
 *        distance of the degree of freedom to the process border.
 */
struct IndexRankDist {
    Index index;
    ProcessRank peerRank;
    BorderDistance borderDistance;
};

/*!
 * \brief This structure stores an index, a process rank, and the
 *        number of processes which "see" the degree of freedom with
 *        the index.
 */
struct IndexDistanceNpeers
{
    Index index;
    BorderDistance borderDistance;
    int numPeers;
};

/*!
 * \brief A single index intersecting with the process boundary.
 */
struct BorderIndex
{
    //! Index of the entity for the local process
    Index localIdx;

    //! Index of the entity for the peer process
    Index peerIdx;

    //! Rank of the peer process
    ProcessRank peerRank;

    //! Distance to the process border for the peer (in hops)
    BorderDistance borderDistance;
};

/*!
 * \brief This class managages a list of indices which are on the
 *        border of a process' partition of the grid
 */
typedef std::list<BorderIndex> BorderList;

/*!
 * \brief The list of indices which are on the process boundary.
 */
class SeedList : public std::list<IndexRankDist>
{
public:
    void update(const BorderList &borderList)
    {
        this->clear();

        auto it = borderList.begin();
        const auto &endIt = borderList.end();
        for (; it != endIt; ++it) {
            IndexRankDist ird;
            ird.index = it->localIdx;
            ird.peerRank = it->peerRank;
            ird.borderDistance = it->borderDistance;

            this->push_back(ird);
        };
    }
};

/*!
 * \brief A set of process ranks
 */
class PeerSet : public std::set<Index>
{
public:
    void update(const BorderList &borderList)
    {
        this->clear();

        auto it = borderList.begin();
        const auto &endIt = borderList.end();
        for (; it != endIt; ++it) {
            this->insert(it->peerRank);
        };
    }
};

/*!
 * \brief The list of indices which overlap with a peer rank.
 */
typedef std::vector<IndexDistanceNpeers> OverlapWithPeer;

/*!
 * \brief A type mapping the process rank to the list of indices
 *        shared with this peer.
 */
typedef std::map<ProcessRank, OverlapWithPeer> OverlapByRank;

/*!
 * \brief Maps each index to a list of processes .
 */
typedef std::vector<std::map<ProcessRank, BorderDistance> > OverlapByIndex;

/*!
 * \brief The list of domestic indices are owned by peer rank.
 */
typedef std::vector<Index> DomesticOverlapWithPeer;

/*!
 * \brief A type mapping the process rank to the list of domestic indices
 *        which are owned by the peer.
 */
typedef std::map<ProcessRank, DomesticOverlapWithPeer> DomesticOverlapByRank;

} // namespace Linear
} // namespace Dumux

#endif
