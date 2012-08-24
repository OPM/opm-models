// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \brief This class maps domestic row indices to and from "global"
 *        indices which is used to construct an algebraic overlap
 *        for the parallel linear solvers.
 */
#ifndef DUMUX_GLOBAL_INDICES_HH
#define DUMUX_GLOBAL_INDICES_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <tuple>

#if HAVE_MPI
#include <mpi.h>
#endif

#include "overlaptypes.hh"

namespace Dumux {
namespace Linear {
/*!
 * \brief This class maps domestic row indices to and from "global"
 *        indices which is used to construct an algebraic overlap
 *        for the parallel linear solvers.
 */
template <class ForeignOverlap>
class GlobalIndices
{
    GlobalIndices(const GlobalIndices &A)
    {}

    typedef std::map<Index, Index> GlobalToDomesticMap;
    typedef std::map<Index, Index> DomesticToGlobalMap;

public:
    GlobalIndices(const ForeignOverlap &foreignOverlap)
        : foreignOverlap_(foreignOverlap)
    {
        myRank_ = 0;
        mpiSize_ = 1;

#if HAVE_MPI
        {
            int tmp;
            MPI_Comm_rank(MPI_COMM_WORLD, &tmp); myRank_ = tmp;
            MPI_Comm_size(MPI_COMM_WORLD, &tmp); mpiSize_ = tmp;
        }
#endif

        // calculate the domestic overlap (i.e. all overlap indices in
        // foreign processes which the current process overlaps.)
        // This requires communication via MPI.
        buildGlobalIndices_();
    }

    /*!
     * \brief Converts a domestic index to a global one.
     */
    int domesticToGlobal(int domesticIdx) const
    {
        assert(domesticToGlobal_.find(domesticIdx) != domesticToGlobal_.end());

        return domesticToGlobal_.find(domesticIdx)->second;
    }

    /*!
     * \brief Converts a global index to a domestic one.
     */
    int globalToDomestic(int globalIdx) const
    {
        assert(globalToDomestic_.find(globalIdx) != globalToDomestic_.end());

        return globalToDomestic_.find(globalIdx)->second;
    }

    /*!
     * \brief Returns the number of indices which are in the interior or
     *        on the border of the current rank.
     */
    int numLocal() const
    { return foreignOverlap_.numLocal(); }

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its copies of indices in the overlap regions
     */
    int numDomestic() const
    {
        return numDomestic_;
    }

    /*!
     * \brief Add an index to the domestic<->global mapping.
     */
    void addIndex(int domesticIdx, int globalIdx)
    {
        domesticToGlobal_[domesticIdx] = globalIdx;
        globalToDomestic_[globalIdx] = domesticIdx;
        numDomestic_ = domesticToGlobal_.size();

        assert(domesticToGlobal_.size() == globalToDomestic_.size());
    }

    /*!
     * \brief Send a border index to a remote process.
     */
    void sendBorderIndex(int peerRank, int domesticIdx, int peerLocalIdx)
    {
#if HAVE_MPI
        PeerIndexGlobalIndex sendBuf;
        sendBuf.peerIdx = peerLocalIdx;
        sendBuf.globalIdx = domesticToGlobal(domesticIdx);
        MPI_Send(&sendBuf, // buff
                 sizeof(PeerIndexGlobalIndex), // count
                 MPI_BYTE, // data type
                 peerRank,
                 0, // tag
                 MPI_COMM_WORLD); // communicator
#endif // HAVE_MPI
    }

    /*!
     * \brief Receive an index on the border from a remote
     *        process and add it the translation maps.
     */
    void receiveBorderIndex(int peerRank)
    {
#if HAVE_MPI
        PeerIndexGlobalIndex recvBuf;
        MPI_Recv(&recvBuf, // buff
                 sizeof(PeerIndexGlobalIndex), // count
                 MPI_BYTE, // data type
                 peerRank,
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int domesticIdx = recvBuf.peerIdx;
        int globalIdx = recvBuf.globalIdx;
        addIndex(domesticIdx, globalIdx);
#endif // HAVE_MPI
    }

    /*!
     * \brief Return true iff a given global index already exists
     */
    bool hasGlobalIndex(int globalIdx) const
    { return globalToDomestic_.find(globalIdx) != globalToDomestic_.end(); }

    /*!
     * \brief Prints the global indices of all domestic indices
     *        for debugging purposes.
     */
    void print() const
    {
        std::cout << "(domestic index, global index, domestic->global->domestic) list for rank " <<
            myRank_ << "\n";

        for (int domIdx = 0; domIdx < domesticToGlobal_.size(); ++ domIdx) {
            std::cout << "(" << domIdx
                      << ", " << domesticToGlobal(domIdx)
                      << ", " << globalToDomestic(domesticToGlobal(domIdx))
                      << ") ";
        };
        std::cout << "\n";
    }

protected:
    // retrieve the offset for the indices where we are master in the
    // global index list
    void buildGlobalIndices_()
    {
#if HAVE_MPI
        numDomestic_ = 0;
#else
        numDomestic_ = foreignOverlap_.numLocal();
#endif

#if HAVE_MPI
        if (myRank_ == 0) {
            // the first rank starts at index zero
            domesticOffset_ = 0;
        }
        else {
            // all other ranks retrieve their offset from the next
            // lower rank
            MPI_Recv(&domesticOffset_, // buffer
                     1, // count
                     MPI_INT, // data type
                     myRank_ - 1,
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);
        }

        // create maps for all indices for which the current process
        // is the master
        int numMaster = 0;
        for (int i = 0; i < foreignOverlap_.numLocal(); ++i) {
            if (!foreignOverlap_.iAmMasterOf(i))
                continue;

            addIndex(i, domesticOffset_ + numMaster);
            ++ numMaster;
        }

        if (myRank_ < mpiSize_ - 1) {
            // send the domestic offset plus the number of master
            // indices to the process which is one rank higher
            int tmp = domesticOffset_ + numMaster;
            MPI_Send(&tmp, // buff
                     1, // count
                     MPI_INT, // data type
                     myRank_ + 1, // peer rank
                     0, // tag
                     MPI_COMM_WORLD); // communicator
        };

        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = peerSet_().end();
        // receive the border indices from the lower ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt < myRank_)
                receiveBorderFrom_(*peerIt);
        }

        // send the border indices to the higher ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt > myRank_)
                sendBorderTo_(*peerIt);
        }

        // receive the border indices from the higher ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt > myRank_)
                receiveBorderFrom_(*peerIt);
        }

        // send the border indices to the lower ranks
        peerIt = peerSet_().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            if (*peerIt < myRank_)
                sendBorderTo_(*peerIt);
        }
#endif // HAVE_MPI
    }

    void sendBorderTo_(ProcessRank peerRank)
    {
#if HAVE_MPI
        // send (local index on myRank, global index) pairs to the
        // peers
        BorderList::const_iterator borderIt = borderList_().begin();
        BorderList::const_iterator borderEndIt = borderList_().end();
        for (; borderIt != borderEndIt; ++borderIt) {
            ProcessRank borderPeer = borderIt->peerRank;
            if (borderPeer != peerRank)
                continue;

            Index localIdx = borderIt->localIdx;
            Index peerIdx = borderIt->peerIdx;
            if (foreignOverlap_.iAmMasterOf(borderIt->localIdx)) {
                sendBorderIndex(borderPeer, localIdx, peerIdx);
            }
        }
#endif // HAVE_MPI
    }

    void receiveBorderFrom_(ProcessRank peerRank)
    {
#if HAVE_MPI
        // retrieve the global indices for which we are not master
        // from the processes with lower rank
        BorderList::const_iterator borderIt = borderList_().begin();
        BorderList::const_iterator borderEndIt = borderList_().end();
        for (; borderIt != borderEndIt; ++borderIt) {
            ProcessRank borderPeer = borderIt->peerRank;
            if (borderPeer != peerRank)
                continue;

            if (foreignOverlap_.masterRank(borderIt->localIdx) == borderPeer) {
                receiveBorderIndex(borderPeer);
            }
        }
#endif // HAVE_MPI
    }

    const PeerSet &peerSet_() const
    { return foreignOverlap_.peerSet(); }

    const BorderList &borderList_() const
    { return foreignOverlap_.borderList(); }

    ProcessRank myRank_;
    size_t mpiSize_;

    size_t domesticOffset_;
    size_t numDomestic_;
    const ForeignOverlap &foreignOverlap_;

    GlobalToDomesticMap globalToDomestic_;
    DomesticToGlobalMap domesticToGlobal_;
};

} // namespace Linear
} // namespace Dumux

#endif
