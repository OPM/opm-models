// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
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
 * \brief Simplifies handling of buffers to be used in conjunction with MPI
 */
#ifndef DUMUX_MPI_BUFFER_HH
#define DUMUX_MPI_BUFFER_HH

#if HAVE_MPI
#include <mpi.h>
#endif

#include <type_traits>
#include <cassert>

namespace Dumux
{
/*!
 * \brief Simplifies handling of buffers to be used in conjunction with MPI
 */
template <class DataType>
class MpiBuffer
{
public:
    MpiBuffer()
    {
        data_ = nullptr;
        dataSize_ = 0;
        
        setMpiDataType_();
        updateMpiDataSize_();
    }

    MpiBuffer(int size)
    {
        data_ = new DataType[size];
        dataSize_ = size;

        setMpiDataType_();
        updateMpiDataSize_();
    }

    ~MpiBuffer()
    {
        delete[] data_;
    }

    /*!
     * \brief Set the size of the buffer
     */
    void resize(size_t newSize)
    {
        delete[] data_;
        data_ = new DataType[newSize];
        dataSize_ = newSize;
        updateMpiDataSize_();
    }

    /*!
     * \brief Send the buffer asyncronously to a peer process.
     */
    void send(int peerRank, bool setNoAccess=true)
    {
#if HAVE_MPI
        MPI_Isend(data_,
                  mpiDataSize_,
                  mpiDataType_,
                  peerRank,
                  0, // tag
                  MPI_COMM_WORLD,
                  &mpiRequest_);
#endif
    }

    /*!
     * \brief Wait until the buffer was send to the peer completely.
     */
    void wait()
    {
#if HAVE_MPI
        MPI_Wait(&mpiRequest_, &mpiStatus_);
#endif // HAVE_MPI
    }

    /*!
     * \brief Receive the buffer syncronously from a peer rank
     */
    void receive(int peerRank)
    {
#if HAVE_MPI
        MPI_Recv(data_,
                 mpiDataSize_,
                 mpiDataType_,
                 peerRank,
                 0, // tag
                 MPI_COMM_WORLD,
                 &mpiStatus_);
        assert(! mpiStatus_.MPI_ERROR);
#endif // HAVE_MPI
    }

#if HAVE_MPI
    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    MPI_Request &request()
    { return mpiRequest_; }
    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    const MPI_Request &request() const
    { return mpiRequest_; }

    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    MPI_Status &status()
    { return mpiStatus_; }
    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    const MPI_Status &status() const
    { return mpiStatus_; }
#endif // HAVE_MPI


    /*!
     * \brief Returns the number of data objects in the buffer
     */
    size_t size() const
    { return dataSize_; }

    /*!
     * \brief Provide access to the buffer data.
     */
    DataType &operator[](size_t i)
    {
        assert(0 <= i && i < dataSize_);
        return data_[i];
    }

    /*!
     * \brief Provide access to the buffer data.
     */
    const DataType &operator[](size_t i) const
    {
        assert(0 <= i && i < dataSize_);
        return data_[i];
    }

private:
    void setMpiDataType_()
    {
#if HAVE_MPI
        // set the MPI data type
        if (std::is_same<DataType, char>::value)
            mpiDataType_ = MPI_CHAR;
        else if (std::is_same<DataType, unsigned char>::value)
            mpiDataType_ = MPI_UNSIGNED_CHAR;
        else if (std::is_same<DataType, short>::value)
            mpiDataType_ = MPI_SHORT;
        else if (std::is_same<DataType, unsigned short>::value)
            mpiDataType_ = MPI_UNSIGNED_SHORT;
        else if (std::is_same<DataType, int>::value)
            mpiDataType_ = MPI_INT;
        else if (std::is_same<DataType, unsigned>::value)
            mpiDataType_ = MPI_UNSIGNED;
        else if (std::is_same<DataType, long>::value)
            mpiDataType_ = MPI_LONG;
        else if (std::is_same<DataType, unsigned long>::value)
            mpiDataType_ = MPI_UNSIGNED_LONG;
        else if (std::is_same<DataType, long long>::value)
            mpiDataType_ = MPI_LONG_LONG;
        else if (std::is_same<DataType, unsigned long long>::value)
            mpiDataType_ = MPI_UNSIGNED_LONG_LONG;
        else if (std::is_same<DataType, float>::value)
            mpiDataType_ = MPI_FLOAT;
        else if (std::is_same<DataType, double>::value)
            mpiDataType_ = MPI_DOUBLE;
        else if (std::is_same<DataType, long double>::value)
            mpiDataType_ = MPI_LONG_DOUBLE;
        else {
            mpiDataType_ = MPI_BYTE;
        }
#endif // HAVE_MPI
    }

    void updateMpiDataSize_()
    {
#if HAVE_MPI
        mpiDataSize_ = dataSize_;
        if (mpiDataType_ == MPI_BYTE)
            mpiDataSize_ *= sizeof(DataType);
#endif // HAVE_MPI
    }

    DataType *data_;
    size_t dataSize_;
#if HAVE_MPI
    size_t mpiDataSize_;
    MPI_Datatype mpiDataType_;
    MPI_Request mpiRequest_;
    MPI_Status mpiStatus_;
#endif // HAVE_MPI
};

}

#endif
