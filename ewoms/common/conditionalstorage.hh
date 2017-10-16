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
 *
 * \copydoc Ewoms::ConditionalStorage
 */
#ifndef EWOMS_CONDITIONAL_STORAGE_HH
#define EWOMS_CONDITIONAL_STORAGE_HH

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <utility>

namespace Ewoms {
/*!
 * \ingroup Common
 *
 * \brief A simple class which only stores a given member attribute if a boolean
 *        condition is true
 *
 * If the condition is false, nothing is stored and an exception is thrown when trying to
 * access the object.
 */
template <bool cond, class T>
class ConditionalStorage
{
public:
    typedef T type;
    static const bool value = cond;

    ConditionalStorage()
    {}

    ConditionalStorage(const T& t)
        : data_(t)
    {};

    ConditionalStorage(T&& t)
        : data_(std::move(t))
    {};

    const T& operator*() const
    { return data_; }
    T& operator*()
    { return data_; }

    const T* operator->() const
    { return &data_; }
    T* operator->()
    { return &data_; }

private:
    T data_;
};

template <class T>
class ConditionalStorage<false, T>
{
public:
    typedef T type;
    static const bool value = false;

    ConditionalStorage()
    {}

    ConditionalStorage(const T&)
    {};

    const T& operator*() const
    { OPM_THROW(std::logic_error, "data member deactivated"); }
    T& operator*()
    { OPM_THROW(std::logic_error, "data member deactivated"); }

    const T* operator->() const
    { OPM_THROW(std::logic_error, "data member deactivated"); }
    T* operator->()
    { OPM_THROW(std::logic_error, "data member deactivated"); }
};

} // namespace Ewoms

#endif
