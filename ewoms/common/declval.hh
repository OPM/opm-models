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
 * \brief Provides a std::declval equivalent.
 *
 * For compilers like GCC 4.4 which do not feature std::declval in their standard
 * library.
 */
#ifndef EWOMS_DECLVAL_HH
#define EWOMS_DECLVAL_HH

namespace Ewoms {
/*!
 * \brief Template function which returns an object of an arbitrary type.
 *
 * This is intended to be used in conjunction with the decltype keyword. It is required
 * because decltype requires an object as an argument and some objects do not provide a
 * default constructor. If you try to call the declval function at run time, you'll get a
 * compiler error.
 *
 * This function is equivalent to std::declval, but it is required for compilers like GCC
 * 4.4 which do not feature std::declval in their standard library.
 */
template <class T>
T& declval();
} // namespace Ewoms

#endif
