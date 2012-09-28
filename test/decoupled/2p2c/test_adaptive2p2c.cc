/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#include "config.h"

#include "test_adaptive2p2cproblem.hh"
#include <dumux/common/start.hh>

// The main function using the standard start procedure
int main(int argc, char** argv)
{
    typedef TTAG(Adaptive2p2c) TypeTag;
    return Dumux::start<TypeTag>(argc, argv);
}
