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
 * \copydoc Opm::FvBaseDiscretization
 */
#ifndef EWOMS_FV_BASE_DISCRETIZATION_FEMADAPT_HH
#define EWOMS_FV_BASE_DISCRETIZATION_FEMADAPT_HH
#include "fvbasediscretization.hh"
#if HAVE_DUNE_FEM
namespace Opm
{

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief The base class for the finite volume discretization schemes.
 */
template <class TypeTag>
class FvBaseDiscretizationFemAdapt : public FvBaseDiscretization<TypeTag>
{
public:
    using Parent = FvBaseDiscretization<TypeTag>;
    // this constructor required to be explicitly specified because
    // we've defined a constructor above which deletes all implicitly
    // generated constructors in C++.
    using DiscreteFunctionSpace = GetPropType<TypeTag, Properties::DiscreteFunctionSpace>;

    // discrete function storing solution data
    using DiscreteFunction = Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace, PrimaryVariables>;

    // problem restriction and prolongation operator for adaptation
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ProblemRestrictProlongOperator = typename Problem ::RestrictProlongOperator;

    // discrete function restriction and prolongation operator for adaptation
    using DiscreteFunctionRestrictProlong = Dune::Fem::RestrictProlongDefault<DiscreteFunction>;
    using RestrictProlong
        = Dune::Fem::RestrictProlongTuple<DiscreteFunctionRestrictProlong, ProblemRestrictProlongOperator>;
    // adaptation classes
    using AdaptationManager = Dune::Fem::AdaptationManager<Grid, RestrictProlong>;
    FvBaseDiscretizationBaseFemAdapt(Simulator& simulator)
        : Parent(simulator)
    {
    }
    void adaptGrid()
    {
        // adapt the grid if enabled and if all dependencies are available
        // adaptation is only done if markForGridAdaptation returns true
        if (enableGridAdaptation_) {
            // check if problem allows for adaptation and cells were marked
            if (simulator_.problem().markForGridAdaptation()) {
                // adapt the grid and load balance if necessary
                adaptationManager().adapt();

                // if the grid has potentially changed, we need to re-create the
                // supporting data structures.
                elementMapper_.update(gridView_);
                vertexMapper_.update(gridView_);
                resetLinearizer();
                // this is a bit hacky because it supposes that Problem::finishInit()
                // works fine multiple times in a row.
                //
                // TODO: move this to Problem::gridChanged()
                finishInit();

                // notify the problem that the grid has changed
                //
                // TODO: come up with a mechanism to access the unadapted data structures
                // outside of the problem (i.e., grid, mappers, solutions)
                simulator_.problem().gridChanged();

                // notify the modules for visualization output
                auto outIt = outputModules_.begin();
                auto outEndIt = outputModules_.end();
                for (; outIt != outEndIt; ++outIt)
                    (*outIt)->allocBuffers();
            }
        }
    }



    AdaptationManager& adaptationManager()
    {
        if (!adaptationManager_) {
            // create adaptation objects here, because when doing so in constructor
            // problem is not yet intialized, aka seg fault
            restrictProlong_.reset(new RestrictProlong(DiscreteFunctionRestrictProlong(*(solution_[/*timeIdx=*/0])),
                                                       simulator_.problem().restrictProlongOperator()));
            adaptationManager_.reset(new AdaptationManager(simulator_.vanguard().grid(), *restrictProlong_));
        }
        return *adaptationManager_;
    }

private:
    std::unique_ptr<RestrictProlong> restrictProlong_;
    std::unique_ptr<AdaptationManager> adaptationManager_;
};
} // namespace Opm
#endif
#endif
