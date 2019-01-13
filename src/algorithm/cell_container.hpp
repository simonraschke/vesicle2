/*  
*   Copyright 2019 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/



#pragma once



#include "cell.hpp"
// #include <tbb/concurrent_vector.h>
#include <tbb/parallel_for_each.h>



namespace ves { struct CellContainer; };



struct ves::CellContainer
    : public enhance::ContainerBase<std::deque<Cell>>
{
    using cartesian = Cell::cartesian;
    using particle_t = Cell::particle_t;

    void setup();

    template<typename CONTAINER>
    void deployParticles(const CONTAINER&);

    void deployParticle(const particle_t&);

    template<CellState::STATE S>
    bool allInState() const;

    template<CellState::STATE S>
    bool noneInState() const;

protected:
    ves::Box<PERIODIC::ON> box;
};



template<CellState::STATE S>
bool ves::CellContainer::allInState() const
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    return std::all_of(begin(), end(), [](const Cell& cell){ return cell.state == S; } );
}



template<CellState::STATE S>
bool ves::CellContainer::noneInState() const
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    return std::none_of(begin(), end(), [](const Cell& cell){ return cell.state == S; } );
}



template<typename CONTAINER>
void ves::CellContainer::deployParticles(const CONTAINER& particles)
{
    vesDEBUG(__PRETTY_FUNCTION__);
    tbb::parallel_for_each(std::begin(particles), std::end(particles), [&](const auto& particle)
    {
        bool done = false;
        std::for_each(begin(), end(), [&](auto & cell)
        {
            if(done) return;
            done = cell.try_add(particle.get());
        });
        assert(done);
    });
}