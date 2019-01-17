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
// #include "enhance/random.hpp"
#include "enhance/iterator_utility.hpp"
#include "enhance/parallel.hpp"
#include <tbb/parallel_for_each.h>
#include <vector>



namespace ves { struct CellContainer; };



struct ves::CellContainer
    : public enhance::ContainerBase<std::deque<Cell>>
{
    using cartesian = Cell::cartesian;
    using particle_t = Cell::particle_t;

    void setup();
    void reorder();
    void preparation();
    std::size_t membersContained() const;

    template<typename FUNCTOR>
    void cellBasedApplyFunctor(FUNCTOR&& func);

    template<typename CONTAINER>
    void deployParticles(const CONTAINER&);

    decltype(data)::iterator getCellOfCartesian(const cartesian&);

    decltype(data)::iterator deployParticle(particle_t&);
    void removeParticle(particle_t&);

    template<CellState::STATE S>
    bool allInState() const;

    template<CellState::STATE S>
    bool noneInState() const;


protected:
    std::vector<std::size_t> iteration_IDs;
    ves::Box<PERIODIC::ON> box;

    thread_local static std::mt19937_64 pseudo_engine;
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



template<typename FUNCTOR>
void ves::CellContainer::cellBasedApplyFunctor(FUNCTOR&& func)
{
    std::shuffle(std::begin(iteration_IDs), std::end(iteration_IDs), pseudo_engine);
    
    enhance::scoped_root_dummy ROOT;

    while( ! (allInState<CellState::STATE::FINISHED>()) )
    {
        // enhance::for_each_from(
        //     begin(), 
        //     end(),
        //     begin() + enhance::random<std::size_t>(0, data.size()-1), 
        //     [&](Cell& cell)
        for(auto ID : iteration_IDs)
        {
            ves::Cell& cell = data[ID];
            if( 
                cell.regionNoneInState<CellState::STATE::BLOCKED>() && 
                cell.state == CellState::STATE::IDLE
            )
            {
                cell.state = CellState::STATE::BLOCKED;
                
                assert( cell.state == CellState::STATE::BLOCKED );
                assert( cell.proximityNoneInState<CellState::STATE::BLOCKED>() );
                
                ROOT.enqueue_child( [&]
                {
                    assert( cell.state == CellState::STATE::BLOCKED );
                    func( cell ); 
                    cell.state = CellState::STATE::FINISHED;
                    assert( cell.state == CellState::STATE::FINISHED );
                } );
            }
        // });
        }
    }
    
    vesDEBUG("now check if all are FINISHED");
    assert( allInState<CellState::STATE::FINISHED>() );
}