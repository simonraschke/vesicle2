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



#include "common/definitions.hpp"
#include "particles/mobile.hpp"
#include "cell_state.hpp"
#include "enhance/container_class_base.hpp"
#include "observer_ptr.hpp"
#include "common/box.hpp"
#include "particles/interaction.hpp"
#include <deque>
#include <numeric>
#include <tbb/spin_rw_mutex.h>
#include <tbb/concurrent_vector.h>



namespace ves { struct Cell; };



struct ves::Cell
    : public enhance::ContainerBase<std::deque<nonstd::observer_ptr<Particle::Base>>>
{
public:
    using particle_ptr_t = element_t;
    using particle_t = element_t::element_type;
    using box3d = Particle::Base::box3d;
    using cartesian = Particle::Base::cartesian;

protected:
    tbb::spin_rw_mutex particles_access_mutex {};
    box3d bounding_box;
    
    ves::Box<PERIODIC::ON> box;
    AngularLennardJonesInteraction interaction;

    std::vector<std::reference_wrapper<Cell>> proximity {};
    std::vector<std::reference_wrapper<Cell>> region {};

public:
    Cell();

    inline bool operator==(const Cell& other) const { return std::addressof(*this) == std::addressof(other); }
    static bool areNeighbours(const Cell&, const Cell&);

    void setBoundaries(const cartesian&, const cartesian&);
    auto getBoundaries() const -> const decltype(bounding_box)&;

    void clear();
    void removeParticle(const particle_t&);
    bool try_add(particle_t*);
    auto getLeavers() -> decltype(data);
    REAL potential(const particle_t&) const;
    bool contains(const cartesian&) const;
    bool contains(const particle_t&);

    template<typename CONTAINER>
    void setupProximityAndRegion(CONTAINER&);

    inline const decltype(proximity)& getProximity() const { return proximity; }
    inline const decltype(region)& getRegion() const { return region; }

    CellState state {};

    template<CellState::STATE S>
    bool proximityAllInState() const;
    template<CellState::STATE S>
    bool proximityNoneInState() const;

    template<CellState::STATE S>
    bool regionAllInState() const;
    template<CellState::STATE S>
    bool regionNoneInState() const;
};



template<CellState::STATE S>
bool ves::Cell::proximityAllInState() const
{
    return std::all_of(std::begin(proximity),std::end(proximity), [](const ves::Cell& cell){ return cell.state == S; } );
}



template<CellState::STATE S>
bool ves::Cell::proximityNoneInState() const
{
    return std::none_of(std::begin(proximity),std::end(proximity), [](const ves::Cell& cell){ return cell.state == S; } );
}



template<CellState::STATE S>
bool ves::Cell::regionAllInState() const
{
    return std::all_of(std::begin(region),std::end(region), [](const ves::Cell& cell){ return cell.state == S; } );
}



template<CellState::STATE S>
bool ves::Cell::regionNoneInState() const
{
    return std::none_of(std::begin(region),std::end(region), [](const ves::Cell& cell){ return cell.state == S; } );
}



template<typename CONTAINER>
void ves::Cell::setupProximityAndRegion(CONTAINER& cells)
{
    // vesDEBUG(__PRETTY_FUNCTION__)
    for(Cell& cell : cells)
    {
        if(cell == *this)
        {
            region.push_back( std::ref(cell) );
        }
        else if( ves::Cell::areNeighbours(*this, cell) )
        {
            proximity.push_back( std::ref(cell) );
            region.push_back( std::ref(cell) );
        }
    }
}