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



#include "cell_container.hpp"



thread_local std::mt19937_64 ves::CellContainer::pseudo_engine(std::random_device().operator()());



void ves::CellContainer::setup()
{
    box.setLengthX(Parameters::getInstance().getOption("system.box.x").as<REAL>());
    box.setLengthY(Parameters::getInstance().getOption("system.box.y").as<REAL>());
    box.setLengthZ(Parameters::getInstance().getOption("system.box.z").as<REAL>());

    const REAL min_edge = Parameters::getInstance().getOption("system.cell_min_edge").as<REAL>();
    const std::size_t max_cells_dim = Parameters::getInstance().getOption("system.max_cells_dim").as<std::size_t>();

    vesDEBUG("box x y z: " << box.getLengthX()  << " " << box.getLengthY()  << " " << box.getLengthZ())
    vesDEBUG("minimum edge for a cell is " << min_edge)
    vesDEBUG("maximum number of cells per dimension is " << max_cells_dim)

    // calculate the amount of zells per dimension
    std::size_t x_helper = std::trunc( box.getLengthX() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( box.getLengthX() / (min_edge) )  ;
    std::size_t y_helper = std::trunc( box.getLengthY() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( box.getLengthY() / (min_edge) )  ;
    std::size_t z_helper = std::trunc( box.getLengthZ() / (min_edge) ) > max_cells_dim ? max_cells_dim : std::trunc( box.getLengthZ() / (min_edge) )  ;

    // if dimension is smaller than min_edge it will be 0, therefor set to 1
    const std::size_t cells_x = x_helper == 0 ? 1 : x_helper;
    const std::size_t cells_y = y_helper == 0 ? 1 : y_helper;
    const std::size_t cells_z = z_helper == 0 ? 1 : z_helper;

    // edge of cell per dimension
    const float x_edge = box.getLengthX()/cells_x;
    const float y_edge = box.getLengthY()/cells_y;
    const float z_edge = box.getLengthZ()/cells_z;
    
    vesLOG("building cell environment " << cells_x << " " << cells_y << " " << cells_z)
    for ( std::size_t x = 0; x < cells_x; ++x ) 
    for ( std::size_t y = 0; y < cells_y; ++y )
    for ( std::size_t z = 0; z < cells_z; ++z )
    {
        data.emplace_back();
        const cartesian from = cartesian( x_edge*x,     y_edge*y,     z_edge*z     ) - cartesian(0.01, 0.01, 0.01);
        const cartesian to   = cartesian( x_edge*(x+1), y_edge*(y+1), z_edge*(z+1) ) + cartesian(0.01, 0.01, 0.01);
        data.back().setBoundaries(from,to);
    }

    assert( data.size() == cells_x*cells_y*cells_z );

    vesLOG("setup cells proximity and region")
    tbb::parallel_for_each(begin(), end(), [&](Cell& cell)
    {
        cell.setupProximityAndRegion(data);
        cell.state = CellState::STATE::IDLE;
    });

    #ifndef NDEBUG
    std::for_each(begin(), end(), [&](Cell& cell)
    {
        assert(cell.getProximity().size() == 26);
        assert(cell.getRegion().size() == 27);
        assert(cell.state == CellState::STATE::IDLE);
    });
    #endif

    iteration_IDs.resize(data.size());
    std::iota(std::begin(iteration_IDs), std::end(iteration_IDs), 0);
}



void ves::CellContainer::reorder()
{
    tbb::parallel_for_each(begin(), end(), [&](Cell& cell)
    {
        auto leavers = cell.getLeavers();
        vesDEBUG("cell got "<< leavers.size() << " leaver");
        
        for(auto& leaver : leavers)
        {
            bool was_added = false;
            for(Cell& proximity_cell : cell.getProximity())
            {
                was_added = proximity_cell.try_add(leaver.get());
                if(was_added) 
                {
                    vesDEBUG("leaver " << leaver->getCoordinates().format(ROWFORMAT) << " DID fit into cell with bounds from" << proximity_cell.getBoundaries().corner(Particle::Base::box3d::BottomLeftFloor).format(ROWFORMAT) << " to " <<  proximity_cell.getBoundaries().corner(Particle::Base::box3d::TopRightCeil).format(ROWFORMAT));
                    break;
                }
                vesDEBUG("leaver " << leaver->getCoordinates().format(ROWFORMAT) << " DID NOT fit into cell with bounds from " << proximity_cell.getBoundaries().corner(Particle::Base::box3d::BottomLeftFloor).format(ROWFORMAT) << " to " <<  proximity_cell.getBoundaries().corner(Particle::Base::box3d::TopRightCeil).format(ROWFORMAT));
            }
            assert(was_added);
            cell.removeParticle(*leaver);
            assert(!cell.contains(*leaver));
        }
    });
}



void ves::CellContainer::preparation()
{
    tbb::parallel_for_each(begin(), end(), [](Cell& cell)
    {
        cell.state = CellState::STATE::IDLE;
    });
}



std::size_t ves::CellContainer::membersContained() const
{
    return std::accumulate(begin(), end(), std::size_t(0), [](auto i, const Cell& cell){ return i + cell.data.size(); });
}