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
#include "parameters.hpp"
#include "particles/base.hpp"
#include "enhance/output_utility.hpp"
#include <filesystem>
#include <type_traits>
#include <boost/multi_array.hpp>

#ifdef __clang_major__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexceptions"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#include "h5xx/h5xx.hpp"
#pragma clang diagnostic pop
#elif  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wexceptions"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "h5xx/h5xx.hpp"
#pragma GCC diagnostic pop
#else
    #error no valid compiler
#endif



namespace ves 
{ 
    struct TrajectoryWriterH5;
    struct MonteCarloSystem;
};



struct ves::TrajectoryWriterH5
{
    using PATH = Parameters::PATH;
    using FSTREAM = std::ofstream;
    using array1d_t = boost::multi_array<std::uint32_t,1>;
    using array2d_t = boost::multi_array<REAL,2>;
    using cartesian = Particle::Base::cartesian;

    void setup();

    template<typename SYSTEM>
    void write(const SYSTEM&);
    

protected:
    template<typename SYSTEM>
    array1d_t getResidueTypes(const SYSTEM&) const;
    template<typename SYSTEM>
    array2d_t getPositions(const SYSTEM&) const;
    template<typename SYSTEM>
    array2d_t getOrientations(const SYSTEM&) const;
    template<typename SYSTEM>
    array2d_t getCoordinatesSphereBoundaries(const SYSTEM&) const;
    template<typename SYSTEM>
    array2d_t getCoordinatesBoxBoundaries(const SYSTEM&) const;
    template<typename SYSTEM>
    array2d_t getOrientationBoundaries(const SYSTEM&) const;

    PATH working_dir {std::filesystem::current_path()};
    PATH file_path {ves::Parameters::getInstance().getOption("output.path").as<PATH>()};
    h5xx::file h5file;

    const std::string filetype = {"h5"};
};



template<typename SYSTEM>
void ves::TrajectoryWriterH5::write(const SYSTEM& sys)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  simulation time: " << sys.getTime());
    
    assert(h5file.valid());
    assert(file_path);
    assert(target_range);

    const std::string group_name("snapshot"+std::to_string(static_cast<std::size_t>(sys.getTime())));
    h5xx::group group(h5file, group_name);

    {
        h5xx::write_attribute(group, "general.acceptance", Parameters::getInstance().getOption("general.acceptance").as<std::string>());
        h5xx::write_attribute(group, "general.ensemble", Parameters::getInstance().getOption("general.ensemble").as<std::string>());
        h5xx::write_attribute(group, "general.simulation_mode", Parameters::getInstance().getOption("general.simulation_mode").as<std::string>());
        h5xx::write_attribute(group, "system.frame_guides_grid_edge", Parameters::getInstance().getOption("system.frame_guides_grid_edge").as<std::size_t>());
        h5xx::write_attribute(group, "system.guiding_elements_each", Parameters::getInstance().getOption("system.guiding_elements_each").as<std::size_t>());
        h5xx::write_attribute(group, "system.plane_edge", Parameters::getInstance().getOption("system.plane_edge").as<REAL>());
        h5xx::write_attribute(group, "system.osmotic_density_inside", Parameters::getInstance().getOption("system.osmotic_density_inside").as<REAL>());
        h5xx::write_attribute(group, "system.osmotic_density_outside", Parameters::getInstance().getOption("system.osmotic_density_outside").as<REAL>());
        h5xx::write_attribute(group, "system.num_mobile_particles", sys.getParticles().get().template numType<ves::Particle::TYPE::MOBILE>());
        h5xx::write_attribute(group, "system.num_frame_particles", sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>());
        h5xx::write_attribute(group, "system.num_osmotic_particles", sys.getParticles().get().template numType<ves::Particle::TYPE::OSMOTIC>());
        h5xx::write_attribute(group, "system.num_all_particles", sys.getParticles().get().data.size());
        h5xx::write_attribute(group, "system.box.x", Parameters::getInstance().getOption("system.box.x").as<REAL>());
        h5xx::write_attribute(group, "system.box.y", Parameters::getInstance().getOption("system.box.y").as<REAL>());
        h5xx::write_attribute(group, "system.box.z", Parameters::getInstance().getOption("system.box.z").as<REAL>());
        h5xx::write_attribute(group, "system.temperature", Parameters::getInstance().getOption("system.temperature").as<REAL>());
        h5xx::write_attribute(group, "system.actual_time", sys.getTime());
        h5xx::write_attribute(group, "system.time_max", Parameters::getInstance().getOption("system.time_max").as<std::size_t>());
        h5xx::write_attribute(group, "system.ljepsilon", Parameters::getInstance().getOption("system.ljepsilon").as<REAL>());
        h5xx::write_attribute(group, "system.ljsigma", Parameters::getInstance().getOption("system.ljsigma").as<REAL>());
        h5xx::write_attribute(group, "system.kappa", Parameters::getInstance().getOption("system.kappa").as<REAL>());
        h5xx::write_attribute(group, "system.gamma", Parameters::getInstance().getOption("system.gamma").as<REAL>());
        h5xx::write_attribute(group, "system.sw_position_min", Parameters::getInstance().getOption("system.sw_position_min").as<REAL>());
        h5xx::write_attribute(group, "system.sw_position_max", Parameters::getInstance().getOption("system.sw_position_max").as<REAL>());
        h5xx::write_attribute(group, "system.sw_position_actual", sys.getSWPosition().get()());
        h5xx::write_attribute(group, "system.acceptance_position_target", Parameters::getInstance().getOption("system.acceptance_position_target").as<REAL>());
        h5xx::write_attribute(group, "system.sw_orientation_min", Parameters::getInstance().getOption("system.sw_orientation_min").as<REAL>());
        h5xx::write_attribute(group, "system.sw_orientation_max", Parameters::getInstance().getOption("system.sw_orientation_max").as<REAL>());
        h5xx::write_attribute(group, "system.sw_orientation_actual", sys.getSWOrientation().get()());
        h5xx::write_attribute(group, "system.acceptance_orientation_target", Parameters::getInstance().getOption("system.acceptance_orientation_target").as<REAL>());
        h5xx::write_attribute(group, "system.cell_min_edge", Parameters::getInstance().getOption("system.cell_min_edge").as<REAL>());
        h5xx::write_attribute(group, "system.max_cells_dim", Parameters::getInstance().getOption("system.max_cells_dim").as<std::size_t>());
        h5xx::write_attribute(group, "output.path", ves::Parameters::getInstance().getOption("output.path").as<PATH>().string());
        h5xx::write_attribute(group, "output.skip", Parameters::getInstance().getOption("output.skip").as<std::size_t>());
        h5xx::write_attribute(group, "input.path", Parameters::getInstance().getOption("input.path").as<PATH>().string());
    }
    
    {
        const std::string dataset_name("position");
        auto dataset = getPositions(sys);
        h5xx::create_dataset(group, dataset_name, dataset);
        h5xx::write_dataset(group, dataset_name, dataset);
    }
    
    {
        const std::string dataset_name("orientation");
        auto dataset = getOrientations(sys);
        h5xx::create_dataset(group, dataset_name, dataset);
        h5xx::write_dataset(group, dataset_name, dataset);
    }
    
    {
        const std::string dataset_name("type");
        auto dataset = getResidueTypes(sys);
        h5xx::create_dataset(group, dataset_name, dataset);
        h5xx::write_dataset(group, dataset_name, dataset);
    }


    if(sys.getTime() == 0) 
    {
        h5xx::group rootgroup(h5file, "/");
        try
        {
            const std::string dataset_name("position_sphere_bounds");
            auto dataset = getCoordinatesSphereBoundaries(sys);
            h5xx::create_dataset(rootgroup, dataset_name, dataset);
            h5xx::write_dataset(rootgroup, dataset_name, dataset);
        }
        catch(...)
        {
            try
            {
                const std::string dataset_name("position_box_bounds");
                auto dataset = getCoordinatesBoxBoundaries(sys);
                h5xx::create_dataset(rootgroup, dataset_name, dataset);
                h5xx::write_dataset(rootgroup, dataset_name, dataset);
            }
            catch(...)
            {
                ;
            }
        }
        
        try
        {
            const std::string dataset_name("orientation_bounds");
            auto dataset = getOrientationBoundaries(sys);
            h5xx::create_dataset(rootgroup, dataset_name, dataset);
            h5xx::write_dataset(rootgroup, dataset_name, dataset);
        }
        catch(...)
        {
            ;
        }
    }
    
    // {
    //     const std::string dataset_name("IDs");
    //     auto IDs = getParticleIDs(sys);
    //     h5xx::create_dataset(group, dataset_name, IDs);
    //     h5xx::write_dataset(group, dataset_name, IDs);
    // }
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getPositions(const SYSTEM& sys) const
{
    typedef array2d_t::index index;

    array2d_t array(boost::extents[sys.getParticles().get().data.size()][3]);

    for(index residue = 0; residue < static_cast<index>(sys.getParticles().get().data.size()); ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const cartesian& coords = sys.getBox().get().scaleDown( target->getCoordinates());

        array[residue][static_cast<index>(0)] = coords(0);
        array[residue][static_cast<index>(1)] = coords(1);
        array[residue][static_cast<index>(2)] = coords(2);
    }
    
    return array;
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getOrientations(const SYSTEM& sys) const
{
    typedef array2d_t::index index;

    array2d_t array(boost::extents[sys.getParticles().get().data.size()][3]);

    for(index residue = 0; residue < static_cast<index>(sys.getParticles().get().data.size()); ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const cartesian orientation = target->getType() == ves::Particle::TYPE::OSMOTIC ? cartesian(0,0,0) : target->getOrientation().normalized();

        array[residue][static_cast<index>(0)] = orientation(0);
        array[residue][static_cast<index>(1)] = orientation(1);
        array[residue][static_cast<index>(2)] = orientation(2);
    }
    
    return array;
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getCoordinatesSphereBoundaries(const SYSTEM& sys) const
{
    typedef array2d_t::index index;

    const index num = sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>();
    if(num == 0)
        throw std::runtime_error("no guiding elements");

    array2d_t array(boost::extents[sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>()][5]);

    for(index residue = 0; residue < num; ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const cartesian origin = *(target->coordinates_bounding.origin);
        const auto fromto = target->coordinates_bounding.getSphereBounds();

        array[residue][static_cast<index>(0)] = origin(0);
        array[residue][static_cast<index>(1)] = origin(1);
        array[residue][static_cast<index>(2)] = origin(2);
        array[residue][static_cast<index>(3)] = fromto.first;
        array[residue][static_cast<index>(4)] = fromto.second;
    }
    
    return array;
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getCoordinatesBoxBoundaries(const SYSTEM& sys) const
{
    typedef array2d_t::index index;

    const index num = sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>();
    if(num == 0)
        throw std::runtime_error("no guiding elements");

    array2d_t array(boost::extents[sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>()][6]);

    for(index residue = 0; residue < num; ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const auto box = target->coordinates_bounding.getBoundingBox();
        const auto mincorner = box.corner(ves::Particle::Base::box3d::BottomLeftFloor);
        const auto maxcorner = box.corner(ves::Particle::Base::box3d::TopRightCeil);

        array[residue][static_cast<index>(0)] = mincorner(0);
        array[residue][static_cast<index>(1)] = mincorner(1);
        array[residue][static_cast<index>(2)] = mincorner(2);
        array[residue][static_cast<index>(3)] = maxcorner(0);
        array[residue][static_cast<index>(4)] = maxcorner(1);
        array[residue][static_cast<index>(5)] = maxcorner(2);
    }
    
    return array;
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getOrientationBoundaries(const SYSTEM& sys) const
{
    typedef array2d_t::index index;

    const index num = sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>();
    if(num == 0)
        throw std::runtime_error("no guiding elements");

    array2d_t array(boost::extents[sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>()][4]);

    for(index residue = 0; residue < num; ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const cartesian origin = *(target->orientation_bounding.origin);
        const REAL value = target->orientation_bounding.value;

        array[residue][static_cast<index>(0)] = origin(0);
        array[residue][static_cast<index>(1)] = origin(1);
        array[residue][static_cast<index>(2)] = origin(2);
        array[residue][static_cast<index>(3)] = value;
    }
    
    return array;
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array1d_t ves::TrajectoryWriterH5::getResidueTypes(const SYSTEM& sys) const
{
    typedef array1d_t::index index;

    array1d_t array(boost::extents[sys.getParticles().get().data.size()]);

    for(index residue = 0; residue < static_cast<index>(sys.getParticles().get().data.size()); ++residue)
    {
        array[residue] = static_cast<std::underlying_type<ves::Particle::TYPE>::type>(sys.getParticles().get().data[residue]->getType());
    }
    
    return array;
}



// template<typename SYSTEM>
// ves::TrajectoryWriterH5::array1d_t ves::TrajectoryWriterH5::getParticleIDs(const SYSTEM& sys) const
// {
//     typedef array1d_t::index index;

//     array1d_t array(boost::extents[sys.getParticles().get().data.size()]);

//     for(index residue = 0; residue < static_cast<index>(sys.getParticles().get().data.size()); ++residue)
//     {
//         array[residue] = sys.getParticles().get().data[residue]->ID;
//     }
    
//     return array;
// }