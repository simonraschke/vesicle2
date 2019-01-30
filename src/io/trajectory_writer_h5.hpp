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
#include <type_traits>
#include <boost/multi_array.hpp>

#ifdef __clang_major__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexceptions"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wdelete-non-virtual-dtor"
#pragma clang diagnostic ignored "-Weffc++"
#include "h5xx/h5xx.hpp"
#pragma clang diagnostic pop
#elif  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wexceptions"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
#pragma GCC diagnostic ignored "-Weffc++"
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
    using FSTREAM = fs::ofstream;
    using array1d_t = boost::multi_array<std::uint32_t,1>;
    using array2d_t = boost::multi_array<REAL,2>;
    using cartesian = Particle::Base::cartesian;

    template<typename SYSTEM>
    void setup(const SYSTEM&);

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

    PATH working_dir {fs::current_path()};
    PATH file_path {ves::Parameters::getInstance().getOption("output.path").as<PATH>()};
    h5xx::file h5file;

    const std::string filetype = {"h5"};
};






template<typename SYSTEM>
void ves::TrajectoryWriterH5::setup([[maybe_unused]] const SYSTEM& sys)
{
    file_path = PATH( *(enhance::splitAtDelimiter(file_path.string(), ".").rbegin()+1)+"."+filetype );
    vesDEBUG(__PRETTY_FUNCTION__<< "  " << file_path);
    try
    {
        h5file.close();
    }
    catch(...)
    {
        vesWARNING("unable to close h5file");
    }


    if(GLOBAL::getInstance().startmode == GLOBAL::STARTMODE::NEW)
    {

        vesLOG("trunc open " << file_path.string());
        h5file.open(file_path.string(), h5xx::file::trunc);

        // if(GLOBAL::getInstance().ensemble == GLOBAL::ENSEMBLE::uVT)
        // {
        //     vesLOG("writing constants")
        //     const REAL temperature = Parameters::getInstance().getOption("system.temperature").as<REAL>();
        //     const REAL thermal_wavelength_cubic = std::pow(std::sqrt(REAL(1) / (PI*2*temperature)), 3);
        //     const REAL mu = -temperature*std::log(sys.getBox().get().getVolume()/(thermal_wavelength_cubic*sys.getParticles().get().data.size()));

        //     // h5xx::group rootgroup(h5file, "/");
        //     // h5xx::write_attribute(rootgroup, "constant.mu", mu);
            
        //     // Parameters::getInstance().mutableAccess().insert(std::make_pair("constant.mu", boost::program_options::variable_value(mu, false)));
        // }

    }
    else
    {
        vesLOG("app open " << file_path.string());
        h5file.open(file_path.string(), h5xx::file::out);
    }

}




template<typename SYSTEM>
void ves::TrajectoryWriterH5::write(const SYSTEM& sys)
{
    vesDEBUG(__PRETTY_FUNCTION__<< "  simulation time: " << sys.getTime());
    
    assert(h5file.valid());

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
        h5xx::write_attribute(group, "system.mu", Parameters::getInstance().getOption("system.mu").as<REAL>());
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
        catch(const std::runtime_error& e)
        {
            vesLOG(e.what());
            vesLOG("no position_sphere_bounds found. will try to get position_box_bounds");
            try
            {
                const std::string dataset_name("position_box_bounds");
                auto dataset = getCoordinatesBoxBoundaries(sys);
                h5xx::create_dataset(rootgroup, dataset_name, dataset);
                h5xx::write_dataset(rootgroup, dataset_name, dataset);
            }
            catch(const std::runtime_error& e)
            {
                vesLOG(e.what());
                vesLOG("no position_box_bounds found");
                if(sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>() > 0)
                    vesCRITICAL(sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>()<< " GUIDING ELEMENTS HAVE NO BOUNDARIES");
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
            vesLOG("no orientation_bounds found");
        }
    }
}



template<typename SYSTEM>
ves::TrajectoryWriterH5::array2d_t ves::TrajectoryWriterH5::getPositions(const SYSTEM& sys) const
{
    typedef array2d_t::index index;
    const index particles_with_positions = sys.getParticles().get().data.size();

    array2d_t array(boost::extents[particles_with_positions][3]);

    for(index residue = 0; residue < particles_with_positions; ++residue)
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
    const index particles_with_orientations = sys.getParticles().get().template numType<ves::Particle::TYPE::MOBILE>() + sys.getParticles().get().template numType<ves::Particle::TYPE::FRAME>();

    array2d_t array(boost::extents[particles_with_orientations][3]);

    for(index residue = 0; residue < particles_with_orientations; ++residue)
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

    array2d_t array(boost::extents[num][5]);

    for(index residue = 0; residue < num; ++residue)
    {
        const auto& target = sys.getParticles().get().data[residue];
        const cartesian origin = *(target->coordinates_bounding.origin);
        
        // vesLOG(target->getCoordinates().format(ROWFORMAT));
        if(! target->coordinates_bounding.isSphereBound())
            throw std::runtime_error("particle has no sphere bounds");

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

        if(! target->coordinates_bounding.isBoxBound())
            throw std::runtime_error("particle has no box bounds");

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