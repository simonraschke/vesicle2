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

#include "particle_container.hpp"



void ves::ParticleContainer::setup()
{
    vesDEBUG(__PRETTY_FUNCTION__);

    if(std::filesystem::exists(Parameters::getInstance().getOption("input.path").as<Parameters::PATH>()))
    {
        GLOBAL::getInstance().startmode.store(GLOBAL::STARTMODE::RESTART);
        setupFromH5();
    }
    else
    {
        GLOBAL::getInstance().startmode.store(GLOBAL::STARTMODE::NEW);
        setupFromNew();
    }
}



void ves::ParticleContainer::setupFromNew()
{
    switch (GLOBAL::getInstance().simulationmode.load())
    {
        case GLOBAL::SIMULATIONMODE::SA:
            try
            {
                REAL x = Parameters::getInstance().getOption("system.box.x").as<REAL>();
                REAL y = Parameters::getInstance().getOption("system.box.y").as<REAL>();
                REAL z = Parameters::getInstance().getOption("system.box.z").as<REAL>();
                box.setLengthX(x);
                box.setLengthY(y);
                box.setLengthZ(z);
            }
            catch(std::logic_error e)
            {
                vesLOG("unable to get box from x|y|z, will calculate cubic box from density and mobile particles");
                try
                {
                    REAL density = Parameters::getInstance().getOption("system.density").as<REAL>();
                    std::size_t mobile = Parameters::getInstance().getOption("system.mobile").as<std::size_t>();
                    REAL volume = static_cast<REAL>(mobile)/density;
                    REAL edge = std::cbrt(volume);
                    box.setLengthX(edge);
                    box.setLengthY(edge);
                    box.setLengthZ(edge);
                }
                catch(std::logic_error e)
                {
                    vesCRITICAL("unable to get box from density|mobile, abort");
                }
                vesLOG("box(x|y|z) = (" << box.getLengthX() << "|" << box.getLengthY() << "|" << box.getLengthZ() << ")");
                
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.x", boost::program_options::variable_value(box.getLengthX(), false)));
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.y", boost::program_options::variable_value(box.getLengthY(), false)));
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.z", boost::program_options::variable_value(box.getLengthZ(), false)));
            }
            for(std::size_t i = 0; i < Parameters::getInstance().getOption("system.mobile").as<std::size_t>(); ++i)
            {
                const auto minimum_offset = Parameters::getInstance().getOption("system.ljsigma").as<REAL>();
                addParticle<Particle::Mobile>();
                do
                {
                    // vesLOG("tries: " << ++k);
                    data.back()->try_setCoordinates(box.randomPointInside());
                }
                while(placement_conflict(*(data.back()), minimum_offset));
                data.back()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized());
                // vesLOG("placed particle: " << data.back()->getCoordinates().format(ROWFORMAT));
            }
            

            break;

        case GLOBAL::SIMULATIONMODE::FGA:
            break;

        case GLOBAL::SIMULATIONMODE::OSMOTIC:
            break;
    
        default:
            break;
    }
}



std::vector<std::string> ves::ParticleContainer::getGroupNames()
{
    vesLOG(__PRETTY_FUNCTION__);
    auto file_info = [](hid_t loc_id, const char *name, const H5L_info_t * __attribute__((unused)) linfo, void *opdata) -> herr_t
    {
        hid_t group;
        auto group_names = reinterpret_cast< std::vector<std::string>* >(opdata);
        group = H5Gopen2(loc_id, name, H5P_DEFAULT);
        //do stuff with group object, if needed
        group_names->push_back(name);
        // cout << "Name : " << name << endl;
        H5Gclose(group);
        return 0;
    };
    
    if(!h5file.valid())
        h5file.open(Parameters::getInstance().getOption("input.path").as<Parameters::PATH>().string(), h5xx::file::in);
    h5xx::group root_group(h5file, "/");

    std::vector<std::string> group_names;
    herr_t  __attribute__((unused)) idx = H5Literate(root_group.hid(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, &group_names);

    const std::string numbers = "0123456789";
    std::sort(std::begin(group_names), std::end(group_names), [&](const std::string& a, const std::string& b) 
    { 
        std::string cmp_a { a.begin() + a.find_first_of(numbers), a.end() };
        std::string cmp_b { b.begin() + b.find_first_of(numbers), b.end() };
        return std::stoi(cmp_a) < std::stoi(cmp_b);
    });
    // group_names.sort();
    assert(std::is_sorted(group_names));
    return group_names;
}



std::size_t ves::ParticleContainer::getTime() const
{
    return time_from_input;
}



void ves::ParticleContainer::setupFromH5()
{
    vesDEBUG(__PRETTY_FUNCTION__);
    const auto filepath = Parameters::getInstance().getOption("input.path").as<Parameters::PATH>();
    vesLOG("distributing particles from " << Parameters::getInstance().getOption("input.path").as<Parameters::PATH>());

    if(!h5file.valid())
    {
            vesDEBUG("will now open " << filepath.string());
            h5file.open(filepath.string(), h5xx::file::in);
    }
    std::string latest_group_name = getGroupNames().back();
    vesLOG("get group" << latest_group_name);
    h5xx::group group(h5file, latest_group_name);

    REAL x = h5xx::read_attribute<REAL>(group, "system.box.x");
    REAL y = h5xx::read_attribute<REAL>(group, "system.box.y");
    REAL z = h5xx::read_attribute<REAL>(group, "system.box.z");
    box.setLengthX(x);
    box.setLengthY(y);
    box.setLengthZ(z);
    Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.x", boost::program_options::variable_value(box.getLengthX(), false)));
    Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.y", boost::program_options::variable_value(box.getLengthY(), false)));
    Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.z", boost::program_options::variable_value(box.getLengthZ(), false)));
    Parameters::getInstance().mutableAccess().insert(std::make_pair("system.sw_position_actual", boost::program_options::variable_value(h5xx::read_attribute<REAL>(group, "system.sw_position_actual"), false)));
    Parameters::getInstance().mutableAccess().insert(std::make_pair("system.sw_orientation_actual", boost::program_options::variable_value(h5xx::read_attribute<REAL>(group, "system.sw_orientation_actual"), false)));
    
    time_from_input = h5xx::read_attribute<decltype(time_from_input)>(group, "system.actual_time")+1;
    
    array1d_t types;
    array1d_t IDs;
    array2d_t positions;
    array2d_t orientations;

    {
        vesLOG("load type");
        h5xx::dataset dataset_type(group, "type");
        h5xx::read_dataset(dataset_type, types);
    }

    {
        vesLOG("load pos");
        h5xx::dataset dataset_positions(group, "position");
        h5xx::read_dataset(dataset_positions, positions);
    }

    {
        vesLOG("load orien");
        h5xx::dataset dataset_orientations(group, "orientation");
        h5xx::read_dataset(dataset_orientations, orientations);
    }

    // generate particles
    for(std::size_t i = 0; i < types.size() ; ++i)
    {
        switch (types[i])
        {
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::MOBILE):
                addParticle<ves::Particle::Mobile>();
                break;
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::FRAME):
                addParticle<ves::Particle::Frame>();
                break;
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::OSMOTIC):
                addParticle<ves::Particle::Osmotic>();
                break;
        
            default:
                throw std::logic_error("undefined particle type in TrajectoryDistributorH5 in row "+std::to_string(i));
                break;
        }
    }

    // set 'em up
    for(std::size_t i = 0; i < types.size() ; ++i)
    {
        vesDEBUG(i << cartesian(positions[i][0], positions[i][1], positions[i][2]).format(ROWFORMAT) << cartesian(orientations[i][0], orientations[i][1], orientations[i][2]).format(ROWFORMAT));
        
        data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);

        switch (types[i])
        {
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::MOBILE):
                break;
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::FRAME):
                break;
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::OSMOTIC):
                break;
        
            default:
                throw std::logic_error("undefined particle type in TrajectoryDistributorH5 in row "+std::to_string(i));
                break;
        }

        data.at(i)->getOrientation() = cartesian(orientations[i][0], orientations[i][1], orientations[i][2]);
    }

    h5file.close();
}



bool ves::ParticleContainer::placement_conflict(const particle_t& p, REAL minimum_distance) const
{
    const REAL minimum_distance_squared = minimum_distance * minimum_distance;
    return std::find_if_not(std::begin(data), std::end(data), [&](const auto& p_ptr)
    {
        return (p == *p_ptr) || box.squared_distance(p, *p_ptr) > minimum_distance_squared;
    }) != std::end(data);
}



void ves::ParticleContainer::removeParticle(const particle_t& p)
{
    tbb::mutex::scoped_lock lock(mutex);
    data.erase(std::remove_if(begin(), end(), [&](const auto& comp){ return p == *comp; }), end());
}



void ves::ParticleContainer::removeParticle(const particle_ptr_t& p)
{
    tbb::mutex::scoped_lock lock(mutex);
    data.erase(std::remove(begin(), end(), p), end());
}