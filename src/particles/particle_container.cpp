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