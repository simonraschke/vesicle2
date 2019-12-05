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

    if(fs::exists(Parameters::getInstance().getOption("input.path").as<Parameters::PATH>()))
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



bool ves::ParticleContainer::placement_conflict(const particle_t& p, REAL minimum_distance) const
{
    const REAL minimum_distance_squared = minimum_distance * minimum_distance;
    return std::find_if_not(std::begin(data), std::end(data), [&](const auto& p_ptr)
    {
        return (p == *p_ptr) || box.squared_distance(p, *p_ptr) > minimum_distance_squared;
    }) != std::end(data);
}



ves::ParticleContainer::cartesian ves::ParticleContainer::getRandomValidPoint(REAL minimum_distance) const
{
    const REAL minimum_distance_squared = minimum_distance * minimum_distance;

    auto random_point = box.randomPointInside();

    while( std::find_if_not(begin(), end(), [&](const auto& p_ptr)
    {
        return box.squared_distance(random_point, p_ptr->getCoordinates()) > minimum_distance_squared;
    }) != end() )
    {
        random_point = box.randomPointInside();
    }
    return random_point;

}



auto ves::ParticleContainer::getClosestParticle(const cartesian& c) const -> decltype(std::begin(data))
{
    tbb::concurrent_vector<REAL> distances(data.size());
    tbb::parallel_for(std::size_t(0), data.size(), std::size_t(1), [&](std::size_t i){ distances[i] = box.squared_distance(c, data[i]->getCoordinates()); });
    return std::begin(data) + std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));
}



void ves::ParticleContainer::shiftAll(const cartesian& shift_vector)
{
    tbb::parallel_for_each(std::begin(data), std::end(data), [&shift_vector](const auto& particle_ptr)
    {
        if(!particle_ptr)
            vesCRITICAL("particle_ptr" << particle_ptr.get());
        particle_ptr->forcefullyShift(shift_vector);
    });
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



void ves::ParticleContainer::setupFromNew()
{

    switch (GLOBAL::getInstance().simulationmode.load())
    {
        case GLOBAL::SIMULATIONMODE::SA:
        {
            vesLOG("GLOBAL::SIMULATIONMODE::SA");
            try
            {
                REAL x = Parameters::getInstance().getOption("system.box.x").as<REAL>();
                REAL y = Parameters::getInstance().getOption("system.box.y").as<REAL>();
                REAL z = Parameters::getInstance().getOption("system.box.z").as<REAL>();
                box.setLengthX(x);
                box.setLengthY(y);
                box.setLengthZ(z);
            }
            catch(std::logic_error& e)
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
                catch(std::logic_error& e)
                {
                    vesCRITICAL("unable to get box from density|mobile, abort");
                }
                vesLOG("box(x|y|z) = (" << box.getLengthX() << "|" << box.getLengthY() << "|" << box.getLengthZ() << ")");
                
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.x", boost::program_options::variable_value(box.getLengthX(), false)));
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.y", boost::program_options::variable_value(box.getLengthY(), false)));
                Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.z", boost::program_options::variable_value(box.getLengthZ(), false)));
            }

            std::size_t mobile;
            if(Parameters::getInstance().getOptions().count("system.mobile"))
            {
                mobile = Parameters::getInstance().getOption("system.mobile").as<std::size_t>();
            }
            else if(Parameters::getInstance().getOptions().count("system.density"))
            {
                mobile = std::round(box.getVolume()*Parameters::getInstance().getOption("system.density").as<REAL>());
            }
            else
            {
                vesCRITICAL( "unable to obtain number of mobile particles");
            }

            for(std::size_t i = 0; i < mobile; ++i)
            {
                const auto minimum_offset = Parameters::getInstance().getOption("system.ljsigma").as<REAL>();
                auto particle_it = addParticle<Particle::TYPE::MOBILE>();
                do
                {
                    particle_it->get()->try_setCoordinates(box.randomPointInside());
                }
                while(placement_conflict(*(particle_it->get()), minimum_offset));
                while(!particle_it->get()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized())){;}
            }
            

            break;
        }

        case GLOBAL::SIMULATIONMODE::FGA:
        {
            vesLOG("GLOBAL::SIMULATIONMODE::FGA");
            std::size_t mobile;
            try
            {
                mobile = Parameters::getInstance().getOption("system.mobile").as<std::size_t>();
            }
            catch(const std::exception& e)
            {
                vesCRITICAL(e.what());
            }
            
            const std::size_t frame_guides_grid_edge = Parameters::getInstance().getOption("system.frame_guides_grid_edge").as<std::size_t>();
            const std::size_t guiding_elements_each = Parameters::getInstance().getOption("system.guiding_elements_each").as<std::size_t>();
            const REAL density = Parameters::getInstance().getOption("system.density").as<REAL>();
            const REAL ljsigma = Parameters::getInstance().getOption("system.ljsigma").as<REAL>();
            const REAL gamma = Parameters::getInstance().getOption("system.gamma").as<REAL>();

            try
            {
                REAL x = Parameters::getInstance().getOption("system.box.x").as<REAL>();
                REAL y = Parameters::getInstance().getOption("system.box.y").as<REAL>();
                REAL z = Parameters::getInstance().getOption("system.box.z").as<REAL>();
                box.setLengthX(x);
                box.setLengthY(y);
                box.setLengthZ(z);
            }
            catch(std::logic_error& e)
            {
                vesLOG("unable to get box from x|y|z, will calculate cubic box from density and particles");
                try
                {
                    switch(GLOBAL::getInstance().fgamode.load())
                    {
                        case GLOBAL::FGAMODE::SPHERE:
                        {
                            vesLOG("GLOBAL::FGAMODE::SPHERE");
                            const std::size_t particles_complete = mobile + frame_guides_grid_edge*guiding_elements_each;
                            const REAL volume = static_cast<REAL>(particles_complete)/density;
                            const REAL edge = std::cbrt(volume);
                            box.setLengthX(edge);
                            box.setLengthY(edge);
                            box.setLengthZ(edge);
                            break;
                        }

                        case GLOBAL::FGAMODE::PLANE:
                        {
                            vesLOG("GLOBAL::FGAMODE::PLANE");
                            const std::size_t particles_complete = mobile + guiding_elements_each;
                            const REAL volume = static_cast<REAL>(particles_complete)/density;
                            const REAL edge = std::cbrt(volume);
                            box.setLengthX(edge);
                            box.setLengthY(edge);
                            box.setLengthZ(edge);
                            break;
                        }
                        
                        case GLOBAL::FGAMODE::PAIR:
                        {
                            vesLOG("GLOBAL::FGAMODE::PAIR");
                            const std::size_t particles_complete = mobile + guiding_elements_each;
                            const REAL volume = static_cast<REAL>(particles_complete)/density;
                            const REAL edge = std::cbrt(volume);
                            box.setLengthX(edge);
                            box.setLengthY(edge);
                            box.setLengthZ(edge);
                            break;
                        }

                        default: 
                            vesCRITICAL("encountered invalde FGAMODE during box setup");
                            break;
                    }
                }
                catch(std::logic_error& e)
                {
                    vesCRITICAL("unable to get box from density|mobile+guiding_elements, abort  ");
                }
            }

            Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.x", boost::program_options::variable_value(box.getLengthX(), false)));
            Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.y", boost::program_options::variable_value(box.getLengthY(), false)));
            Parameters::getInstance().mutableAccess().insert(std::make_pair("system.box.z", boost::program_options::variable_value(box.getLengthZ(), false)));

            if(Parameters::getInstance().getOptions().count("system.mobile"))
            {
                mobile = Parameters::getInstance().getOption("system.mobile").as<std::size_t>();
            }
            else if(Parameters::getInstance().getOptions().count("system.density"))
            {
                mobile = std::round(box.getVolume()*Parameters::getInstance().getOption("system.density").as<REAL>());
            }
            else
            {
                vesCRITICAL( "unable to obtain number of mobile particles");
            }

            switch(GLOBAL::getInstance().fgamode.load())
            {
                case GLOBAL::FGAMODE::SPHERE:
                {
                    vesLOG("GLOBAL::FGAMODE::SPHERE");

                    const REAL radius = std::pow(ljsigma,1.0/6.0)/(2.0*std::sin(gamma));
                    const REAL dist_x = box.getLengthX()/frame_guides_grid_edge;
                    const REAL dist_y = box.getLengthY()/frame_guides_grid_edge;
                    const REAL dist_z = box.getLengthZ()/frame_guides_grid_edge;

                    vesLOG("generate sphere grid with edge " << frame_guides_grid_edge << " and " << guiding_elements_each <<" points each" );
                    auto spheregrid = ves::SphereGridGeometry(frame_guides_grid_edge, radius, guiding_elements_each);
                    spheregrid.scale(cartesian(dist_x, dist_y, dist_z));
                    spheregrid.shift(cartesian(dist_x/2, dist_y/2, dist_z/2));

                    for(const auto& sphere : spheregrid.spheres)
                    {
                        for(const auto& point : sphere.points)
                        {
                            auto particle_it = addParticle<Particle::TYPE::FRAME>();
                            particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(point));
                            particle_it->get()->coordinates_bounding.setBoundingSphere(0, ljsigma);
                            bool worked = particle_it->get()->try_setCoordinates(point);
                            if(!worked)
                                vesCRITICAL("setting Particle::Frame coordinates to point did not work");

                            // particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(sphere.origin));
                            worked = particle_it->get()->try_setOrientation(point-sphere.origin);
                            particle_it->get()->orientation_bounding.value = PI_4;
                            if(!worked)
                                vesCRITICAL("setting Particle::Frame orientation to point did not work");
                        }
                    }
                    for(std::size_t i = 0; i < mobile; ++i)
                    {
                        const auto minimum_offset = ljsigma;
                        auto particle_it = addParticle<Particle::TYPE::MOBILE>();
                        do
                        {
                            particle_it->get()->try_setCoordinates(box.randomPointInside());
                        }
                        while(placement_conflict(*(particle_it->get()), minimum_offset));
                        while(!particle_it->get()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized())){;}
                        // vesLOG("placed particle: " << particle_it->get()->getCoordinates().format(ROWFORMAT));
                    }
                    break;
                }

                case GLOBAL::FGAMODE::PLANE:
                {
                    vesLOG("GLOBAL::FGAMODE::PLANE");

                    const std::size_t guiding_elements_per_dimension = std::sqrt(guiding_elements_each);
                    const auto plane_edge = Parameters::getInstance().getOption("system.plane_edge").as<REAL>();
                    const auto plane_edge_half = plane_edge/2;
                    const auto scaling_factor = plane_edge / (guiding_elements_per_dimension);
                    // const auto scaling_factor = plane_edge / (guiding_elements_per_dimension-1);
                    const auto shift_vec = cartesian(box.getLengthX()/2-plane_edge_half+scaling_factor/2, box.getLengthY()/2-plane_edge_half+scaling_factor/2, box.getLengthZ()/2);
                    // const auto shift_vec = cartesian(box.getLengthX()/2-plane_edge_half, box.getLengthY()/2-plane_edge_half, box.getLengthZ()/2);

                    if(guiding_elements_per_dimension*guiding_elements_per_dimension != guiding_elements_each)
                    {
                        vesCRITICAL("guiding_elements_per_dimension != guiding_elements_each");
                    }

                    ves::PlaneGeometry plane(guiding_elements_per_dimension, guiding_elements_per_dimension);
                    plane.scale(cartesian(scaling_factor, scaling_factor, 0));
                    plane.shift(shift_vec);
                    
                    vesLOG("will setup plane with " << plane.points.size() << " particles and initial distance of " << scaling_factor);

                    for(const auto& point : plane.points)
                    {
                        vesLOG(point.format(ROWFORMAT));
                        auto particle_it = addParticle<Particle::TYPE::FRAME>();
                        if(Parameters::getInstance().getOption("system.guiding_elements_restriction").as<std::string>() == "inplace")
                        {
                            vesDEBUG("binding guiding element to sphere from " << 0 << " to " << ljsigma);
                            particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(point));
                            particle_it->get()->coordinates_bounding.setBoundingSphere(0, ljsigma);
                        }
                        else if(Parameters::getInstance().getOption("system.guiding_elements_restriction").as<std::string>() == "structure")
                        {
                            const auto box_min = cartesian(box.getLengthX()/2-plane_edge_half, box.getLengthY()/2-plane_edge_half, box.getLengthZ()/2-ljsigma/2);
                            const auto box_max = cartesian(box.getLengthX()/2+plane_edge_half, box.getLengthY()/2+plane_edge_half, box.getLengthZ()/2+ljsigma/2);
                            const ves::Particle::Base::box3d frame_bounds(box_min, box_max);
                            vesDEBUG("binding guiding element to box from " << box_min.format(ROWFORMAT) << " to " << box_max.format(ROWFORMAT));
                            particle_it->get()->coordinates_bounding.setBoundingBox(ves::Particle::Base::box3d(frame_bounds));
                        }
                        else
                            vesCRITICAL("could not set bounding to guiding element");

                        bool worked = particle_it->get()->try_setCoordinates(point);
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame coordinates to point did not work");

                        // particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(sphere.origin));
                        worked = particle_it->get()->try_setOrientation(cartesian::UnitZ());
                        particle_it->get()->orientation_bounding.value = PI_4;
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame orientation to point did not work");
                    }
                    for(std::size_t i = 0; i < mobile; ++i)
                    {
                        const auto minimum_offset = ljsigma;
                        auto particle_it = addParticle<Particle::TYPE::MOBILE>();

                        do particle_it->get()->try_setCoordinates(box.randomPointInside());
                        while(placement_conflict(*(particle_it->get()), minimum_offset));

                        while(!particle_it->get()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized())){;}
                        // vesLOG("placed particle: " << particle_it->get()->getCoordinates().format(ROWFORMAT));
                    }

                    break;
                }

                case GLOBAL::FGAMODE::TUBE:
                {
                    vesLOG("GLOBAL::FGAMODE::TUBE");

                    auto tube = ves::TubeGeometry();
                    tube.height = box.getLengthZ()/(guiding_elements_each+1)*guiding_elements_each;
                    tube.size = guiding_elements_each;
                    tube.radius = ljsigma * std::pow(2,1.f/6) / enhance::deg_to_rad(15.f);
                    tube.generate();
                    tube.shift(cartesian(0, 0, box.getLengthZ()/(guiding_elements_each)/2));

                    tube.shift(cartesian(box.getLengthX()/2, box.getLengthY()/2, 0));
                    for(const auto& point : tube.points)
                    {
                        auto particle_it = addParticle<Particle::TYPE::FRAME>();
                        particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(point);
                        particle_it->get()->coordinates_bounding.setBoundingSphere(0, ljsigma);
                        bool worked = particle_it->get()->try_setCoordinates(point);
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame coordinates to point did not work");

                        const auto center = cartesian(box.getLengthX()/2, box.getLengthY()/2, point(2));
                        worked = particle_it->get()->try_setOrientation(point-center);
                        particle_it->get()->orientation_bounding.value = PI_4;
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame orientation to point did not work");
                    }

                    for(std::size_t i = 0; i < mobile; ++i)
                    {
                        const auto minimum_offset = ljsigma;
                        auto particle_it = addParticle<Particle::TYPE::MOBILE>();
                        do
                        {
                            particle_it->get()->try_setCoordinates(box.randomPointInside());
                        }
                        while(placement_conflict(*(particle_it->get()), minimum_offset));
                        while(!particle_it->get()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized())){;}
                    }
                    
                    break;
                }

                default:
                    vesCRITICAL("encountered invalde FGAMODE during particle generation");
                    break;
                case GLOBAL::FGAMODE::PAIR:
                {
                    vesLOG("GLOBAL::FGAMODE::TUBE");

                    const auto plane_edge = Parameters::getInstance().getOption("system.plane_edge").as<REAL>();
                    const auto center = box.getCenter();
                    const auto orientation_bound = PI_4/2;

                    {
                        const auto left = cartesian(center(0)-plane_edge/2, center(1), center(2));
                        auto particle_it = addParticle<Particle::TYPE::FRAME>();
                        particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(left);
                        particle_it->get()->coordinates_bounding.setBoundingSphere(0, 0);
                        bool worked = particle_it->get()->try_setCoordinates(left);
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame coordinates to point did not work");

                        worked = particle_it->get()->try_setOrientation(cartesian(0,0,1));
                        particle_it->get()->orientation_bounding.value = orientation_bound;
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame orientation to point did not work");
                    }

                    {
                        const auto right = cartesian(center(0)+plane_edge/2, center(1), center(2));
                        auto particle_it = addParticle<Particle::TYPE::FRAME>();
                        particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(right);
                        particle_it->get()->coordinates_bounding.setBoundingSphere(0, 0);
                        bool worked = particle_it->get()->try_setCoordinates(right);
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame coordinates to point did not work");

                        worked = particle_it->get()->try_setOrientation(cartesian(0,0,1));
                        particle_it->get()->orientation_bounding.value = orientation_bound;
                        if(!worked)
                            vesCRITICAL("setting Particle::Frame orientation to point did not work");
                    }

                    for(std::size_t i = 0; i < mobile; ++i)
                    {
                        const auto minimum_offset = ljsigma;
                        auto particle_it = addParticle<Particle::TYPE::MOBILE>();
                        do
                        {
                            particle_it->get()->try_setCoordinates(box.randomPointInside());
                        }
                        while(placement_conflict(*(particle_it->get()), minimum_offset));
                        while(!particle_it->get()->try_setOrientation(ves::Particle::Base::cartesian::Random().normalized())){;}
                    }
                    
                    break;
                }
            }
            break;
        }

        case GLOBAL::SIMULATIONMODE::OSMOTIC:
        {
            vesLOG("GLOBAL::SIMULATIONMODE::OSMOTIC");

            const std::size_t mobile = Parameters::getInstance().getOption("system.mobile").as<std::size_t>();
            const std::size_t guiding_elements_each = Parameters::getInstance().getOption("system.guiding_elements_each").as<std::size_t>();
            const REAL ljsigma = Parameters::getInstance().getOption("system.ljsigma").as<REAL>();
            const REAL gamma = Parameters::getInstance().getOption("system.gamma").as<REAL>();
            const REAL kappa = Parameters::getInstance().getOption("system.kappa").as<REAL>();

            try
            {
                REAL x = Parameters::getInstance().getOption("system.box.x").as<REAL>();
                REAL y = Parameters::getInstance().getOption("system.box.y").as<REAL>();
                REAL z = Parameters::getInstance().getOption("system.box.z").as<REAL>();
                box.setLengthX(x);
                box.setLengthY(y);
                box.setLengthZ(z);
            }
            catch(std::logic_error& e)
            {
                vesCRITICAL("unable to get box from x|y|z");
            }

            const REAL radius = enhance::nth_root<6>(ljsigma*2)/(2.0*std::sin(gamma)) + kappa;
            vesLOG("construct sphere: center " << box.getCenter().format(ROWFORMAT) << "  | radius " << radius << "  | " <<  mobile+guiding_elements_each);
            ves::SphereGeometry sphere(box.getCenter(), radius, mobile+guiding_elements_each+1);

            auto guiding_elements_to_place = guiding_elements_each;
            auto mobile_to_place = mobile;
            vesLOG("construct sphere from " << guiding_elements_to_place << " guiding elements and " << mobile_to_place << " mobile particles");
            std::size_t second_skipper = 0;
            for(const auto& point : sphere.points)
            {
                if(second_skipper++ == 1)
                    continue;

                if(guiding_elements_to_place > 0)
                {
                    --guiding_elements_to_place;
                    auto particle_it = addParticle<Particle::TYPE::FRAME>();
                    particle_it->get()->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(point));
                    particle_it->get()->coordinates_bounding.setBoundingSphere(0, ljsigma);
                    bool worked = particle_it->get()->try_setCoordinates(point);
                    // vesLOG(particle_it->get()->getCoordinates().format(ROWFORMAT));
                    // vesLOG(particle_it->get()->coordinates_bounding.origin->format(ROWFORMAT));
                    if(!worked)
                    {
                        vesCRITICAL("setting Particle::Frame coordinates to point did not work");
                    }

                    worked = particle_it->get()->try_setOrientation(point-sphere.origin);
                    particle_it->get()->orientation_bounding.value = PI_4;
                    if(!worked)
                        vesCRITICAL("setting Particle::Frame orientation to point did not work");
                }
                else if(mobile_to_place > 0)
                {
                    --mobile_to_place;
                    auto particle_it = addParticle<Particle::TYPE::MOBILE>();
                    bool worked = particle_it->get()->try_setCoordinates(point);
                    if(!worked)
                        vesCRITICAL("setting Particle::Mobile coordinates to point did not work");

                    worked = particle_it->get()->try_setOrientation(point-sphere.origin);
                    if(!worked)
                        vesCRITICAL("setting Particle::Mobile orientation to point did not work");
                }
                else
                    vesCRITICAL("nothing to place");
            }

            const REAL inner_radius = radius-ljsigma;
            const REAL outer_radius = radius+ljsigma;
            const std::size_t osmotic_inside = Parameters::getInstance().getOption("system.osmotic_density_inside").as<REAL>() * (enhance::sphere_volume(inner_radius));
            const std::size_t osmotic_outside = Parameters::getInstance().getOption("system.osmotic_density_outside").as<REAL>() * (box.getVolume() - enhance::sphere_volume(outer_radius));
            vesLOG("place osmotic_inside " << osmotic_inside << " and osmotic_outside " << osmotic_outside);
            vesLOG("inner radius " << inner_radius);
            vesLOG("outer radius " << outer_radius);

            for(std::size_t i = 0; i < osmotic_inside; ++i)
            {
                auto particle_it = addParticle<Particle::TYPE::OSMOTIC>();
                // auto point = box.randomPointInside();
                auto point = cartesian(box.getCenter() + cartesian(enhance::random<REAL>()(-inner_radius, inner_radius), enhance::random<REAL>()(-inner_radius, inner_radius), enhance::random<REAL>()(-inner_radius, inner_radius)));
                do
                {
                    point = cartesian(box.getCenter() + cartesian(enhance::random<REAL>()(-inner_radius, inner_radius), enhance::random<REAL>()(-inner_radius, inner_radius), enhance::random<REAL>()(-inner_radius, inner_radius)));
                    particle_it->get()->try_setCoordinates(point);
                }
                while((point-box.getCenter()).norm() > inner_radius || placement_conflict(*(particle_it->get()), ljsigma));
                // vesLOG("placed particle inside" << particle_it->get()->getCoordinates().format(ROWFORMAT) << " with distance from center " << (point-box.getCenter()).norm());
            }

            for(std::size_t i = 0; i < osmotic_outside; ++i)
            {
                auto particle_it = addParticle<Particle::TYPE::OSMOTIC>();
                auto point = box.randomPointInside();
                do
                {
                    point = box.randomPointInside();
                    particle_it->get()->try_setCoordinates(point);
                }
                while((point-box.getCenter()).norm() < outer_radius || placement_conflict(*(particle_it->get()), ljsigma));
                // vesLOG("placed particle outside" << particle_it->get()->getCoordinates().format(ROWFORMAT) << " with distance from center " << (point-box.getCenter()).norm());
            }

            Parameters::getInstance().mutableAccess().insert(std::make_pair("system.density", boost::program_options::variable_value(REAL(data.size())/box.getVolume(), false)));
            break;
        }

        default:
            vesCRITICAL("encountered invalid SIMULATIONMODE");
            break;
    }
}



std::vector<std::string> ves::ParticleContainer::getGroupNames()
{
    vesLOG(__PRETTY_FUNCTION__);
    auto file_info = [](hid_t loc_id, const char *name, [[maybe_unused]] const H5L_info_t * linfo, void *opdata) -> herr_t
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
    try{ group_names.erase(std::remove_if(group_names.begin(), group_names.end(), [](auto a){ return a == std::string("position_sphere_bounds");}), group_names.end()); } catch(...){;}
    try{ group_names.erase(std::remove_if(group_names.begin(), group_names.end(), [](auto a){ return a == std::string("position_box_bounds");}), group_names.end()); } catch(...){;}
    try{ group_names.erase(std::remove_if(group_names.begin(), group_names.end(), [](auto a){ return a == std::string("orientation_bounds");}), group_names.end()); } catch(...){;}

    const std::string numbers = "0123456789";
    std::sort(std::begin(group_names), std::end(group_names), [&](const std::string& a, const std::string& b) 
    { 
        std::string cmp_a { a.begin() + a.find_first_of(numbers), a.end() };
        std::string cmp_b { b.begin() + b.find_first_of(numbers), b.end() };
        return std::stoi(cmp_a) < std::stoi(cmp_b);
    });
    // group_names.sort();
    assert(std::is_sorted(std::begin(group_names), std::end(group_names)));
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
    vesLOG("get group " << latest_group_name);
    h5xx::group group(h5file, latest_group_name);
    h5xx::group rootgroup(h5file, "/");

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
    
    // if(GLOBAL::getInstance().ensemble == GLOBAL::ENSEMBLE::uVT)
    //     Parameters::getInstance().mutableAccess().insert(std::make_pair("constant.mu", boost::program_options::variable_value(h5xx::read_attribute<REAL>(rootgroup, "constant.mu"), false)));
    
    time_from_input = h5xx::read_attribute<decltype(time_from_input)>(group, "system.actual_time")+1;
    
    array1d_t types;
    array1d_t IDs;
    array2d_t positions;
    array2d_t orientations;
    array2d_t position_sphere_bounds;
    array2d_t position_box_bounds;
    array2d_t orientation_bounds;

    {
        vesLOG("load type");
        h5xx::dataset dataset(group, "type");
        h5xx::read_dataset(dataset, types);
    }

    {
        vesLOG("load pos");
        h5xx::dataset dataset(group, "position");
        h5xx::read_dataset(dataset, positions);
    }

    {
        vesLOG("load orien");
        h5xx::dataset dataset(group, "orientation");
        h5xx::read_dataset(dataset, orientations);
    }

    try
    {
        vesLOG("try load position_sphere_bounds");
        h5xx::dataset dataset(rootgroup, "position_sphere_bounds");
        h5xx::read_dataset(dataset, position_sphere_bounds);
    }catch(...){ vesLOG("couldnt find position_sphere_bounds");}

    try
    {
        vesLOG("try load position_box_bounds");
        h5xx::dataset dataset(rootgroup, "position_box_bounds");
        h5xx::read_dataset(dataset, position_box_bounds);
    }catch(...){ vesLOG("couldnt find position_box_bounds");}

    try
    {
        vesLOG("try load orientation_bounds");
        h5xx::dataset dataset(rootgroup, "orientation_bounds");
        h5xx::read_dataset(dataset, orientation_bounds);
    }catch(...){ vesLOG("couldnt find orientation_bounds");}

    // generate particles
    for(std::size_t i = 0; i < types.size() ; ++i)
    {
        switch (types[i])
        {
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::MOBILE):
            {
                [[maybe_unused]] auto particle_it = addParticle<ves::Particle::TYPE::MOBILE>();
                break;
            }
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::FRAME):
            {
                [[maybe_unused]] auto particle_it = addParticle<ves::Particle::TYPE::FRAME>();
                break;
            }
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::OSMOTIC):
            {
                [[maybe_unused]] auto particle_it = addParticle<ves::Particle::TYPE::OSMOTIC>();
                break;
            }
        
            default:
            {
                throw std::logic_error("undefined particle type in TrajectoryDistributorH5 in row "+std::to_string(i));
                break;
            }
        }
    }

    // set 'em up
    for(std::size_t i = 0; i < types.size() ; ++i)
    {
        vesDEBUG(i << cartesian(positions[i][0], positions[i][1], positions[i][2]).format(ROWFORMAT) << cartesian(orientations[i][0], orientations[i][1], orientations[i][2]).format(ROWFORMAT));
        

        switch (types[i])
        {
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::MOBILE):
            {
                data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);
                data.at(i)->getOrientation() = cartesian(orientations[i][0], orientations[i][1], orientations[i][2]);
                break;
            }
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::FRAME):
            {
                if(position_sphere_bounds.size() > 0)
                {
                    data.at(i)->coordinates_bounding.origin = std::make_unique<cartesian>(cartesian(position_sphere_bounds[i][0], position_sphere_bounds[i][1], position_sphere_bounds[i][2]));
                    data.at(i)->coordinates_bounding.setBoundingSphere(position_sphere_bounds[i][3], position_sphere_bounds[i][4]);
                    data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);
                }
                else if(position_box_bounds.size() > 0)
                {
                    const auto min = cartesian(position_box_bounds[i][0], position_box_bounds[i][1], position_box_bounds[i][2]);
                    const auto max = cartesian(position_box_bounds[i][3], position_box_bounds[i][4], position_box_bounds[i][5]);
                    data.at(i)->coordinates_bounding.setBoundingBox(ves::Particle::Base::box3d(min, max));
                    data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);

                }
                else
                {
                    vesWARNING("setting Frame without position bounds");
                    data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);
                }

                if(orientation_bounds.size() > 0)
                {   
                    data.at(i)->orientation_bounding.origin = std::make_unique<cartesian>(cartesian(orientation_bounds[i][0], orientation_bounds[i][1], orientation_bounds[i][2]));
                    data.at(i)->orientation_bounding.value = orientation_bounds[i][3];
                    data.at(i)->getOrientation() = cartesian(orientations[i][0], orientations[i][1], orientations[i][2]);
                }
                else
                {
                    vesWARNING("setting Frame without orientation bounds");
                    data.at(i)->getOrientation() = cartesian(orientations[i][0], orientations[i][1], orientations[i][2]);
                }

                break;
            }
            case static_cast<std::underlying_type<ves::Particle::TYPE>::type>(ves::Particle::TYPE::OSMOTIC):
            {
                data.at(i)->getCoordinates() = cartesian(positions[i][0], positions[i][1], positions[i][2]);
                break;
            }
        
            default:
                throw std::logic_error("undefined particle type in TrajectoryDistributorH5 in row "+std::to_string(i));
                break;
        }

    }

    h5file.close();
}