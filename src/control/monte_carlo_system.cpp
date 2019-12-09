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

#include "monte_carlo_system.hpp"



thread_local std::mt19937_64 ves::MonteCarloSystem::pseudo_engine(std::random_device().operator()());



ves::MonteCarloSystem::MonteCarloSystem()
    : time_max(ves::Parameters::getInstance().getOption("system.time_max").as<std::size_t>())
    , output_skip(ves::Parameters::getInstance().getOption("output.skip").as<std::size_t>())
    , widom_skip(ves::Parameters::getInstance().getOption("system.widom_skip").as<std::size_t>())
{

}



auto ves::MonteCarloSystem::status()
{
    static const auto milliseconds_per_day = (24*60*60*1000);
    static auto start = std::chrono::high_resolution_clock::now();
    const auto now = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration<double>(std::chrono::duration_cast<std::chrono::milliseconds>(now - start)).count();
    const auto time_per_step = duration/output_skip;

    std::stringstream ss;
    ss << std::fixed;
    ss << "step: " << std::setw(10) << std::setprecision(3) << std::right << getTime()  << " | ";
    ss << "particles: " << std::setw(7) << std::setprecision(3) << std::right << particles.data.size() << "  | ";
    ss << "energy: " << std::setw(10) << std::setprecision(3) << std::right << potential() << "  | ";
    ss << "time: " << std::setw(7) << std::setprecision(3) << std::right << duration << " s  | ";
    ss << " /particle/step: " << std::setw(7) << std::setprecision(3) << std::right << time_per_step/particles.data.size()*1e6 << " ns  | ";
    ss << " /cell/step: " << std::setw(7) << std::setprecision(3) << std::right << time_per_step/cells.data.size()*1e6 << " ns  | ";
    ss << " per day: " << std::setw(10) << std::setprecision(2) << std::right << std::scientific << milliseconds_per_day/duration << " steps  | ";
    ss << " to goal: " << std::setw(6) << std::setprecision(2) << std::right << std::fixed << (time_max-getTime())*time_per_step/24/60/60 << " d  | ";
    ss << " sw_coord: " << std::setw(7) << std::setprecision(3) << std::right << sw_position() << " (" << std::setprecision(0) << std::round(sw_position.getRatio()*100) << "%)  | ";
    ss << " sw_orien: " << std::setw(7) << std::setprecision(3) << std::right << sw_orientation() << " (" << std::setprecision(0) <<std::round(sw_orientation.getRatio()*100) << "%)  | ";

    start = now;
    return ss;
}



void ves::MonteCarloSystem::setup()
{
    particles.setup();
    time = particles.getTime();
    box.setLengthX(Parameters::getInstance().getOption("system.box.x").as<REAL>());
    box.setLengthY(Parameters::getInstance().getOption("system.box.y").as<REAL>());
    box.setLengthZ(Parameters::getInstance().getOption("system.box.z").as<REAL>());
    cells.setup();
    cells.deployParticles(particles);
    interaction = std::make_unique<AngularLennardJonesInteraction>();

    if(GLOBAL::getInstance().startmode == GLOBAL::STARTMODE::NEW)
    {
        sw_position.setup
        (
            1000*particles.data.size(), 
            ves::Parameters::getInstance().getOption("system.sw_position_min").as<REAL>(), 
            ves::Parameters::getInstance().getOption("system.sw_position_max").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.acceptance_position_target").as<REAL>()
        );
        sw_orientation.setup
        (
            1000*particles.data.size(), 
            ves::Parameters::getInstance().getOption("system.sw_orientation_min").as<REAL>(), 
            ves::Parameters::getInstance().getOption("system.sw_orientation_max").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.acceptance_orientation_target").as<REAL>()
        );
    }
    else if(GLOBAL::getInstance().startmode == GLOBAL::STARTMODE::RESTART)
    {
        sw_position.setup
        (
            1000*particles.data.size(), 
            ves::Parameters::getInstance().getOption("system.sw_position_min").as<REAL>(), 
            ves::Parameters::getInstance().getOption("system.sw_position_max").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.sw_position_actual").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.acceptance_position_target").as<REAL>()
        );
        sw_orientation.setup
        (
            1000*particles.data.size(), 
            ves::Parameters::getInstance().getOption("system.sw_orientation_min").as<REAL>(), 
            ves::Parameters::getInstance().getOption("system.sw_orientation_max").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.sw_orientation_actual").as<REAL>(),
            ves::Parameters::getInstance().getOption("system.acceptance_orientation_target").as<REAL>()
        );
    }

    if(! Parameters::getInstance().getOption("output.suppress_gro").as<bool>())
    {
        traj_gro.setup();
    }

    traj_h5.setup(*this);

    if(time == 0)
    {
        if(! Parameters::getInstance().getOption("output.suppress_gro").as<bool>())
        {
            traj_gro.write(*this);
        }
        
        traj_h5.write(*this);
    }
}



void ves::MonteCarloSystem::run()
{
    GLOBAL::getInstance().simulationstatus.store(GLOBAL::SIMULATIONSTATUS::RUNNING);

    for(; time <= time_max; ++time)
    {
        if(ves::Controller::SIGNAL.load() != 0)
        {
            break;
        }
        vesDEBUG("prepare cells");
        cells.preparation();
        vesDEBUG("do MC step");
        cells.cellBasedApplyFunctor([&](const auto& c){ cellStep(c);});
        
        // const auto shift_vector = ves::Particle::Base::cartesian::Random()*Parameters::getInstance().getOption("system.ljsigma").as<REAL>();
        // if(GLOBAL::getInstance().simulationmode.load(std::memory_order_relaxed) == GLOBAL::SIMULATIONMODE::SA)
        // {
        //     particles.shiftAll(shift_vector);
        //     // vesLOG(shift_vector.format(ROWFORMAT));
        // }

        vesDEBUG("reorder cells");
        cells.reorder();
        assert(cells.membersContained()==particles.data.size());

        const auto num_members_is = cells.membersContained();
        const auto num_members_should = particles.data.size();
        if(num_members_is!=num_members_should)
        {
            vesCRITICAL("lost particles while reordering: " << std::boolalpha << (num_members_is==num_members_should) << " found " << num_members_is << " should be " << num_members_should )
        }

        if(time % widom_skip == 0 and time != 0 && GLOBAL::getInstance().ensemble.load(std::memory_order_relaxed) == GLOBAL::ENSEMBLE::uVT)
        {
            grandCanonicalStep();
        }

        if(time % output_skip == 0 or time == 0)
        {
            vesLOG(status().str());
        }

        if(time % output_skip == 0 and time != 0)
        {
            
            if(! Parameters::getInstance().getOption("output.suppress_gro").as<bool>())
            {
                traj_gro.write(*this);
            }
            traj_h5.write(*this);
        }
    }

}



void ves::MonteCarloSystem::cellStep(const ves::Cell& cell)
{
    // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
    // do not change order or modify any function call

    REAL last_energy_value;
    REAL energy_after;
    std::uniform_real_distribution<REAL> dist_coords(-sw_position(),sw_position());
    std::uniform_real_distribution<REAL> dist_orientation(-sw_orientation(),sw_orientation());
    Particle::Base::cartesian translation;
    Particle::Base::cartesian orientation_before;

    // std::vector<std::uint32_t> iteration_IDs(cell.data.size());
    // std::iota(std::begin(iteration_IDs), std::end(iteration_IDs), 0);
    // std::shuffle(std::begin(iteration_IDs), std::end(iteration_IDs), pseudo_engine);

    std::vector<ves::Cell::container_t::const_iterator> iterators(cell.data.size());
    std::iota(std::begin(iterators), std::end(iterators), std::begin(cell));
    std::shuffle(std::begin(iterators), std::end(iterators), pseudo_engine);


    // for(const auto num : iteration_IDs)
    // for(const Cell::particle_ptr_t& particle : iterators)
    for(const auto& iterator : iterators)
    {
        // const Cell::particle_ptr_t& particle = std::cref(cell.data[num]);
        const Cell::particle_ptr_t& particle = *iterator;
        
        // coordinates move
        {
            translation = Particle::Base::cartesian
            (
                dist_coords(pseudo_engine),
                dist_coords(pseudo_engine),
                dist_coords(pseudo_engine)
            );

            last_energy_value = cell.potentialOfSingleParticle(*particle);
            if(particle->try_setCoordinates(particle->getCoordinates()+translation))
            {
                energy_after = cell.potentialOfSingleParticle(*particle);

                // rejection
                if(!acceptance.isValid(energy_after - last_energy_value))
                {
                    sw_position.rejected();
                    particle->getCoordinates() -= translation;
                }
                // acctance
                else
                {
                    sw_position.accepted();
                    last_energy_value = energy_after;
                }
            }
            else
            {
                sw_position.rejected();
            }
        }

    // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
    // do not change order or modify any function call

        // orientation move
        if(particle->getType() != ves::Particle::TYPE::OSMOTIC)
        {
            // rotation = Eigen::AngleAxis<REAL>(dist_orientation(pseudo_engine), Particle::Base::cartesian::Random());

            if( 
                orientation_before = particle->getOrientation(); 
                particle->try_setOrientation(Eigen::AngleAxis<REAL>(dist_orientation(pseudo_engine), Particle::Base::cartesian::Random()) * particle->getOrientation())
            )
            {
                energy_after = cell.potentialOfSingleParticle(*particle);

                // rejection
                if(!acceptance.isValid(energy_after - last_energy_value))
                {
                    sw_orientation.rejected();
                    particle->getOrientation() = orientation_before;
                }
                // acctance
                else
                {
                    sw_orientation.accepted();
                    last_energy_value = energy_after;
                }
            }
            else
            {
                sw_orientation.rejected();
            }
        }
    }

    // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
    // do not change order or modify any function call
}



REAL ves::MonteCarloSystem::potential() const
{
    auto energy = std::make_shared<REAL>(0);
    tbb::parallel_for(std::size_t(0), particles.data.size(), std::size_t(1), [&](const std::size_t& i)
    {
        REAL pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
        {
            pre_sum += interaction->calculate(*particles.data[i], *particles.data[j]);
        }

        *energy += pre_sum;
    });
    return !std::isnan(*energy) ? *energy : throw std::runtime_error("potential Energy is NAN");
}



void ves::MonteCarloSystem::grandCanonicalStep()
{
    static std::uniform_real_distribution<REAL> dist(REAL(0), REAL(1));
    bool something_happened = false;

    if(dist(pseudo_engine) < REAL(0.5))
        something_happened = try_addParticle();
    else
        something_happened = try_removeParticle();

    if(something_happened)
    {
        sw_position.setAlignmentEvery(1000*particles.data.size());
        sw_orientation.setAlignmentEvery(1000*particles.data.size());
    }
}



bool ves::MonteCarloSystem::try_addParticle()
{
    auto random_point = box.randomPointInside();
    auto new_particle = ves::Particle::Mobile();
    new_particle.getCoordinates() = random_point;
    new_particle.getOrientation() = ves::Particle::Base::cartesian::Random().normalized();

    auto cell_it = cells.getCellOfCartesian(new_particle.getCoordinates());
    const auto energy_before = cell_it->potential();
    const auto energy_after  = cell_it->potentialWithPhantomParticle(new_particle);

    GrandCanonicalInsertion insertion;
    if(!insertion.isValid(*this, energy_after-energy_before))
    {
        return false;
    }
    else
    {
        vesDEBUG("added");
        auto particle_it = particles.addParticle<ves::Particle::TYPE::MOBILE>();
        particle_it->get()->getCoordinates() = new_particle.getCoordinates();
        particle_it->get()->getOrientation() = new_particle.getOrientation();
        cells.deployParticle(*particle_it->get());
        return true;
    }
}



bool ves::MonteCarloSystem::try_removeParticle()
{
    static const auto optimum_distance_squared = std::pow(enhance::nth_root<6>(Parameters::getInstance().getOption("system.ljsigma").as<REAL>()*2)/2, 2);

    const auto random_point = box.randomPointInside();
    const auto nearest = particles.getClosestParticle(random_point);

    if(nearest->get()->getType() != ves::Particle::TYPE::MOBILE)
        return false;

    bool is_a_particle = box.squared_distance(random_point, nearest->get()->getCoordinates()) <= optimum_distance_squared;
    if(is_a_particle)
    {
        auto cell_it = cells.getCellOfCartesian(nearest->get()->getCoordinates());
        const auto energy_before = cell_it->potential();
        const auto energy_after  = cell_it->potentialIgnoreParticle(*nearest->get());

        GrandCanonicalDeletion deletion;
        if(!deletion.isValid(*this, energy_after-energy_before))
        {
            return false;
        }
        else
        {
            vesDEBUG("removed");
            cells.removeParticle(*nearest->get());
            particles.removeParticle(*nearest->get());
            return true;
        }
    }
    else
        return false;
    
}

// bool ves::MonteCarloSystem::try_addParticle()
// {
//     // auto random_point = particles.getRandomValidPoint(ves::Parameters::getInstance().getOption("system.ljsigma").as<REAL>());
//     auto random_point = box.randomPointInside();
//     auto new_particle = ves::Particle::Mobile();
//     new_particle.getCoordinates() = random_point;
//     new_particle.getOrientation() = ves::Particle::Base::cartesian::Random().normalized();
//     auto cell_it = cells.getCellOfCartesian(new_particle.getCoordinates());
//     const auto energy_before = cell_it->potential();
//     const auto energy_after = cell_it->potentialWithPhantomParticle(std::ref(new_particle));
//     const REAL cavity_prob = REAL(1) - static_cast<REAL>(particles.data.size()) * enhance::sphere_volume(enhance::nth_root<6>(Parameters::getInstance().getOption("system.ljsigma").as<REAL>()*2)/2) / box.getVolume();

//     GrandCanonicalInsertion insertion;
//     if(!insertion.isValid(*this, cell_it->chemicalPotential(new_particle), cavity_prob, energy_after-energy_before))
//     {
//         return false;
//     }
//     else
//     {
//         vesDEBUG("added");
//         auto particle_it = particles.addParticle<ves::Particle::TYPE::MOBILE>();
//         particle_it->get()->getCoordinates() = new_particle.getCoordinates();
//         particle_it->get()->getOrientation() = new_particle.getOrientation();
//         cells.deployParticle(*particle_it->get());
//         return true;
//     }
// }



// bool ves::MonteCarloSystem::try_removeParticle()
// {
//     auto particle = particles.getRandomParticle<ves::Particle::TYPE::MOBILE>();
//     const auto particle_coords = particle->getCoordinates();
//     auto cell_it = cells.getCellOfCartesian(particle_coords);
//     const auto energy_before = cell_it->potential();
//     const auto energy_after = cell_it->potentialIgnoreParticle(*particle.get());
//     const REAL cavity_prob = REAL(1) - static_cast<REAL>(particles.data.size()) * enhance::sphere_volume(enhance::nth_root<6>(Parameters::getInstance().getOption("system.ljsigma").as<REAL>()*2)/2) / box.getVolume();

//     GrandCanonicalDeletion deletion;
//     if(!deletion.isValid(*this, cell_it->chemicalPotential(*particle.get()), cavity_prob, energy_after-energy_before))
//     {
//         return false;
//     }
//     else
//     {
//         vesDEBUG("removed");
//         cells.removeParticle(*particle.get());
//         particles.removeParticle(*particle.get());
//         return true;
//     }
// }