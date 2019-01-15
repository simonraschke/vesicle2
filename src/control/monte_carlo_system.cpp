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
    box.setLengthX(Parameters::getInstance().getOption("system.box.x").as<REAL>());
    box.setLengthY(Parameters::getInstance().getOption("system.box.y").as<REAL>());
    box.setLengthZ(Parameters::getInstance().getOption("system.box.z").as<REAL>());
    cells.setup();
    cells.deployParticles(particles);
    interaction = std::make_unique<AngularLennardJonesInteraction>();
    sw_position.setup
    (
        1000, 
        ves::Parameters::getInstance().getOption("system.sw_position_min").as<REAL>(), 
        ves::Parameters::getInstance().getOption("system.sw_position_max").as<REAL>(),
        ves::Parameters::getInstance().getOption("system.acceptance_position_target").as<REAL>()
    );
    sw_orientation.setup
    (
        1000, 
        ves::Parameters::getInstance().getOption("system.sw_orientation_min").as<REAL>(), 
        ves::Parameters::getInstance().getOption("system.sw_orientation_max").as<REAL>(),
        ves::Parameters::getInstance().getOption("system.acceptance_orientation_target").as<REAL>()
    );

    traj_gro.setup();
    traj_h5.setup();
}



void ves::MonteCarloSystem::run()
{
    GLOBAL::getInstance().simulationstatus.store(GLOBAL::SIMULATIONSTATUS::RUNNING);

    for(; time < time_max; ++time)
    {
        if(ves::Controller::SIGNAL.load() != 0)
        {
            break;
        }
        vesDEBUG("prepare cells");
        cells.preparation();
        vesDEBUG("do MC step");
        cells.cellBasedApplyFunctor([&](const auto& c){ cellStep(c);});
        vesDEBUG("reorder cells");
        cells.reorder();
        assert(cells.membersContained()==particles.data.size());


        const auto num_members_is = cells.membersContained();
        const auto num_members_should = particles.data.size();
        if(num_members_is!=num_members_should)
        {
            // Backtrace();
            vesCRITICAL("lost particles while reordering: " << std::boolalpha << (num_members_is==num_members_should) << " found " << num_members_is << " should be " << num_members_should )
        }

        if(time % output_skip == 0)
        {
            vesLOG(status().str());
            // vesLOG(time << " took " << duration << " s,  acceptance: " << sw_position.getRatio() << "  " << sw_orientation.getRatio() << "  sw: " << sw_position() << "  " << sw_orientation());
            traj_gro.write(*this);
            traj_h5.write(*this);
        }
    }

}



void ves::MonteCarloSystem::cellStep(const ves::Cell& cell)
{
    // vesDEBUG(__PRETTY_FUNCTION__ << std::addressof(cell));
    REAL last_energy_value;
    
    for(const Cell::particle_ptr_t& particle : cell)
    {
        last_energy_value = 0;
        // coordinates move
        {
            const auto stepwidth = sw_position();
            std::uniform_real_distribution<REAL> dist(-stepwidth,stepwidth);
            const auto translation = Particle::Base::cartesian
            (
                dist(pseudo_engine),
                dist(pseudo_engine),
                dist(pseudo_engine)
            );

            last_energy_value = cell.potential(*particle);
            if(particle->try_setCoordinates(particle->getCoordinates()+translation))
            {
                const REAL energy_after = cell.potential(*particle);

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

        // orientation move
        {
            const REAL stepwidth = sw_orientation();
            const auto random_vec = Particle::Base::cartesian::Random();
            const Eigen::AngleAxis<REAL> rotation (stepwidth, random_vec);

            if(const auto orientation_before = particle->getOrientation(); particle->try_setOrientation(rotation * particle->getOrientation()))
            {
                const REAL energy_after = cell.potential(*particle);

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
