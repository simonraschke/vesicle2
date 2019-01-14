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



ves::MonteCarloSystem::MonteCarloSystem()
    : time_max(ves::Parameters::getInstance().getOption("system.time_max").as<std::size_t>())
    , output_skip(ves::Parameters::getInstance().getOption("output.skip").as<std::size_t>())
{

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
            vesLOG(time << " " << potential() << "  acceptance: " << sw_position.getRatio() << "  " << sw_orientation.getRatio() << "  sw: " << sw_position() << "  " << sw_orientation());
            traj_gro.write(*this);
        }
    }

}



void ves::MonteCarloSystem::cellStep(const ves::Cell& cell)
{
    // vesDEBUG(__PRETTY_FUNCTION__ << std::addressof(cell));
    for(const Cell::particle_ptr_t& particle : cell)
    {
        // coordinates move
        {
            const auto stepwidth = sw_position();
            const auto translation = Particle::Base::cartesian
            (
                enhance::random<REAL>(-stepwidth,stepwidth),
                enhance::random<REAL>(-stepwidth,stepwidth),
                enhance::random<REAL>(-stepwidth,stepwidth)
            );

            const REAL energy_before = cell.potential(*particle);
            // vesLOG("before " << particle->getCoordinates().format(ROWFORMAT));
            particle->try_setCoordinates(particle->getCoordinates()+translation);
            const REAL delta_energy = cell.potential(*particle) - energy_before;

            // rejection
            if(!acceptance.isValid(delta_energy))
            {
                sw_position.rejected();
                // particle->getCoordinates() = particle->getCoordinates() - translation;
                particle->try_setCoordinates(particle->getCoordinates() - translation);
            }
            // acctance
            else
            {
                // vesLOG("after " << particle->getCoordinates().format(ROWFORMAT));
                sw_position.accepted();
            }
        }

        // orientation move
        {
            const REAL stepwidth = sw_orientation();
            const auto random_vec = Particle::Base::cartesian::Random();

            const Eigen::AngleAxisf rotation (stepwidth, random_vec);

            const REAL energy_before = cell.potential(*particle);
            particle->try_setOrientation(rotation * particle->getOrientation());
            const REAL delta_energy = cell.potential(*particle) - energy_before;

            // rejection
            if(!acceptance.isValid(delta_energy))
            {
                sw_orientation.rejected();
                const Eigen::AngleAxisf rotation_back (-stepwidth, random_vec);
                // particle->getOrientation() = (rotation_back * particle->getOrientation());
                particle->try_setOrientation(rotation_back * particle->getOrientation());
            }
            // acctance
            else
                sw_orientation.accepted();
        }
    }
}



REAL ves::MonteCarloSystem::potential() const
{
    // vesDEBUG(__PRETTY_FUNCTION__);
    auto energy = std::make_shared<REAL>(0);
    tbb::parallel_for(std::size_t(0), particles.data.size(), std::size_t(1), [&](const std::size_t& i)
    {
        float pre_sum = 0;
        for(std::size_t j = 0; j<i; ++j)
        {
            pre_sum += interaction->calculate(*particles.data[i], *particles.data[j]);
        }

        // this works for an atomic float, but may busy wait to compare_exchange like hell if too many particles
        *energy += pre_sum;
    });
    return !std::isnan(*energy) ? *energy : throw std::runtime_error("potential Energy is NAN");
}