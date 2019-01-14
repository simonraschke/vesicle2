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
#include "particles/particle_container.hpp"
#include "algorithm/cell_container.hpp"
#include "algorithm/metropolis.hpp"
#include "particles/interaction.hpp"
#include "stepwidth_alignment_unit.hpp"
#include "io/parameters.hpp"
#include "io/trajectory_writer_gro.hpp"



namespace ves { struct MonteCarloSystem; }



struct ves::MonteCarloSystem
{
    MonteCarloSystem();

    void setup();
    void run();
    REAL potential() const;
    
    inline auto getTime() const { return time; }
    inline auto getParticles() const { return std::cref(particles); }
    inline auto getBox() const { return std::cref(box); }

protected:
    std::size_t time {0};
    const std::size_t time_max;
    const std::size_t output_skip;
    
    void cellStep(const ves::Cell&);

    ves::CellContainer cells;
    ves::ParticleContainer particles;
    std::unique_ptr<ves::AngularLennardJonesInteraction> interaction {nullptr};

    ves::StepwidhtAlignmentUnit sw_position;
    ves::StepwidhtAlignmentUnit sw_orientation;

    ves::MetropolisAcceptance acceptance;

    ves::Box<PERIODIC::ON> box;
    ves::TrajectoryWriterGro traj_gro;
};



#include "control/controller.hpp"