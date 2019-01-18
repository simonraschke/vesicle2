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
#include "algorithm/grand_canonical.hpp"
#include "particles/interaction.hpp"
#include "stepwidth_alignment_unit.hpp"
#include "io/parameters.hpp"
#include "io/trajectory_writer_gro.hpp"
#include "io/trajectory_writer_h5.hpp"
#include <chrono>
#include <random>



namespace ves { struct MonteCarloSystem; }



struct ves::MonteCarloSystem
{
    MonteCarloSystem();

    void setup();
    void run();
    REAL potential() const;

    auto status();
    
    inline auto getTime() const { return time; }
    inline auto getParticles() const { return std::cref(particles); }
    inline auto getBox() const { return std::cref(box); }
    inline auto getCells() const { return std::cref(cells); }
    inline auto getInteraction() const { return std::cref(*interaction); }
    inline auto getSWPosition() const { return std::cref(sw_position); }
    inline auto getSWOrientation() const { return std::cref(sw_orientation); }

protected:
    std::size_t time {0};
    const std::size_t time_max;
    const std::size_t output_skip;
    
    void cellStep(const ves::Cell&);
    void grandCanonicalStep();
    bool try_addParticle();
    bool try_removeParticle();

    ves::CellContainer cells;
    ves::ParticleContainer particles;
    std::unique_ptr<ves::AngularLennardJonesInteraction> interaction {nullptr};
    ves::StepwidhtAlignmentUnit sw_position;
    ves::StepwidhtAlignmentUnit sw_orientation;
    ves::MetropolisAcceptance acceptance;
    ves::Box<PERIODIC::ON> box;
    ves::TrajectoryWriterGro traj_gro;
    ves::TrajectoryWriterH5 traj_h5;

    thread_local static std::mt19937_64 pseudo_engine;
};



#include "control/controller.hpp"