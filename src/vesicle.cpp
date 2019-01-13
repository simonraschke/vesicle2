/*  
*   Copyright 2017-2018 Simon Raschke
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


#include "common/definitions.hpp"
#include "particles/mobile.hpp"
#include "particles/frame.hpp"
#include "particles/osmotic.hpp"
#include "io/parameters.hpp"
#include "control/controller.hpp"
#include <tbb/task_scheduler_init.h>



auto main(int argc, const char *argv[]) -> int
{
    // register important signals in Controller base class
    // allowing civilized shutdown
    std::signal( SIGHUP,  ves::Controller::signal );
    std::signal( SIGINT,  ves::Controller::signal );
    std::signal( SIGQUIT, ves::Controller::signal );
    std::signal( SIGILL,  ves::Controller::signal );
    std::signal( SIGTRAP, ves::Controller::signal );
    std::signal( SIGABRT, ves::Controller::signal );
    std::signal( SIGIOT,  ves::Controller::signal );
    std::signal( SIGBUS,  ves::Controller::signal );
    std::signal( SIGFPE,  ves::Controller::signal );
    std::signal( SIGKILL, ves::Controller::signal );
    std::signal( SIGTERM, ves::Controller::signal );

    ves::Parameters prms;
    prms.read(argc,argv);
    vesDEBUG( prms.getOption("system.time_max").as<float>() );

    // // create a task arena 
    // tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
    // tbb::task_arena limited(prms.cpu_threads);
    
    // // execute in limited task arena
    // SimulationControl control;
    // control.setParameters(prms);
    // control.setup();
    // control.start();

    std::vector<std::unique_ptr<ves::Particle::Base>> particles;
    particles.emplace_back(std::make_unique<ves::Particle::Mobile>());
    particles.emplace_back(std::make_unique<ves::Particle::Frame>());
    particles.emplace_back(std::make_unique<ves::Particle::Osmotic>());

    return EXIT_SUCCESS;
}