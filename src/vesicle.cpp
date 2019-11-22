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


#include "common/definitions.hpp"
#include "particles/mobile.hpp"
#include "particles/frame.hpp"
#include "particles/osmotic.hpp"
#include "io/parameters.hpp"
#include "control/controller.hpp"
#include "enhance/output_utility.hpp"
#include <tbb/task_scheduler_init.h>


#include "geometries/helix.hpp"
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
    std::signal( SIGUSR1, ves::Controller::signal );
    std::signal( SIGUSR2, ves::Controller::signal );

    ves::Parameters::getInstance().read(argc,argv);
    vesDEBUG( ves::Parameters::getInstance().getOption("system.time_max").as<std::size_t>() );

    tbb::task_scheduler_init init(ves::Parameters::getInstance().getOption("general.cpu_threads").as<std::size_t>() );
    
    ves::Controller controller;
    controller.setup();

    // run with wall and cpu time measurements
    double wall0 = enhance::get_wall_time();
    double cpu0  = enhance::get_cpu_time();
    controller.start();
    double wall1 = enhance::get_wall_time();
    double cpu1  = enhance::get_cpu_time();

    vesLOG("WALLTIME         " << wall1 - wall0);
    vesLOG("CPUTIME          " << cpu1 - cpu0);
    vesLOG("CPUTIMEPERTHREAD " << (cpu1 - cpu0)/ves::Parameters::getInstance().getOption("general.cpu_threads").as<std::size_t>());

    if(ves::Controller::SIGNAL.load() != 0)
        return ves::Controller::SIGNAL.load();
    else
        return EXIT_SUCCESS;
}