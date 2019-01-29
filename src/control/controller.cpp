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

#include "controller.hpp"



std::atomic<int> ves::Controller::SIGNAL = {0};
tbb::mutex ves::Controller::signal_mutex {};



void ves::Controller::signal(int SIG)
{
    static std::size_t got_called = 0;
    ++got_called;

    if(SIG == SIGABRT)
    {
        std::exit(SIG);
    }
    
    if(got_called == 1)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
    }
    else if (got_called == 2)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
        vesWARNING(__PRETTY_FUNCTION__ << " still trying civilized shutdown...")
    }
    else if (got_called == 3)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieverd SIGNAL " << SIG)
        vesWARNING(__PRETTY_FUNCTION__ << " TERMINATING!")
        std::exit(SIG);
    }

    tbb::mutex::scoped_lock lock(Controller::signal_mutex);
    SIGNAL.store(SIG, std::memory_order_relaxed);
}



void ves::Controller::setup()
{
    vesDEBUG(__PRETTY_FUNCTION__);

    system.setup();
    
    vesLOG("");
    vesLOG("SIZE OF DATA STRUCTURES");
    vesLOG("int                             " << sizeof(int));
    vesLOG("float                           " << sizeof(float));
    vesLOG("REAL                            " << sizeof(REAL));
    vesLOG("ves::Particle::Base::cartesian  " << sizeof(ves::Particle::Base::cartesian));
    vesLOG("ves::Particle::Base             " << sizeof(ves::Particle::Base));
    vesLOG("ves::Particle::Mobile           " << sizeof(ves::Particle::Mobile));
    vesLOG("ves::Particle::Frame            " << sizeof(ves::Particle::Frame));
    vesLOG("ves::Particle::Osmotic          " << sizeof(ves::Particle::Osmotic));
    vesLOG("ves::ParticleContainer          " << sizeof(ves::ParticleContainer));
    vesLOG("ves::Cell                       " << sizeof(ves::Cell));
    vesLOG("ves::CellContainer              " << sizeof(ves::CellContainer));
    vesLOG("ves::Box<PERIODIC::ON>          " << sizeof(ves::Box<PERIODIC::ON>));
    vesLOG("ves::Box<PERIODIC::OFF>         " << sizeof(ves::Box<PERIODIC::OFF>));
    vesLOG("");

}


void ves::Controller::start()
{
    vesDEBUG(__PRETTY_FUNCTION__);
    system.run();
}