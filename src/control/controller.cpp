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
        vesWARNING(__PRETTY_FUNCTION__ << " recieved SIGNAL " << SIG)
    }
    else if (got_called == 2)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieved SIGNAL " << SIG)
        vesWARNING(__PRETTY_FUNCTION__ << " still trying civilized shutdown...")
    }
    else if (got_called == 3)
    {
        vesWARNING(__PRETTY_FUNCTION__ << " recieved SIGNAL " << SIG)
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
    vesLOG(std::setw(30) << std::left << "int" << std::setw(4) << std::right << sizeof(int) << std::setw(4) << std::right << alignof(int));
    vesLOG(std::setw(30) << std::left << "float" << std::setw(4) << std::right << sizeof(float) << std::setw(4) << std::right << alignof(float));
    vesLOG(std::setw(30) << std::left << "REAL" << std::setw(4) << std::right << sizeof(REAL) << std::setw(4) << std::right << alignof(REAL));
    vesLOG(std::setw(30) << std::left << "ves::Particle::Base::cartesian" << std::setw(4) << std::right << sizeof(ves::Particle::Base::cartesian) << std::setw(4) << std::right << alignof(ves::Particle::Base::cartesian));
    vesLOG(std::setw(30) << std::left << "ves::CoordinatesBounding" << std::setw(4) << std::right << sizeof(ves::CoordinatesBounding) << std::setw(4) << std::right << alignof(ves::CoordinatesBounding));
    vesLOG(std::setw(30) << std::left << "ves::OrientationBounding" << std::setw(4) << std::right << sizeof(ves::OrientationBounding) << std::setw(4) << std::right << alignof(ves::OrientationBounding));
    vesLOG(std::setw(30) << std::left << "ves::Particle::Base" << std::setw(4) << std::right << sizeof(ves::Particle::Base) << std::setw(4) << std::right << alignof(ves::Particle::Base));
    vesLOG(std::setw(30) << std::left << "ves::Particle::Mobile" << std::setw(4) << std::right << sizeof(ves::Particle::Mobile) << std::setw(4) << std::right << alignof(ves::Particle::Mobile));
    vesLOG(std::setw(30) << std::left << "ves::Particle::Frame" << std::setw(4) << std::right << sizeof(ves::Particle::Frame) << std::setw(4) << std::right << alignof(ves::Particle::Frame));
    vesLOG(std::setw(30) << std::left << "ves::Particle::Osmotic" << std::setw(4) << std::right << sizeof(ves::Particle::Osmotic) << std::setw(4) << std::right << alignof(ves::Particle::Osmotic));
    vesLOG(std::setw(30) << std::left << "ves::ParticleContainer" << std::setw(4) << std::right << sizeof(ves::ParticleContainer) << std::setw(4) << std::right << alignof(ves::ParticleContainer));
    vesLOG(std::setw(30) << std::left << "ves::Cell" << std::setw(4) << std::right << sizeof(ves::Cell) << std::setw(4) << std::right << alignof(ves::Cell));
    vesLOG(std::setw(30) << std::left << "ves::CellContainer" << std::setw(4) << std::right << sizeof(ves::CellContainer) << std::setw(4) << std::right << alignof(ves::CellContainer));
    vesLOG(std::setw(30) << std::left << "ves::Box<PERIODIC::ON>" << std::setw(4) << std::right << sizeof(ves::Box<PERIODIC::ON>) << std::setw(4) << std::right << alignof(ves::Box<PERIODIC::ON>));
    vesLOG(std::setw(30) << std::left << "ves::Box<PERIODIC::OFF>" << std::setw(4) << std::right << sizeof(ves::Box<PERIODIC::OFF>) << std::setw(4) << std::right << alignof(ves::Box<PERIODIC::OFF>));
    vesLOG("");

}


void ves::Controller::start()
{
    vesDEBUG(__PRETTY_FUNCTION__);
    system.run();
}