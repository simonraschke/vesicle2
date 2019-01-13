/*  
*  
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

// #include "systems/system.hpp"
#include "common/definitions.hpp"
#include "monte_carlo_system.hpp"
#include "io/parameters.hpp"
#include <tbb/mutex.h>
#include <csignal>
#include <chrono>
#include <atomic>


namespace ves { struct Controller; }



// a simulation control base class to derive from
// is parameter dependent component
// can catch signals if registered by std::signal 
struct ves::Controller
    // : public ParameterDependentComponent
{
    void setup();
    void start();
    // void pause();
    
    // static member function to catch signal
    // store in atomic which is accesible by derived
    static void signal(int SIG);
    static std::atomic<int> SIGNAL;

protected:

    // the actual system
    ves::MonteCarloSystem system;

    // signal handling
    static tbb::mutex signal_mutex;
};
