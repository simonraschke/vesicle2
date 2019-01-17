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

#pragma once

#include "io/parameters.hpp"
#include "enhance/random.hpp"
#include <cmath>



namespace ves { struct MetropolisAcceptance; }



struct ves::MetropolisAcceptance
{
    inline bool isValid(float energy_difference) const
    {
        #ifndef NDEBUG
            const REAL exr = std::exp(-energy_difference/temperature);
            const REAL ran = enhance::random<REAL>(0.f,1.f);
            const bool acc = exr > ran;
            // vesDEBUG("energy difference: "<< energy_difference << "  temp: " << temperature << "  exp: " << exr <<"  random: " << ran << "  accepted: " << std::boolalpha << acc )
            assert( energy_difference < 0.f ? acc : true);
            return acc;
        #else
            return energy_difference < 0.f ? true : std::exp(-energy_difference/temperature) > enhance::random<REAL>(0.0,1.0);
        #endif
    }

protected:
    const REAL temperature = Parameters::getInstance().getOption("system.temperature").as<REAL>();
};


