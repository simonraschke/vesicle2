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
#include "enhance/math_utility.hpp"
#include <cmath>



namespace ves 
{ 
    struct GrandCanonicalInsertion; 
    struct GrandCanonicalDeletion;
    struct MonteCarloSystem;
}



struct ves::GrandCanonicalInsertion
{
    bool isValid(const MonteCarloSystem& sys, REAL delta_U) const;

protected:
    const REAL temperature = Parameters::getInstance().getOption("system.temperature").as<REAL>();
    const REAL mu = Parameters::getInstance().getOption("system.mu").as<REAL>();
};



struct ves::GrandCanonicalDeletion
{
    bool isValid(const MonteCarloSystem& sys, REAL delta_U) const;

protected:
    const REAL temperature = Parameters::getInstance().getOption("system.temperature").as<REAL>();
    const REAL mu = Parameters::getInstance().getOption("system.mu").as<REAL>();
};



#include "control/monte_carlo_system.hpp"