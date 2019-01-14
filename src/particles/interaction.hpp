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
#include "io/parameters.hpp"
#include "particles/base.hpp"
#include "common/box.hpp"
#include "enhance/math_utility.hpp"
#include <memory>



namespace ves { struct AngularLennardJonesInteraction; };



struct ves::AngularLennardJonesInteraction
{
public:
    using particle_t = Particle::Base;
    using box3d = particle_t::box3d;
    using cartesian = particle_t::cartesian;

    AngularLennardJonesInteraction();

    REAL calculate(const particle_t&, const particle_t&) const;


protected:
    ves::Box<PERIODIC::ON> box;

    const REAL kappa;
    const REAL a;
    const REAL b;
    const REAL c;
    const REAL sigma;
    const REAL epsilon;
    const REAL cutoff_rez_sq;
};
