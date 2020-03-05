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
#include "common/component.hpp"
#include "base.hpp"

#include <memory>

#if __has_include(<eigen3/Eigen/Core>)
    #include <eigen3/Eigen/Core>
#elif __has_include(<Eigen/Core>)
    #include <Eigen/Core>
#else
    #pragma error "no eigen include"
#endif

#if __has_include(<eigen3/Eigen/Geometry>)
    #include <eigen3/Eigen/Geometry>
#elif __has_include(<Eigen/Geometry>)
    #include <Eigen/Geometry>
#else
    #pragma error "no eigen include"
#endif




struct ves::Particle::Osmotic
    : public Base
{
    Osmotic() : Base()
    {
        setLJAttraction(0);
        setLJRejection(1);
    }

    inline virtual auto getType() const -> TYPE override { return TYPE::OSMOTIC; }
    // inline virtual Osmotic* getDerived() override { return this; };
    inline virtual std::string getName() const override { return "OSMOT"; };
};