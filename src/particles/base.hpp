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
#include "boundings.hpp"

#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>



namespace ves
{ 
    namespace Particle
    {
        enum class TYPE : std::int16_t
        {
            ANY = -1,
            UNDEFINED = -1,
            MOBILE = 0,
            FRAME = 1,
            OSMOTIC = 2
        };
        
        struct Base;
        struct Mobile;
        struct Frame;
        struct Osmotic;
    }
}



struct ves::Particle::Base
    // : public ves::Component<Base>
{

    using cartesian = Eigen::Matrix<REAL,3,1>; 
    using box3d = Eigen::AlignedBox<REAL,3>; 

protected:
    Base() = default;

    cartesian coordinates = {0,0,0};
    cartesian orientation = {0,0,0};
    REAL LJ_attraction_intensity = {1};
    REAL LJ_rejection_intensity = {1};

    // nonstd::observer_ptr<Cell> paren_cell {nullptr};

public:
    virtual ~Base() = default;
    // inline virtual Base* getDerived() override { return this; };

    inline bool operator==(Base& other) const { return std::addressof(*this) == std::addressof(other); };
    inline bool operator!=(Base& other) const { return std::addressof(*this) != std::addressof(other); };

    inline virtual auto getType() const -> TYPE = 0;

    auto getCoordinates() -> decltype(coordinates)&;
    auto getCoordinates() const -> const decltype(coordinates)&;

    auto getOrientation() -> decltype(orientation)&;
    auto getOrientation() const -> const decltype(orientation)&;

    auto getLJAttraction() const -> decltype(LJ_attraction_intensity);
    auto getLJRejection() const -> decltype(LJ_rejection_intensity);

    void setLJAttraction(const decltype(LJ_attraction_intensity));
    void setLJRejection(const decltype(LJ_attraction_intensity));

    CoordinatesBounding coordinates_bounding;
    OrientationBounding orientation_bounding;
    
    virtual bool try_setCoordinates(const cartesian&);
    virtual bool try_setOrientation(const cartesian&);

    virtual std::string getName() const = 0;
};