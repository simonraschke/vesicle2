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



#include "base.hpp"



auto ves::Particle::Base::getCoordinates() -> decltype(coordinates)&
{
    return coordinates;
}



auto ves::Particle::Base::getCoordinates() const -> const decltype(coordinates)&
{
    return coordinates;
}



auto ves::Particle::Base::getOrientation() -> decltype(orientation)&
{
    return orientation;
}



auto ves::Particle::Base::getOrientation() const -> const decltype(orientation)&
{
    return orientation;
}



// void ves::Particle::Base::setBoundingBox(box3d&& b)
// {
//     coordinates_bounding.setBoundingBox(std::move(b));
// }



void ves::Particle::Base::setLJAttraction(const decltype(LJ_attraction_intensity) var)
{
    LJ_attraction_intensity = var;
}



void ves::Particle::Base::setLJRejection(const decltype(LJ_attraction_intensity) var)
{
    LJ_rejection_intensity = var;
}



auto ves::Particle::Base::getLJAttraction() const -> decltype(LJ_attraction_intensity)
{
    return LJ_attraction_intensity;
}



auto ves::Particle::Base::getLJRejection() const -> decltype(LJ_attraction_intensity)
{
    return LJ_rejection_intensity;
}



bool ves::Particle::Base::try_setCoordinates(const cartesian& c)
{
    if(coordinates_bounding.isAllowed(c))
    {
        coordinates = c;
        return true;
    }
    else return false;
}



bool ves::Particle::Base::try_setOrientation(const cartesian& o)
{
    if(orientation_bounding.isAllowed(o.normalized()))
    {
        orientation = o.normalized();
        return true;
    }
    else return false;
}



void ves::Particle::Base::forcefullyShift(const cartesian& s)
{
    coordinates += s;
    coordinates_bounding.forcefullyShift(s);
    orientation_bounding.forcefullyShift(s);
}