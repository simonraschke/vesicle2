/*  
*   Copyright 2017-2019 Simon Raschke
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

#include "tube.hpp"



std::ostream& ves::TubeGeometry::print(std::ostream& os) const
{   
    os << "TubeGeometry of height " << height << " with radius " << radius << " and " << points.size() << " points\n";
    for(const auto& point : points)
        os << "point: " << point.format(ROWFORMAT) << '\n';
    return os;
}



ves::TubeGeometry::TubeGeometry()
{

}



ves::TubeGeometry::TubeGeometry(float _height, float _radius, std::size_t _size)
    : height(_height)
    , radius(_radius)
    , size(_size)
{
    generate();
}



void ves::TubeGeometry::generate()
{
    std::size_t best_num_height = size;
    std::size_t best_num_circle = 1;

    for(decltype(size) circle_size = size/4; circle_size >= std::max(size/25, decltype(size)(3)); --circle_size)
    {
        vesLOG("circle size " << circle_size << " height size " << size / circle_size << "   remainder: " << size % circle_size);
        if(size%circle_size == 0)
        {
            best_num_circle = circle_size;
            best_num_height = size / circle_size;
        }
    }
    if(best_num_height == size && best_num_circle == 1)
        vesCRITICAL("unable to determine " << *this);

    vesLOG("TubeGeometry with " << best_num_circle << " points per circle and a height of " << best_num_height << " circles");
    
    for(decltype(best_num_height) i = 0; i < best_num_height; ++i)
    {
        cartesian origin(0, 0, height/(size-1)*i*best_num_circle);
        auto circle = CircleGeometry(origin, radius, best_num_circle);
        circle.rotate(TWOPI/(2*best_num_circle)*i);
        for(auto point : circle.points)
            points.emplace_back(point);
    }
}



void ves::TubeGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void ves::TubeGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}