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

#include "circle.hpp"



std::ostream& ves::CircleGeometry::print(std::ostream& os) const
{   
    os << "CircleGeometry at " << origin.format(ROWFORMAT) << " with radius " << radius << " and " << points.size() << " points\n";
    for(const auto& point : points)
        os << "point: " << point.format(ROWFORMAT) << '\n';
    return os;
}



ves::CircleGeometry::CircleGeometry()
{

}



ves::CircleGeometry::CircleGeometry(cartesian _origin, float _radius, std::size_t _size)
    : origin(_origin)
    , radius(_radius)
    , size(_size)
{
    generate();
}



void ves::CircleGeometry::generate()
{
    points.resize(size);
    const float angle_per_point = TWOPI / points.size();

    for(decltype(size) i = 0; i < points.size(); ++i)
    {
        const Eigen::AngleAxis<float> rotation(angle_per_point*i, cartesian::UnitZ());
        points[i] = (rotation * cartesian::UnitX())*radius;
    }

    shift(origin);
}



void ves::CircleGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void ves::CircleGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}



void ves::CircleGeometry::rotate(float angle)
{
    shift(-origin);
    for(auto& point : points)
    {
        point = Eigen::AngleAxis<float>(angle, cartesian::UnitZ())*point;
    }
    shift(origin);
}