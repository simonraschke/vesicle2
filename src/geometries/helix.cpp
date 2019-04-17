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

#include "helix.hpp"



std::ostream& ves::HelixGeometry::print(std::ostream& os) const
{   
    os << "HelixGeometry of height " << height << " with radius " << radius << " and " << points.size() << " points\n";
    for(const auto& point : points)
        os << "point: " << point.format(ROWFORMAT) << '\n';
    return os;
}



ves::HelixGeometry::HelixGeometry()
{

}



ves::HelixGeometry::HelixGeometry(float _height, float _radius, std::size_t _size)
    : height(_height)
    , radius(_radius)
    , size(_size)
{
    generate();
}



void ves::HelixGeometry::generate()
{
    points.resize(size);

    float phi, phiold;
    const std::size_t size_old = size++;

    for(decltype(size) i = 0; i < points.size(); ++i)
    {
        // Distributing many points on a sphere
        if ( i == 0 )
        { 
            // first particle fix
            phi = 0.0; 
        }
        // else if ( i == size - 1)
        // { 
        //     // last particle fix
        //     phi = 0.0; 
        // }
        else
        { 
            phiold = phi;
            // claculate the angles in between
            float hk = -1.0 + (float)(2*i)/(size-1);
            phi = phiold + 3.6/(std::sqrt(size)*std::sqrt(1.0-hk*hk));
        }
        
        const Eigen::AngleAxisf rotateZ(phi, cartesian::UnitZ());
        points[i] = (rotateZ * cartesian::UnitY()).normalized() * radius;
        points[i](2) = height/(points.size()-1)*i;
        // points[i](2) = height/points.size()*i;
    }
    size = size_old;
}



void ves::HelixGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void ves::HelixGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}