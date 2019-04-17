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

#include "grid.hpp"



ves::GridGeometry::GridGeometry()
{

}



ves::GridGeometry::GridGeometry(std::size_t _x, std::size_t _y, std::size_t _z)
    : x(_x)
    , y(_y)
    , z(_z)
{
    generate();
}


void ves::GridGeometry::generate()
{
    points.reserve(x*y*z);

    // tbb::parallel_for( tbb::blocked_range3d<std::size_t>(0, x, 0, y, 0, z),
    // [&]( const tbb::blocked_range3d<std::size_t> &r ) {
    //     for(std::size_t i=r.pages().begin(), i_end=r.pages().end(); i<i_end; i++){
    //         for(std::size_t j=r.rows().begin(), j_end=r.rows().end(); j<j_end; j++){
    //             for(std::size_t k=r.cols().begin(), k_end=r.cols().end(); k<k_end; k++){
    //                 #ifdef __clang_major__
    //                 points.push_back(cartesian(i,j,k));
    //                 #elif  __GNUC__
    //                 points.emplace_back(cartesian(i,j,k));
    //                 #else
    //                     #error no valid compiler
    //                 #endif
    //             }
    //         }
    //     }
    // });

    for(std::size_t i = 0; i<x; ++i)
        for(std::size_t j = 0; j<y; ++j)
            for(std::size_t k = 0; k<z; ++k)    
            {
                points.emplace_back(cartesian(i,j,k));
            }
}



ves::Geometry::cartesian& ves::GridGeometry::matrixView(const std::size_t& i, const std::size_t& j, const std::size_t& k)
{
    return points[ i*x*y + j*y + k ];
}



const ves::Geometry::cartesian& ves::GridGeometry::matrixView(const std::size_t& i, const std::size_t& j, const std::size_t& k) const
{
    return points[ i*x*y + j*y + k ];
}



void ves::GridGeometry::scale(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) *= c(0);
        point(1) *= c(1);
        point(2) *= c(2);
    }
}



void ves::GridGeometry::shift(const cartesian& c)
{
    for(auto& point : points)
    {
        point(0) += c(0);
        point(1) += c(1);
        point(2) += c(2);
    }
}